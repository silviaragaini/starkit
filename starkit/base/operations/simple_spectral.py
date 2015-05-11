from astropy import constants as const, units as u
from specutils import Spectrum1D
from astropy import modeling
from scipy.signal import fftconvolve
import numpy as np

from starkit.base.operations.base import SpectralOperationModel

class RotationalBroadening(SpectralOperationModel):

    vrot = modeling.Parameter()
    limb_darkening = modeling.Parameter(fixed=True, default=0.6)

    resolution = (20 * u.km / u.s / const.c).to(1)
    limb_darkening = 0.6
    param_names = ['vrot']

    def __init__(self, velocity_per_pix=None, **kwargs):
        super(RotationalBroadening, self).__init__(**kwargs)

        self.c_in_kms = const.c.to(u.km / u.s).value

        if velocity_per_pix is not None:
            self.log_sampling = True
            self.velocity_per_pix = u.Quantity(velocity_per_pix, u.km / u.s).value
        else:
            self.log_sampling = False

    def rotational_profile(self, vrot):
        vrot_by_c = np.maximum(0.0001, np.abs(vrot)) / self.c_in_kms

        half_width_pix = np.round((self.c_in_kms /
                                   self.velocity_per_pix)).astype(int)
        profile_velocity = (np.linspace(-half_width_pix, half_width_pix,
                                       2 * half_width_pix + 1)
                            * self.velocity_per_pix)
        profile = np.maximum(0.,
                             1. - (profile_velocity / vrot_by_c) ** 2)

        profile = ((2 * (1-self.limb_darkening) * np.sqrt(profile) +
                    0.5 * np.pi * self.limb_darkening * profile) /
                   (np.pi * vrot_by_c * (1. - self.limb_darkening / 3.)))
        return profile/profile.sum()


    def evaluate(self, wavelength, flux, v_rot):
        if self.velocity_per_pix is None:
            raise NotImplementedError('Regridding not implemented yet')

        profile = self.rotational_profile(v_rot)

        return wavelength, fftconvolve(flux, profile, mode='same')

class DopplerShift(SpectralOperationModel):

    vrad = modeling.Parameter()

    def __init__(self, *args, **kwargs):
        super(DopplerShift, self).__init__(*args, **kwargs)
        self.c_in_kms = const.c.to(u.km / u.s).value


    def evaluate(self, wavelength, flux, vrad):

        beta = self.vrad / self.c_in_kms
        doppler_factor = np.sqrt((1+beta) / (1-beta))
        return wavelength * doppler_factor, flux



class CCM89Extinction(object):
    param_names = ['a_v', 'r_v']


    @property
    def a_v(self):
        return self._a_v

    @a_v.setter
    def a_v(self, value):
        self._a_v = np.abs(value)

    @property
    def r_v(self):
        return self._r_v

    @r_v.setter
    def r_v(self, value):
        self._r_v = np.abs(value)


    def __init__(self, a_v=0.0, r_v=3.1):
        self.a_v = a_v
        self.r_v = r_v

    def __call__(self, spectrum):

        from specutils import extinction
        extinction_factor = np.ones_like(spectrum.wavelength.value)
        valid_wavelength = ((spectrum.wavelength > 910 * u.angstrom) &
                            (spectrum.wavelength < 33333 * u.angstrom))
        extinction_factor[valid_wavelength] = 10 ** (-0.4 * extinction.extinction_ccm89(
            spectrum.wavelength[valid_wavelength], a_v=self.a_v,
            r_v=self.r_v).to(u.angstrom).value)


        return Spectrum1D.from_array(spectrum.wavelength,
                                     extinction_factor * spectrum.flux)