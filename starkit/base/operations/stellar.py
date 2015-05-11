from astropy import constants as const, units as u
from specutils import Spectrum1D
from astropy import modeling
from scipy.signal import fftconvolve
import numpy as np

from starkit.base.operations.base import SpectralOperationModel

class StellarOperationModel(SpectralOperationModel):
    pass

class RotationalBroadening(StellarOperationModel):
    operation_name = 'rotation'
    vrot = modeling.Parameter()
    limb_darkening = modeling.Parameter(fixed=True, default=0.6)

    def __init__(self, velocity_per_pix=None, vrot=0):
        super(RotationalBroadening, self).__init__(vrot=vrot)

        self.c_in_kms = const.c.to(u.km / u.s).value

        if velocity_per_pix is not None:
            self.log_sampling = True
            self.velocity_per_pix = u.Quantity(velocity_per_pix, u.km / u.s).value
        else:
            self.log_sampling = False

    def rotational_profile(self, vrot, limb_darkening):
        vrot_by_c = np.maximum(0.0001, np.abs(vrot)) / self.c_in_kms

        half_width_pix = np.round((self.c_in_kms /
                                   self.velocity_per_pix)).astype(int)
        profile_velocity = (np.linspace(-half_width_pix, half_width_pix,
                                       2 * half_width_pix + 1)
                            * self.velocity_per_pix)
        profile = np.maximum(0.,
                             1. - (profile_velocity / vrot_by_c) ** 2)

        profile = ((2 * (1-limb_darkening) * np.sqrt(profile) +
                    0.5 * np.pi * limb_darkening * profile) /
                   (np.pi * vrot_by_c * (1. - limb_darkening / 3.)))
        return profile/profile.sum()


    def evaluate(self, wavelength, flux, v_rot, limb_darkening):
        if self.velocity_per_pix is None:
            raise NotImplementedError('Regridding not implemented yet')

        if np.abs(v_rot) < 1e-5:
            return wavelength, flux

        profile = self.rotational_profile(v_rot, limb_darkening)

        return wavelength, fftconvolve(flux, profile, mode='same')

class DopplerShift(StellarOperationModel):

    operation_name = 'doppler'

    vrad = modeling.Parameter()

    def __init__(self, vrad=0):
        super(DopplerShift, self).__init__(vrad=0)
        self.c_in_kms = const.c.to(u.km / u.s).value


    def evaluate(self, wavelength, flux, vrad):

        beta = self.vrad / self.c_in_kms
        doppler_factor = np.sqrt((1+beta) / (1-beta))
        return wavelength * doppler_factor, flux



class CCM89Extinction(StellarOperationModel):

    operation_name = 'ccm89_extinction'

    a_v = modeling.Parameter(default=0.0)
    r_v = modeling.Parameter(default=3.1)

    @property
    def ebv(self):
        return self.a_v / self.r_v

    def __init__(self, a_v=0.0, r_v=3.1):
        super(CCM89Extinction, self).__init__(a_v=a_v, r_v=r_v)


    def evaluate(self, wavelength, spectrum, a_v, r_v):
        from specutils import extinction
        extinction_factor = np.ones_like(wavelength)
        valid_wavelength = ((wavelength > 910 * u.angstrom) &
                            (wavelength < 33333 * u.angstrom))

        extinction_factor[valid_wavelength] = 10 ** (-0.4 * extinction.extinction_ccm89(
            spectrum.wavelength[valid_wavelength], a_v=self.a_v,
            r_v=self.r_v).to(u.angstrom).value)


        return Spectrum1D.from_array(spectrum.wavelength,
                                     extinction_factor * spectrum.flux)