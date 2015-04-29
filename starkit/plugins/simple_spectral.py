import numpy as np
from astropy import constants as const
from astropy.modeling import Parameter
from specutils import Spectrum1D

from starkit.plugins.base import StarkitModel

c_kms = const.c.to('km/s').value


class DopplerShift(StarkitModel):
    inputs = ('spectrum', )
    outputs = ('spectrum', )

    v_rad = Parameter(default=0.0)

    def evaluate(spectrum, v_rad):
        beta = v_rad / c_kms
        doppler_factor = np.sqrt((1+beta) / (1-beta))
        Spectrum1D.from_array(spectrum.wavelength * doppler_factor,
                                     spectrum.flux)



class WavelengthSlice(StarkitModel):
    inputs = ('spectrum', )
    outputs = ('spectrum', )

    slice_min_wavelength = Parameter(default=-np.inf)
    slice_max_wavelength = Parameter(default=np.inf)

    def evaluate(self, wavelength, flux, uncertainty):
        if not hasattr(self, 'wavelength'):
            self.wavelength = spectrum.wavelength
            self.slice_min_idx = spectrum.wavelength.value.searchsorted(
            self.slice_min_wavelength)

            self.slice_max_idx = spectrum.wavelength.value.searchsorted(
            self.slice_max_wavelength)



        return Spectrum1D.from_array(
            spectrum.wavelength[self.slice_min_idx:self.slice_max_idx],
            spectrum.flux[self.slice_min_idx:self.slice_max_idx])





