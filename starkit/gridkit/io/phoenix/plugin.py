
import numpy as np
from scipy import ndimage as nd
from scipy.interpolate import interp1d
from specutils import Spectrum1D
from astropy import units as u


class PhoenixProcess(object):

    uv_wavelength = (500, 3000)
    uv_R = 3000 / 0.1
    oir_wavelength = (3000, 25000)
    oir_R = 5e5
    nir_wavelength = (25000, 55000)
    nir_R = 1e5

    def __init__(self, input_wavelength, R, wavelength_range = (0, np.inf),
                 pre_sampling=2, sampling=4):


        self.wavelength_start = input_wavelength.searchsorted(
            wavelength_range[0])
        self.wavelength_end = input_wavelength.searchsorted(wavelength_range[1])

        self.cut_wavelength = input_wavelength[
                              self.wavelength_start:self.wavelength_end]

        self.initial_output_wavelength = self._logrange(
            input_wavelength[self.wavelength_start],
            input_wavelength[self.wavelength_end - 1], R, sampling)
        self.output_wavelength = None

        self.R = R
        self.wavelength_range = wavelength_range
        self.sampling = sampling
        self.pre_sampling = pre_sampling
        self.wavelength_interp = []
        for wl_region in ['uv', 'oir', 'nir']:
            current_wl = getattr(self, '{0}_wavelength'.format(wl_region))
            current_R = getattr(self, '{0}_R'.format(wl_region))
            wavelength_interp = self._logrange(current_wl[0], current_wl[1],
                                               current_R, pre_sampling)
            setattr(self, '{0}_wavelength_interp'.format(wl_region),
                    wavelength_interp)
            self.wavelength_interp.append(wavelength_interp)

        self.wavelength_interp = np.hstack(tuple(self.wavelength_interp))


    @property
    def R(self):
        return self._R

    @R.setter
    def R(self, value):
        self._R = float(value)

    @staticmethod
    def _logrange(wavelength_start, wavelength_end, R, sampling):
        return np.exp(
            np.arange(np.log(wavelength_start), np.log(wavelength_end),
                      1 / (sampling * float(R))))

    def interp_wavelength(self, flux):
        cut_flux = flux[self.wavelength_start:self.wavelength_end]
        processed_wavelengths = []
        processed_fluxes = []
        for wl_region in ['uv', 'oir', 'nir']:
            current_R = getattr(self, '{0}_R'.format(wl_region))
            current_wl_interp = getattr(
                self, '{0}_wavelength_interp'.format(wl_region))
            rescaled_R = 1 / np.sqrt((1/self.R)**2 - (1/current_R)**2 )
            sigma = ((current_R / rescaled_R) * self.sampling /
                     (2 * np.sqrt(2 * np.log(2))))

            interp_flux = interp1d(self.cut_wavelength, cut_flux,
                                   bounds_error=False)(
                current_wl_interp
            )
            nan_filter = ~np.isnan(interp_flux)

            if np.all(~nan_filter):
                processed_wavelengths.append(np.array([]))
                processed_fluxes.append(np.array([]))
            else:
                processed_wavelengths.append(current_wl_interp[nan_filter])
                processed_fluxes.append(nd.gaussian_filter1d(
                    interp_flux[nan_filter], sigma))


        processed_wavelength = np.hstack(tuple(processed_wavelengths))
        processed_flux = np.hstack(tuple(processed_fluxes))

        if self.output_wavelength is None:
            initial_output_flux = interp1d(processed_wavelength, processed_flux,
                     bounds_error=False)(self.initial_output_wavelength)
            self.output_wavelength = self.initial_output_wavelength[~np.isnan(
                initial_output_flux)]
            output_flux = initial_output_flux[~np.isnan(initial_output_flux)]
        else:
            output_flux = interp1d(processed_wavelength, processed_flux)(
                self.output_wavelength
            )
        return output_flux


    def __call__(self, flux):
        new_flux =  self.interp_wavelength(flux)
        return new_flux