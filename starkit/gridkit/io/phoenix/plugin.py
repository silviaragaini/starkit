from astropy.modeling import FittableModel, Parameter
import numpy as np
from scipy import ndimage as nd

class PhoenixResolution(object):

    uv_wavelength = (500, 3000)
    uv_R = 0.1/3000
    oir_wavelength = (3000, 25000)
    oir_R = 5e5
    mir_wavelength = (25000, 55000)
    mir_R = 1e5

    def __init__(self, wavelength, R, wavelength_range = (0, np.inf)):

        self.wavelength = wavelength

        self.R = R
        self.wavelength_range = wavelength_range


        uv_wavelength_interp = np.arange(np.log(self.uv_wavelength[0]),
                                             np.log(self.uv_wavelength),
                                             1 /(2 * self.uv_R))
        oir_wavelength_interp = np.arange(np.log(self.oir_wavelength[0]),
                                             np.log(self.oir_wavelength[1]),
                                             1 / 2 * self.oir_R)
        mir_wavelength_interp = np.arange(np.log(self.oir_wavelength[0]),
                                             np.log(self.oir_wavelength[1]),
                                             1 / 2 * self.oir_R)

        self.wavelength_interp = np.hstack((uv_wavelength_interp,
                                            oir_wavelength_interp,
                                            mir_wavelength_interp))


    def __call__(self, flux):

        wavelength, flux = spectrum.wavelength.value, spectrum.flux
        log_grid_log_wavelength = np.arange(np.log(wavelength.min()),
                                            np.log(wavelength.max()),
                                            1 / (self.sampling *
                                                 self.R.to(1).value))
        log_grid_wavelength = np.exp(log_grid_log_wavelength)
        log_grid_flux = np.interp(log_grid_wavelength, wavelength, flux)
        sigma = self.sampling / (2 * np.sqrt(2 * np.log(2)))
        log_grid_convolved = nd.gaussian_filter1d(log_grid_flux, sigma)
        convolved_flux = np.interp(wavelength, log_grid_wavelength,
                                   log_grid_convolved)

