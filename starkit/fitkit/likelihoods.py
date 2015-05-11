from astropy import modeling
import numpy as np

class Chi2Likelihood(modeling.Model):
    inputs = ('wavelength', 'flux')
    outputs = ('loglikelihood')

    def __init__(self, observed):
        self.observed_wavelength = observed.wavelength.to(u.angstrom).value
        self.observed_flux = observed.flux.value
        self.observed_uncertainty = getattr(observed.uncertainty, None)
        if self.observed_uncertainty is not None:
            self.observed_uncertainty = self.observed_uncertainty.value


    def evaluate(self, wavelength, flux):
        return -0.5 * np.sum(
            ((self.observed_flux - flux) / self.observed_uncertainty)**2)

