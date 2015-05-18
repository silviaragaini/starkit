from astropy import units as u, constants as const
from astropy import modeling
import numpy as np

class Chi2Likelihood(modeling.Model):
    inputs = ('wavelength', 'flux')
    outputs = ('loglikelihood', )

    def __init__(self, observed):
        super(Chi2Likelihood, self).__init__()
        self.observed_wavelength = observed.wavelength.to(u.angstrom).value
        self.observed_flux = observed.flux.value
        self.observed_uncertainty = getattr(observed, 'uncertainty', None)
        if self.observed_uncertainty is not None:
            self.observed_uncertainty = self.observed_uncertainty.value
        else:
            self.observed_uncertainty = np.ones_like(self.observed_wavelength)


    def evaluate(self, wavelength, flux):
        loglikelihood =  -0.5 * np.sum(
            ((self.observed_flux - flux) / self.observed_uncertainty)**2)
        return loglikelihood

class PhotometryColorLikelihood(modeling.Model):
    inputs = ('photometry',)
    outputs = ('loglikelihood',)

    def __init__(self, magnitude_set):
        super(PhotometryColorLikelihood, self).__init__()
        self.colors = (magnitude_set.magnitudes[:-1] -
                       magnitude_set.magnitudes[1:])
        self.color_uncertainties = np.sqrt(magnitude_set.magnitudes[:-1]**2 +
                                           magnitude_set.magnitudes[1:]**2)

    def evaluate(self, photometry):
        synth_colors = photometry[:-1] - photometry[1:]
        loglikelihood = -0.5 * np.sum(((self.colors - synth_colors)
                                       / self.color_uncertainties)**2)
        return loglikelihood

class Addition(modeling.Model):

    inputs = ('a', 'b')
    outputs = ('x', )

    def __init__(self):
        super(Addition, self).__init__()

    @staticmethod
    def evaluate(a, b):
        return a + b
