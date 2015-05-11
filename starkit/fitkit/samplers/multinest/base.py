from astropy.modeling import fitting
import pymultinest

class MultiNest(fitting.Fitter):
    def __init__(self):
        super(MultiNest, self).__init__()

    def __call__(self, model, spectrum):
        pass
