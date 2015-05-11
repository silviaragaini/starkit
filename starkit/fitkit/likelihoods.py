from astropy import modeling

class Chi2Likelihood(modeling.Model):
    inputs = ('wavelength', 'flux')
    outputs = ('loglikelihood')
