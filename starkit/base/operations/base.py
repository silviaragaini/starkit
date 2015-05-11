from astropy import modeling


class SpectralOperationModel(modeling.FittableModel):

    inputs = ('wavelength', 'flux')
    outputs = ('wavelength', 'flux')