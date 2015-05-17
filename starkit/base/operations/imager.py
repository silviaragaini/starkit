import numpy as np

from astropy import units as u
from starkit.fix_spectrum1d import Spectrum1D
from starkit.base.operations.base import InstrumentOperationModel

class ImagerInstrumentOperation(InstrumentOperationModel):
    pass

class Photometry(ImagerInstrumentOperation):
    inputs = ('wavelength', 'flux')
    outputs = ('photometry',)

    def __init__(self, filter_set, mag_type='vega'):
        super(Photometry, self).__init__()
        try:
            from wsynphot import FilterSet
        except ImportError:
            raise ImportError('The photometry plugin needs wsynphot')

        if hasattr(filter_set, 'calculate_{0}_magnitudes'.format(mag_type)):
            self.filter_set = filter_set
        else:
            self.filter_set = FilterSet(filter_set)

        self.calculate_magnitudes = getattr(
            self.filter_set, 'calculate_{0}_magnitudes'.format(mag_type))


    def evaluate(self, wavelength, flux):
        spec = Spectrum1D.from_array(wavelength * u.angstrom,
                                     flux * u.erg/u.s/u.cm**2/u.angstrom)
        return np.array(u.Quantity(self.calculate_magnitudes(spec)).value)


