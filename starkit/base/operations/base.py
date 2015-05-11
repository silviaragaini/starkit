from astropy import modeling
from collections import OrderedDict


class SpectralOperationModel(modeling.FittableModel):

    inputs = ('wavelength', 'flux')
    outputs = ('wavelength', 'flux')


from instrument import InstrumentOperationModel
from stellar import StellarOperationModel
stellar_operations = OrderedDict([(item.operation_name, item)
                         for item in StellarOperationModel.__subclasses__()])
instrument_operations = OrderedDict([(item.operation_name, item)
                         for item in InstrumentOperationModel.__subclasses__()])