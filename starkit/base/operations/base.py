from astropy import modeling
from collections import OrderedDict
from itertools import chain


class SpectralOperationModel(modeling.FittableModel):

    inputs = ('wavelength', 'flux')
    outputs = ('wavelength', 'flux')


from instrument import InstrumentOperationModel
from stellar import StellarOperationModel

stellar_operations = OrderedDict([(item.operation_name, item)
                         for item in StellarOperationModel.__subclasses__()])
stellar_parameter2model = OrderedDict(chain(*[[(name, model) for name in model.param_names]
                           for model in stellar_operations.values()]))

instrument_operations = OrderedDict([(item.operation_name, item)
                         for item in InstrumentOperationModel.__subclasses__()])

instrument_parameter2model = OrderedDict(chain(*[[(name, model)
                                      for name in model.param_names]
                                     for model in instrument_operations.values()]))
