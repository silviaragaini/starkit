from itertools import chain
from astropy.modeling import FittableModel


class StarkitModel(FittableModel):
    def __call__(self, *inputs, **kwargs):

        parameters = self._param_sets(raw=True)
        return self.evaluate(*chain(inputs, parameters))