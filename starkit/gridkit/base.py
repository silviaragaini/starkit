from astropy import modeling
from astropy import units as u, constants as const
import pandas as pd
import h5py
from scipy import interpolate
from starkit.fitkit.samplers.priors import UniformPrior

import numpy as np



class BaseSpectralGrid(modeling.FittableModel):
    inputs = tuple()
    outputs = ('wavelength', 'flux')

    def __init__(self, wavelength, index, fluxes, **kwargs):
        self.R = kwargs.pop('R', None)
        self.R_sampling = kwargs.pop('R_sampling', None)
        self.flux_unit = kwargs.pop('flux_unit', None)

        super(BaseSpectralGrid, self).__init__(**kwargs)
        self.interpolator = self._generate_interpolator(index, fluxes)
        self.wavelength = wavelength


#    def prepare_inputs(self, *inputs, **kwargs):
#        return inputs, None

#    def prepare_outputs(self, format_info, *outputs, **kwargs):
#        return outputs

    def get_grid_extent(self):
        extents = []
        for i, param_name in enumerate(self.param_names):
            extents.append((self.interpolator.points[:,i].min(),
                            self.interpolator.points[:,i].max()))
        return extents

    def get_grid_uniform_priors(self):
        extents = self.get_grid_extent()
        priors = []
        for extent in extents:
            priors.append(UniformPrior(*extent))

        return priors

    def evaluate(self, *args):

        return self.wavelength, self.interpolator(np.array(
            args).reshape(len(self.param_names)))[0]

    @staticmethod
    def _generate_interpolator(index, fluxes):
        return interpolate.LinearNDInterpolator(index, fluxes)

    @property
    def velocity_per_pix(self):
        if self.R is None:
            raise ValueError('R and R_sampling for the current grid not known')

        else:
            return const.c / self.R / self.R_sampling

def _get_interpolate_parameters(index):
    interpolate_parameters = []

    for column in index.columns:
        if len(index[column].unique()) > 1:
            interpolate_parameters.append(column)
    return interpolate_parameters


def load_grid(hdf_fname):
    index = pd.read_hdf(hdf_fname, 'index')
    interpolate_parameters = _get_interpolate_parameters(index)

    with h5py.File(hdf_fname) as fh:
        fluxes = fh['fluxes'].__array__()
        flux_unit = u.Unit(fh['fluxes'].attrs['unit'])
        wavelength = fh['wavelength'].__array__()
        data_set_type = fh['wavelength'].attrs['grid']
        if data_set_type == 'log':
            R = fh['wavelength'].attrs.get('R', None)
            R_sampling = fh['wavelength'].attrs.get('R_sampling', 4)

        else:
            R = None
            R_sampling = None

        wavelength_unit = u.Unit(fh['wavelength'].attrs['unit'])

    class_dict = {item:modeling.Parameter() for item in interpolate_parameters}
    class_dict['__init__'] = BaseSpectralGrid.__init__

    SpectralGrid = type('SpectralGrid', (BaseSpectralGrid, ), class_dict)

    initial_parameters = {item:index[item].iloc[0] for item in interpolate_parameters}

    return SpectralGrid(wavelength, index[interpolate_parameters], fluxes,
                        R=R, R_sampling=R_sampling, flux_unit=flux_unit,
                        **initial_parameters)








