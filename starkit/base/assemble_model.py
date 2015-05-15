from starkit.base.operations import (stellar_operations, instrument_operations,
                                     stellar_parameter2model,
                                     instrument_parameter2model)

from starkit.base.operations.instrument import Interpolate, Normalize


def fit_parameters_property(self):
    return [getattr(self, param_name) for param_name in self.param_names
            if getattr(self, param_name).fixed is False]


def assemble_model(spectral_grid, spectrum=None,
                   normalize_npol=None, **kwargs):
    """

    Parameters
    ----------

    spectral_grid: ~specgrid.SpectralGrid
        spectral grid to be used in observation
    spectrum: ~specutils.Spectrum1D
        spectrum to be used for interpolation, if None neither interpolation nor
            will be performed [default None]

    normalize_npol: int
        degree of polynomial to be used for interpolation, only if not None and
        spectrum is not None will the normalization plugin be used [default None]

    plugin_names: ~list of ~str
        select between the following available plugin choices:
        {stellar_operations}
        {instrument_operations}

    Returns
    -------
        : ~Model
        Model with the requested operations


    """

    ObservationModel = spectral_grid
    parameters = kwargs.copy()

    for operation in (stellar_operations.values() +
                          instrument_operations.values()):
        param_values = {}
        for param_name in operation.param_names:
            if param_name in parameters:
                param_values[param_name] = parameters.pop(param_name)

        if param_values != {}:
            if getattr(operation, 'requires_observed', False):
                param_values['observed'] = spectrum

            if hasattr(operation, 'from_grid'):
                current_stellar_operation = operation.from_grid(spectral_grid,
                                                                **param_values)
            else:
                current_stellar_operation = operation(**param_values)
            ObservationModel = ObservationModel | current_stellar_operation

    if parameters != {}:
        raise ValueError('Given parameters {0} not understood'.format(
            ','.join(parameters.join())))

    if spectrum is not None:
        ObservationModel = ObservationModel | Interpolate(spectrum)
        if normalize_npol is not None:
            ObservationModel = ObservationModel | Normalize(spectrum,
                                                            normalize_npol)


    ObservationModel.__class__.fit_parameters = property(fit_parameters_property)
    ObservationModel.rename('ObservationModel')

    return ObservationModel





assemble_model.__doc__ = assemble_model.__doc__.format(
    stellar_operations=','.join(stellar_operations),
    instrument_operations=','.join(instrument_operations))