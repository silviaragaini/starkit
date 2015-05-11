from starkit.base.operations import stellar_operations, instrument_operations

def assemble_model(spectral_grid, spectral_operations, spectrum=None,
                   normalize_npol=None):
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
    for operation in spectral_operations:
        if operation in stellar_operations:
            ObservationModel = ObservationModel | stellar_operations[operation]

        elif operation in instrument_operations:
            ObservationModel = ObservationModel | instrument_operations[
                                                     operation]

    return ObservationModel




assemble_model.__doc__ = assemble_model.__doc__.format(
    stellar_operations=','.join(stellar_operations),
    instrument_operations=','.join(instrument_operations))