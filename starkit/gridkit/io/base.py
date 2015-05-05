import os
import logging
from abc import ABCMeta, abstractmethod

import numpy as np

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

import pandas as pd
import h5py
from specutils import Spectrum1D

from astropy import units as u

logger = logging.getLogger(__name__)

class BaseSpectralGridIO(object):
    """
    Base class for reading in the spectral grid information
    """

    __metaclass__ = ABCMeta

    GridBase = None

    def __init__(self, db_url, base_dir,
                 wavelength_fname='WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'):
        self.engine = create_engine(db_url)
        self.GridBase.metadata.create_all(self.engine)
        self.GridBase.metadata.bind = self.engine

        self.session = sessionmaker(bind=self.engine)()

        self._set_grid_base_dir(base_dir)

        self.wavelength = self._load_wavelength_solution(wavelength_fname)

        self._set_spectrum_wavelength(self.wavelength)



    def _set_grid_base_dir(self, base_dir):
        setattr(self.spectrum_table, 'base_dir', base_dir)
        self.base_dir = base_dir


    def _set_spectrum_wavelength(self, wavelength):
        setattr(self.spectrum_table, 'wavelength', wavelength)


    @staticmethod
    def _load_wavelength_solution():
        raise NotImplementedError('This needs to be implemented in subclasses')


    @staticmethod
    def _create_compound_model(models):
        if len(models) == 0:
            return lambda spectrum: spectrum
        else:
            compound_model = models[0]
            for model in models[1:]:
                compound_model = compound_model | model
            return compound_model

    def get_spectrum_query(self, filter_tuple):
            return self.session.query(self.spectrum_table).join(
            self.parameter_set_table).filter(*filter_tuple)


    def get_query_data(self, filter_tuple, plugin,
                       warning_threshold=1 * u.gigabyte):
        """
        Write spectra to disk
        :param filter_tuple:
        :param models:
        :return:
        """

        query = self.get_spectrum_query(filter_tuple)

        sample_spectrum_row = query.first()
        sample_spectrum_flux = plugin(sample_spectrum_row.get_spectrum1d().flux)

        no_spectra = query.count()

        size_of_spectra = (query.count() *
                           len(sample_spectrum_flux)) * 8 * u.byte

        if size_of_spectra > warning_threshold:
            continue_query = raw_input('The size of the spectra are {0:.2f}. '
                                       'Continue [y/N]'.format(
                size_of_spectra.to(warning_threshold.unit)))
            if continue_query.strip().lower() != 'y':
                raise ValueError('Size of requested grid ({:.2f}) to '
                                 'large for user ... aborting'.format(
                    size_of_spectra.to(warning_threshold.unit)))

        fluxes = np.empty((query.count(),
                          len(sample_spectrum_flux)))
        parameters = []
        param_names = [item.name
                       for item in sample_spectrum_row.parameter_set.parameters]

        for i, spectrum_row in enumerate(query):
            print "{0}/{1}".format(i, no_spectra)
            spectrum = spectrum_row.get_spectrum1d()
            fluxes[i] = plugin(spectrum.flux)
            parameters.append([getattr(spectrum_row.parameter_set, key)
                               for key in param_names])

        parameters = pd.DataFrame(parameters, columns= param_names)
        output_sample_spectrum = Spectrum1D.from_array(
            plugin.output_wavelength * u.angstrom, sample_spectrum_flux)

        return output_sample_spectrum, parameters, fluxes

    def to_hdf(self, fname, filter_tuple, plugin, clobber=False):
        sample_spectrum, parameters, fluxes = self.get_query_data(filter_tuple,
                                                                  plugin)

        if os.path.exists(fname):
            if clobber:
                os.remove(fname)
            else:
                raise IOError('File {0} exists - '
                              'if you want overwrite set clobber=True'.format(fname))

        parameters.to_hdf(fname, 'index')



        with h5py.File(fname, 'a') as fh:
            fh['fluxes'] = fluxes
            fh['fluxes'].attrs['unit'] = str(self.spectrum_table.flux_unit)
            fh['wavelength'] = sample_spectrum.wavelength.value
            fh['wavelength'].attrs['unit'] = str(
                sample_spectrum.wavelength.unit)
