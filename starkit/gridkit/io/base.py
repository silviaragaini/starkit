import os
import logging


import yaml
import numpy as np

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker


import pandas as pd
import h5py

from astropy import units as u

logger = logging.getLogger(__name__)

class BaseSpectralGridIO(object):
    """
    Base class for reading in the spectral grid information
    """

    def __init__(self, yaml_fname):
        config_dict = yaml.load(open(yaml_fname))
        self.engine = self.load_engine(**config_dict['database'])
        if self.engine is None:
            self.initialize_database(config_dict)
        else:
            self.initialize_session()
        self._set_grid_base_dir(config_dict['base_dir'])

        self.wavelength = self._load_wavelength_solution(
            **config_dict['wavelength'])

        self._set_spectrum_wavelength(self.wavelength)

    @staticmethod
    def _load_wavelength_solution():
        raise NotImplementedError('This needs to be implemented in subclasses')

    @staticmethod
    def load_engine(fname, type='sqlite', init_db=False):
        if type!='sqlite3':
            raise ValueError('Currently only type "sqlite" is supported')

        if not os.path.exists(fname) and init_db is not True:
            logger.info('Database {0} does not exist'.format(fname))
            return None
        else:
            return create_engine('sqlite:///{0}'.format(fname))



    def initialize_session(self):
        self.session = sessionmaker(bind=self.engine)()

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
            self.parameter_table).filter(*filter_tuple)


    def get_query_data(self, filter_tuple, models=[],
                       warning_threshold=1 * u.gigabyte):
        """
        Write spectra to disk
        :param filter_tuple:
        :param models:
        :return:
        """
        compound_model = BaseSpectralGridIO._create_compound_model(models)

        query = self.get_spectrum_query(filter_tuple)

        sample_spectrum_row = query.first()
        sample_spectrum = compound_model(sample_spectrum_row.get_spectrum1d())

        no_spectra = query.count()

        size_of_spectra = (query.count() *
                           len(sample_spectrum.flux)) * 8 * u.byte

        if size_of_spectra > warning_threshold:
            continue_query = raw_input('The size of the spectra are {0:.2f}. '
                                       'Continue [y/N]'.format(
                size_of_spectra.to(warning_threshold.unit)))
            if continue_query.strip().lower() != 'y':
                raise ValueError('Size of requested grid ({:.2f}) to '
                                 'large for user ... aborting'.format(
                    size_of_spectra.to(warning_threshold.unit)))

        fluxes = np.empty((query.count(),
                          len(sample_spectrum.flux)))
        parameters = []
        param_names = sample_spectrum_row.parameters.param_names
        for i, spectrum_row in enumerate(query):
            print "{0}/{1}".format(i, no_spectra)
            spectrum = spectrum_row.get_spectrum1d()
            fluxes[i] = compound_model(spectrum).flux
            parameters.append([getattr(spectrum_row.parameters, key)
                               for key in param_names])

        parameters = pd.DataFrame(parameters, columns= param_names)

        return sample_spectrum, parameters, fluxes

    def to_hdf5(self, fname, filter_tuple, models=[], clobber=False):
        sample_spectrum, parameters, fluxes = self.get_query_data(filter_tuple, models)

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
