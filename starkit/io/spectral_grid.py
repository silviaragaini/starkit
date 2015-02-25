import os
import logging


import yaml
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker


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
