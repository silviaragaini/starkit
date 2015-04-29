import os
import fnmatch
import re

from astropy import units as u
from astropy.io import fits
from specutils import Spectrum1D
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import backref, relationship
from sqlalchemy import Column, ForeignKey
from sqlalchemy import Integer, Float, String

from starkit.gridkit.io.spectral_grid import BaseSpectralGridIO

PhoenixBase = declarative_base()



class Spectrum(PhoenixBase):
    __tablename__ = 'spectra'

    id = Column(Integer, primary_key=True)
    fpath = Column(String)
    fname = Column(String)

    base_dir = None
    wavelength = None
    flux_unit = 'erg/s/cm2/angstrom'

    @property
    def full_path(self):
        if self.base_dir is None:
            raise AttributeError('base_dir attribute needs to be set')
        else:
            return os.path.join(self.base_dir, self.fpath, self.fname)

    def __repr__(self):
        return '<Spectrum {0} with {1}>'.format(
            self.fname, ' '.join(self.parameters.to_string()))

    def get_spectrum1d(self):
        flux = self._read_flux()
        return Spectrum1D.from_array(self.wavelength, flux,
                                     unit=u.Unit(self.flux_unit))

    def _read_flux(self):
        return fits.getdata(self.full_path)


class Parameter(PhoenixBase):
    __tablename__ = 'parameters'

    id = Column(Integer, ForeignKey('spectra.id'), primary_key=True)
    meta_id = Column(Integer, ForeignKey('meta_parameters.id'))

    teff = Column(Float)
    logg = Column(Float)
    feh = Column(Float)

    spectrum = relationship('Spectrum', uselist=False,
                            backref=backref('parameters', uselist=False))

    param_names = ['teff', 'logg', 'feh']

    def to_string(self):
        return ['{0}={1}'.format(param_name, getattr(self, param_name))
         for param_name in ['teff', 'logg', 'feh']]


class MetaParameter(PhoenixBase):
    __tablename__ = 'meta_parameters'

    id = Column(Integer, primary_key=True)


class PhoenixSpectralGridIO(BaseSpectralGridIO):

    spectrum_table = Spectrum
    parameter_table = Parameter

    @staticmethod
    def _set_grid_base_dir(base_dir):
        setattr(Spectrum, 'base_dir', base_dir)

    @staticmethod
    def _set_spectrum_wavelength(wavelength):
        setattr(Spectrum, 'wavelength', wavelength)

    @staticmethod
    def _load_wavelength_solution(fname, unit, type):
        return fits.getdata(fname, ext=0) * u.Unit(unit)



    def initialize_database(self, config_dict):
        self.engine = self.load_engine(init_db=True, **config_dict['database'])
        self.initialize_session()
        PhoenixBase.metadata.create_all(self.engine)
        spectral_files = self.get_spectral_files(config_dict['base_dir'])

        self.read_spectra(spectral_files, config_dict['base_dir'])

    def read_spectra(self, spectral_files, base_dir):
        no_of_spectra = len(spectral_files)
        for i, relative_path in enumerate(spectral_files):
            fname = os.path.basename(relative_path)
            fpath = os.path.relpath(os.path.dirname(relative_path), base_dir)
            print "{0} / {1} importing {2}".format(i, no_of_spectra, fname)
            spectrum = Spectrum(fname=fname, fpath=fpath)
            parameter = self.get_parameter_object(fpath, fname)
            parameter.spectrum = spectrum
            self.session.add_all([spectrum, parameter])
        self.session.commit()


    @staticmethod
    def get_parameter_object(fpath, fname):
        parameter_pattern = re.compile('lte(\d+)-(\d+\.\d+)([+-]\d+\.\d+)'
                                       '\.PHOENIX.')
        teff, logg, feh = map(float, parameter_pattern.match(fname).groups())

        current_param = Parameter(teff=teff, logg=logg, feh=feh)
        return current_param


    @staticmethod
    def get_spectral_files(base_dir):
        spectral_files = []
        for root, dirs, files in os.walk(base_dir):
            if 'grid/Z' not in root:
                continue
            for filename in fnmatch.filter(files, '*.fits'):
                spectral_files.append(os.path.join(root, filename))

        return spectral_files


