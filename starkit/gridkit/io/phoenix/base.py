import os
import re
import fnmatch

from astropy.io import fits
from astropy import units as u

from starkit.gridkit.io.base import BaseSpectralGridIO
from starkit.gridkit.io.phoenix.alchemy import (Spectrum, ParameterSet,
                                                PhoenixBase)

class PhoenixSpectralGridIO(BaseSpectralGridIO):

    spectrum_table = Spectrum
    parameter_set_table = ParameterSet

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

        current_param = ParameterSet(teff=teff, logg=logg, feh=feh)
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
