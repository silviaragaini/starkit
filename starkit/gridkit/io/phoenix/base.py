import os

import fnmatch

import h5py
from astropy.io import fits
from astropy import units as u

from starkit.gridkit.io.base import BaseSpectralGridIO
from starkit.gridkit.io.phoenix.alchemy import (Spectrum, ParameterSet,
                                                PhoenixBase)

from starkit.gridkit.io.phoenix.plugin import PhoenixProcess

class PhoenixSpectralGridIO(BaseSpectralGridIO):

    GridBase = PhoenixBase

    spectrum_table = Spectrum
    parameter_set_table = ParameterSet

    wavelength_unit = u.angstrom



    def _load_wavelength_solution(self, fname):
        return fits.getdata(fname, ext=0) * self.wavelength_unit

    def ingest(self):
        self.read_spectra()


    def read_spectra(self):
        spectral_files = self.get_spectral_files()
        no_of_spectra = len(spectral_files)
        for i, relative_path in enumerate(spectral_files):
            fname = os.path.basename(relative_path)
            fpath = os.path.relpath(os.path.dirname(relative_path),
                                    self.base_dir)
            print "{0} / {1} importing {2}".format(i, no_of_spectra, fname)
            spectrum = Spectrum(fname=fname, fpath=fpath)
            parameter_set = ParameterSet.from_file(spectrum.full_path)
            parameter_set.spectrum = spectrum
            self.session.add_all([spectrum, parameter_set])
        self.session.commit()




    def get_spectral_files(self):
        spectral_files = []
        for root, dirs, files in os.walk(self.base_dir):
            if 'grid/Z' not in root:
                continue
            for filename in fnmatch.filter(files, '*.fits'):
                spectral_files.append(os.path.join(root, filename))

        return spectral_files

    def to_hdf(self, fname, filter_tuple, R, wavelength_range, clobber=False):
        plugin = PhoenixProcess(self.wavelength.value, R, wavelength_range)
        super(PhoenixSpectralGridIO, self).to_hdf(fname, filter_tuple, plugin,
                                                  clobber)

        with h5py.File(fname, mode='a') as fh:
            fh['wavelength'].attrs['grid'] = 'log'
            fh['wavelength'].attrs['R'] =  R





