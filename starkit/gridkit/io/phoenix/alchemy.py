from astropy.io import fits

from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import backref, relationship
from sqlalchemy import Column, ForeignKey
from sqlalchemy import Integer, Float, String

from astropy import units as u
from astropy.io import fits
from specutils import Spectrum1D

from starkit.gridkit.io.phoenix import parameters

PhoenixBase = declarative_base()



import os


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


class ParameterSetMixin(object):
    ___tablename__ = 'parameter_sets'

def make_parameter_set():
    type_dict = {'float':Float}
    parameter_classes = parameters.BasePhoenixParameter.__subclasses__()
    class_dict = {'parameters': [item.name for item in parameter_classes]}

    for param in parameter_classes:
        class_dict[param.name] = Column(type_dict(param.type))
    class_dict = class_dict.update({item.name:Column()})

    return type('ParameterSet', (PhoenixBase, ParameterSetMixin), class_dict)

ParameterSet = make_parameter_set()