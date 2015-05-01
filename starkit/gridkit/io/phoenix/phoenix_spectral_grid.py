import os
import fnmatch
import re



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

    id = Column(Integer, primary_key=True)
    parameter_set_id = Column(Integer, ForeignKey('parameter_sets.id'))
    parameter_type_id = Column(Integer, ForeignKey('parameter_type_id'))
    parameter_value = Column()




