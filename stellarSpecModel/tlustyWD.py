from . import config
from .stellarSpecModel import StellarSpecModel
import h5py
import scipy.interpolate as spinterp
from astropy import units as u
import numpy as np
import os
import spectool


class TlustyWDModel(StellarSpecModel):

    def __init__(self):
        self._flux_units = u.erg / u.s / u.cm ** 2 / u.AA
        self._wavelength_units = u.AA
        grid_name = 'TLUSTYWD'
        fname, url, md5_value = config.grid_names[grid_name]
        abs_filename = os.path.join(config.grid_data_dir, fname)
        h5grids = h5py.File(abs_filename, 'r')
        grid_name = 'default'
        grid = h5grids[grid_name]
        wave = grid['wave'].astype(float)[:]
        teff_grid = grid['tgrid'].astype(float)[:]
        logg_grid = grid['ggrid'].astype(float)[:]
        flux_grid = grid['flux'].astype(float)[:]
        model = spinterp.RegularGridInterpolator((teff_grid, logg_grid), np.log10(flux_grid.T))
        self._model = model
        self._wavelength = spectool.spec_func.air2vac(wave)
        self._teff_grid = teff_grid
        self._logg_grid = logg_grid

    def get_flux(self, teff, logg):
        logflux = self._model((teff, logg))
        flux = 10 ** logflux
        return flux