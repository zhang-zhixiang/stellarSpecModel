from . import config
from .stellarSpecModel import StellarSpecModel
import h5py
import scipy.interpolate as spinterp
from astropy import units as u
import numpy as np
import os


class TlustyModel(StellarSpecModel):

    def __init__(self):
        self._flux_units = u.erg / u.s / u.cm ** 2 / u.AA
        self._wavelength_units = u.AA
        grid_name = 'TLUSTY'
        fname, url, md5_value = config.grid_names[grid_name]
        abs_filename = os.path.join(config.grid_data_dir, fname)
        h5grids = h5py.File(abs_filename, 'r')
        grids = []
        loggs_left = []
        loggs_right = []
        models = []
        all_teffs = np.array([])
        all_loggs = np.array([])
        all_fehs = np.array([])
        for gname in h5grids:
            grid = h5grids[gname]
            waves = grid['wave'].astype(float)[:]
            teffs = grid['teff'].astype(float)[:]
            loggs = grid['logg'].astype(float)[:] / 100
            fehs = grid['z'].astype(float)[:]
            spec_grids = grid['spec_grid'].astype(float)[:]
            model = spinterp.RegularGridInterpolator((teffs, fehs, loggs), spec_grids)
            models.append(model)
            loggs_left.append(loggs.min())
            loggs_right.append(loggs.max())
            grids.append([teffs, fehs, loggs])
            all_teffs = np.concatenate((all_teffs, teffs))
            all_loggs = np.concatenate((all_loggs, loggs))
            all_fehs = np.concatenate((all_fehs, fehs))
        self._teff_grid = np.unique(all_teffs)
        self._logg_grid = np.unique(all_loggs)
        self._feh_grid = np.unique(all_fehs)
        np.sort(self._teff_grid)
        np.sort(self._logg_grid)
        np.sort(self._feh_grid)
        h5grids.close()
        self._wavelength = waves
        self._grids = grids
        self._models = models
        self._loggs_left = np.array(loggs_left)
        self._loggs_right = np.array(loggs_right)

    def get_flux(self, teff, feh, logg):
        arg = (logg >= self._loggs_left) & (logg < self._loggs_right)
        if not arg.any():
            raise ValueError(f'logg {logg} out of range')
        ind = np.where(arg)[0][0]
        model = self._models[ind]
        grid = self._grids[ind]
        teffs, fehs, loggs = grid
        if teff < teffs.min() or teff > teffs.max():
            raise ValueError(f'teff {teff} out of range')
        if feh < fehs.min() or feh > fehs.max():
            raise ValueError(f'feh {feh} out of range')
        log_flux = model((teff, feh, logg))
        flux = 10 ** log_flux
        return flux