from . import config
import scipy.interpolate as spinterp
from astropy import units as u
import h5py


class StellarSpecModel:
    """
    Class to hold the stellar spectral model
    """

    def __init__(self, grid_name):
        """
        Constructor for StellarSpecModel

        Parameters
        ----------
        kwargs : dict
            Dictionary of keyword arguments

        Returns
        -------
        StellarSpecModel object
        """
        self._grid_name = grid_name
        grids = h5py.File(self._grid_name, 'r')
        grid = grids['default']
        wave = grid['wave'].astype(float)[:]
        teff_grid = grid['teff'].astype(float)[:]
        feh_grid = grid['feh'].astype(float)[:]
        logg_grid = grid['logg'].astype(float)[:]
        spec_grid = grid['spec_grid'].astype(float)[:]
        grids.close()
        self._wavelength = wave
        self._teff_grid = teff_grid
        self._feh_grid = feh_grid
        self._logg_grid = logg_grid
        self._spec_grid = spec_grid
        self._model = spinterp.RegularGridInterpolator((teff_grid, feh_grid, logg_grid), spec_grid)
        self._flux_units = u.erg / u.s / u.cm ** 2 / u.AA
        self._wavelength_units = u.AA

    @property
    def wavelength(self):
        """
        Wavelength array

        Returns
        -------
        numpy.ndarray
        """
        return self._wavelength

    @property
    def teff_grid(self):
        """
        Teff grid

        Returns
        -------
        numpy.ndarray
        """
        return self._teff_grid

    @property
    def teff_min(self):
        """
        Minimum Teff

        Returns
        -------
        float
        """
        return self._teff_grid.min()

    @property
    def teff_max(self):
        """
        Maximum Teff

        Returns
        -------
        float
        """
        return self._teff_grid.max()

    @property
    def feh_grid(self):
        """
        FeH grid

        Returns
        -------
        numpy.ndarray
        """
        return self._feh_grid

    @property
    def feh_min(self):
        """
        Minimum FeH

        Returns
        -------
        float
        """
        return self._feh_grid.min()
    
    @property
    def feh_max(self):
        """
        Maximum FeH

        Returns
        -------
        float
        """
        return self._feh_grid.max()

    @property
    def logg_grid(self):
        """
        logg grid

        Returns
        -------
        numpy.ndarray
        """
        return self._logg_grid

    @property
    def logg_min(self):
        """
        Minimum logg

        Returns
        -------
        float
        """
        return self._logg_grid.min()

    @property
    def logg_max(self):
        """
        Maximum logg

        Returns
        -------
        float
        """
        return self._logg_grid.max()

    def get_flux(self, teff, feh, logg):
        """
        Get flux for given Teff, FeH, logg

        Parameters
        ----------
        teff : float
            Teff
        feh : float
            FeH
        logg : float
            logg

        Returns
        -------
        numpy.ndarray
        """
        if teff < self.teff_min or teff > self.teff_max:
            raise ValueError('Teff outside of grid range')
        if feh < self.feh_min or feh > self.feh_max:
            raise ValueError('FeH outside of grid range')
        if logg < self.logg_min or logg > self.logg_max:
            raise ValueError('logg outside of grid range')
        log_flux = self._model((teff, feh, logg))
        return 10.0 ** log_flux

    @property
    def flux_units(self):
        """
        Flux units

        Returns
        -------
        astropy.units
        """
        return self._flux_units

    @property
    def wavelength_units(self):
        """
        Wavelength units

        Returns
        -------
        astropy.units
        """
        return self._wavelength_units


class MARCS_Model(StellarSpecModel):
    """
    Class to hold the MARCS stellar spectral model
    """

    def __init__(self):
        """
        Constructor for MARCS_Model

        Parameters
        ----------
        kwargs : dict
            Dictionary of keyword arguments

        Returns
        -------
        MARCS_Model object
        """
        grid_name = config.grid_data_dir + config.grid_names['MARCS']
        super().__init__(grid_name)


class BTCond_Model(StellarSpecModel):
    """
    Class to hold the BTCond stellar spectral model
    """

    def __init__(self):
        """
        Constructor for BTCond_Model

        Parameters
        ----------
        kwargs : dict
            Dictionary of keyword arguments

        Returns
        -------
        BTCond_Model object
        """
        grid_name = config.grid_data_dir + config.grid_names['BTCond']
        super().__init__(grid_name)