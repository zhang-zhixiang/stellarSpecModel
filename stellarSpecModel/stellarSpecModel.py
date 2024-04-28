from . import config
import os
import scipy.interpolate as spinterp
from astropy import units as u
import h5py


class StellarSpecModel:
    """
    Represents a general stellar spectral model.

    This class provides functionality to handle and interpolate stellar spectral data.

    Attributes:
        None
    """

    def __init__(self, grid_name):
        """
        Initialize the StellarSpecModel.

        Args:
            grid_name (str): Name of the spectral grid.

        Returns:
            StellarSpecModel: An instance of the StellarSpecModel class.
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
        """Get the wavelength array."""
        return self._wavelength

    @property
    def teff_grid(self):
        """Get the effective temperature (Teff) grid."""
        return self._teff_grid

    @property
    def min_teff(self):
        """Get the minimum Teff in the grid."""
        return self._teff_grid.min()

    @property
    def max_teff(self):
        """Get the maximum Teff in the grid."""
        return self._teff_grid.max()

    @property
    def feh_grid(self):
        """Get the metallicity (FeH) grid."""
        return self._feh_grid

    @property
    def min_feh(self):
        """Get the minimum FeH in the grid."""
        return self._feh_grid.min()
    
    @property
    def max_feh(self):
        """Get the maximum FeH in the grid."""
        return self._feh_grid.max()

    @property
    def logg_grid(self):
        """Get the surface gravity (logg) grid."""
        return self._logg_grid

    @property
    def min_logg(self):
        """Get the minimum logg in the grid."""
        return self._logg_grid.min()

    @property
    def max_logg(self):
        """Get the maximum logg in the grid."""
        return self._logg_grid.max()

    def get_flux(self, teff, feh, logg):
        """
        Get the flux for a given set of Teff, FeH, and logg values.

        Args:
            teff (float): Effective temperature (Teff).
            feh (float): Metallicity (FeH).
            logg (float): Surface gravity (logg).

        Returns:
            numpy.ndarray: Flux array.
        """
        if teff < self.min_teff or teff > self.max_teff:
            raise ValueError('Teff outside of grid range')
        if feh < self.min_feh or feh > self.max_feh:
            raise ValueError('FeH outside of grid range')
        if logg < self.min_logg or logg > self.max_logg:
            raise ValueError('logg outside of grid range')
        log_flux = self._model((teff, feh, logg))
        return 10.0 ** log_flux

    @property
    def flux_units(self):
        """Get the flux units."""
        return self._flux_units

    @property
    def wavelength_units(self):
        """Get the wavelength units."""
        return self._wavelength_units


class MARCS_Model(StellarSpecModel):
    """
    Represents the MARCS stellar spectral model.

    This class specializes the StellarSpecModel for MARCS model.

    Attributes:
        None
    """

    def __init__(self):
        """
        Initialize the MARCS_Model.

        Args:
            None

        Returns:
            MARCS_Model: An instance of the MARCS_Model class.
        """
        grid_name = 'MARCS'
        file_name, url, md5_value = config.grid_names[grid_name]
        abs_filename = os.path.join(config.grid_data_dir, file_name)
        if not os.path.exists(abs_filename):
            config.fetch_grid(grid_name)
        super().__init__(abs_filename)


class MARCS_Model_hiRes(StellarSpecModel):
    """
    Represents the MARCS stellar spectral model.

    This class specializes the StellarSpecModel for MARCS model with higher resolution.
    The resolution of this grid is corresponding to a velocity dispersion of 4 km/s.

    Attributes:
        None
    """
    def __init__(self):
        """
        Initialize the MARCS_Model_hiRes.

        Args:
            None

        Returns:
            MARCS_Model_hiRes: An instance of the MARCS_Model_hiRes class.
        """
        grid_name = 'MARCS_hiRes'
        file_name, url, md5_value = config.grid_names[grid_name]
        abs_filename = os.path.join(config.grid_data_dir, file_name)
        if not os.path.exists(abs_filename):
            config.fetch_grid(grid_name)
        super().__init__(abs_filename)


class BTCond_Model(StellarSpecModel):
    """
    Represents the BTCond stellar spectral model.

    This class specializes the StellarSpecModel for BTCond model.

    Attributes:
        None
    """

    def __init__(self):
        """
        Initialize the BTCond_Model.

        Args:
            None

        Returns:
            BTCond_Model: An instance of the BTCond_Model class.
        """
        grid_name = 'BTCond'
        file_name, url, md5_value = config.grid_names[grid_name]
        abs_filename = os.path.join(config.grid_data_dir, file_name)
        if not os.path.exists(abs_filename):
            config.fetch_grid(grid_name)
        super().__init__(abs_filename)


class BTCond_Model_hiRes(StellarSpecModel):
    """
    Represents the BTCond stellar spectral model.

    This class specializes the StellarSpecModel for BTCond model.

    Attributes:
        None
    """

    def __init__(self):
        """
        Initialize the BTCond_Model.

        Args:
            None

        Returns:
            BTCond_Model: An instance of the BTCond_Model class.
        """
        grid_name = 'BTCond_hiRes'
        file_name, url, md5_value = config.grid_names[grid_name]
        abs_filename = os.path.join(config.grid_data_dir, file_name)
        if not os.path.exists(abs_filename):
            config.fetch_grid(grid_name)
        super().__init__(abs_filename)


class BTCond_Model_R7500(StellarSpecModel):
    """
    Represents the BTCond stellar spectral model.

    This class specializes the StellarSpecModel for BTCond model.

    Attributes:
        None
    """

    def __init__(self):
        """
        Initialize the BTCond_Model.

        Args:
            None

        Returns:
            BTCond_Model: An instance of the BTCond_Model class.
        """
        grid_name = 'BTCond_R7500'
        file_name, url, md5_value = config.grid_names[grid_name]
        abs_filename = os.path.join(config.grid_data_dir, file_name)
        if not os.path.exists(abs_filename):
            config.fetch_grid(grid_name)
        super().__init__(abs_filename)


class BTCond_Model_R1800(StellarSpecModel):
    """
    Represents the BTCond stellar spectral model.

    This class specializes the StellarSpecModel for BTCond model.

    Attributes:
        None
    """

    def __init__(self):
        """
        Initialize the BTCond_Model.

        Args:
            None

        Returns:
            BTCond_Model: An instance of the BTCond_Model class.
        """
        grid_name = 'BTCond_R1800'
        file_name, url, md5_value = config.grid_names[grid_name]
        abs_filename = os.path.join(config.grid_data_dir, file_name)
        if not os.path.exists(abs_filename):
            config.fetch_grid(grid_name)
        super().__init__(abs_filename)


class BTCond_Model_R500(StellarSpecModel):
    """
    Represents the BTCond stellar spectral model.

    This class specializes the StellarSpecModel for BTCond model.

    Attributes:
        None
    """

    def __init__(self):
        """
        Initialize the BTCond_Model.

        Args:
            None

        Returns:
            BTCond_Model: An instance of the BTCond_Model class.
        """
        grid_name = 'BTCond_R500'
        file_name, url, md5_value = config.grid_names[grid_name]
        abs_filename = os.path.join(config.grid_data_dir, file_name)
        if not os.path.exists(abs_filename):
            config.fetch_grid(grid_name)
        super().__init__(abs_filename)


class BTCond_Model_R100(StellarSpecModel):
    """
    Represents the BTCond stellar spectral model.

    This class specializes the StellarSpecModel for BTCond model.

    Attributes:
        None
    """

    def __init__(self):
        """
        Initialize the BTCond_Model.

        Args:
            None

        Returns:
            BTCond_Model: An instance of the BTCond_Model class.
        """
        grid_name = 'BTCond_R100'
        file_name, url, md5_value = config.grid_names[grid_name]
        abs_filename = os.path.join(config.grid_data_dir, file_name)
        if not os.path.exists(abs_filename):
            config.fetch_grid(grid_name)
        super().__init__(abs_filename)