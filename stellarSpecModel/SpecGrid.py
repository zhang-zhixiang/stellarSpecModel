import numpy as np
import h5py
import datetime


class SpecGrid:
    def __init__(self, wave, axes, axis_names, flux_tensor, valid_mask=None, 
                 grid_parameters=None, metadata=None, h5_file=None):
        """_summary_
        Args:
            wave (np.ndarray): 1D array, wavelength
            axes (dict): info of each axis, example: {'teff': [3000, 3150...], 'logg': [2.0, ...]}
            axis_names (list): names of each axis, example: ['teff', 'feh', 'logg']
            flux_tensor (np.Ndarray): N+1 dimension flux tensor
            valid_mask (np.ndarray, optional): N dimension bool tensor, sign whether the flux of given subscript valid. Defaults to None.
            metadata (dict, optional): meta data of the grid, the grid info should like this:
            metadata = 
            {   'model_name': 'MARCS', 'PHOENIX', or 'Kurucz'
                'reference': DOI or URL
                'wave_sampling': 'log', 'linear', or 'custom'
                'parent_hash': if 'is_derived' is True, this key record the parent info of the grid
                'creation_date': 'YYYY-MM-DDTHH:MM:SS',
                'is_derived': False
            }. Defaults to None.
        """
        self.wave = np.asarray(wave)
        self.axes = axes
        self.axis_names = tuple(axis_names)
        self.flux_tensor = flux_tensor

        self._h5_file = h5_file
        
        if valid_mask is not None:
            self.valid_mask = np.asarray(valid_mask)
        else:
            self.valid_mask = np.ones(self.flux_tensor.shape[:-1], dtype=bool)

        self.grid_parameters = grid_parameters or {}

        self.metadata = metadata if metadata is not None else {}
        self._init_default_metadata()
        self._validate_dimensions()

    def _init_default_metadata(self):
        """inject creation time"""
        if 'creation_date' not in self.metadata:
            self.metadata['creation_date'] = datetime.datetime.now().isoformat()
        if 'is_derived' not in self.metadata:
            self.metadata['is_derived'] = False

    def _validate_dimensions(self):
        expected_grid_shape = tuple(len(self.axes[p]) for p in self.axis_names)
        expected_flux_shape = expected_grid_shape + (len(self.wave),)
        
        if self.valid_mask.shape != expected_grid_shape:
            raise ValueError(f"Mask shape {self.valid_mask.shape} mismatches axes shape {expected_grid_shape}")
            
        if self.flux_tensor.shape != expected_flux_shape:
            raise ValueError(f"Flux tensor shape {self.flux_tensor.shape} mismatches expected {expected_flux_shape}")

    @classmethod
    def from_hdf5(cls, filepath, lazy=True):
        """load grid and the meta data from a hdf5 file从 HDF5"""
        f = h5py.File(filepath, 'r')
        try:
            wave = f['wave'][:]
            
            axis_names = tuple(f.attrs['axis_names']) 
            
            axes = {}
            for name in axis_names:
                axes[name] = f[f'axes/{name}'][:]

            grid_parameters = {}
            if "grid_parameters" in f:
                group = f["grid_parameters"]
                for key in group.keys():
                    grid_parameters[key] = group[key][:]

            flux_tensor = f['flux_tensor'] if lazy else f['flux_tensor'][:]
            valid_mask = f['valid_mask'][:] if 'valid_mask' in f else None

            metadata = {}
            for key, value in f.attrs.items():
                if key not in ['axis_names']:
                    metadata[key] = value

            if not lazy:
                f.close()
                file_handle = None
            else:
                file_handle = f

            return cls(wave, axes, axis_names, flux_tensor, valid_mask, 
                       grid_parameters, metadata, h5_file=file_handle)
        
        except Exception as e:
            f.close()
            raise e

    def close(self):
        if getattr(self, '_h5_file', None) is not None:
            self._h5_file.close()
            self._h5_file = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def to_hdf5(self, filepath):
        """save grid and metadata to hdf5 file"""
        with h5py.File(filepath, 'w') as f:
            f.create_dataset('wave', data=self.wave)
            f.create_dataset('flux_tensor', data=self.flux_tensor)
            if self.valid_mask is not None:
                f.create_dataset('valid_mask', data=self.valid_mask)
            
            for name, array in self.axes.items():
                f.create_dataset(f'axes/{name}', data=array)
            
            f.attrs['axis_names'] = self.axis_names
            
            for key, value in self.metadata.items():
                f.attrs[key] = value

    @property
    def shape(self):
        return self.flux_tensor.shape[:-1]

    @property
    def ndim(self):
        return len(self.axis_names)

    @property
    def n_wave(self):
        return len(self.wave)