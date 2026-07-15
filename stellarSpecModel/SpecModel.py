import os
import json
import hashlib
import numpy as np
from pathlib import Path
from scipy.interpolate import RegularGridInterpolator
from tqdm.auto import tqdm
from spectool import pyrebin
from SpecGrid import SpecGrid
# from .excepts import AliasAlreadyExistsError
# from . import config
import config


class SpecModel:
    def __init__(self, grid):
        """
        :param grid: SpecGrid 实例，包含底层的数据结构和元信息
        """
        self.grid = grid

    @classmethod
    def load(cls, filepath):
        """
        类方法：直接从 HDF5 文件实例化模型，方便研究者之间的数据流转分享。
        """
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Derived grid file not found at {filepath}")
        
        # 依赖 SpecGrid 自身的 HDF5 读取能力
        loaded_grid = SpecGrid.from_hdf5(filepath)
        return cls(loaded_grid)

    def _generate_cache_key(self, select: dict, wavelength: dict) -> str:
        """
        基于规范化后的 select 和 wavelength 生成唯一哈希指纹
        """
        # 1. 规范化 select
        norm_select = {}
        for k in sorted(select.keys()):
            val = select[k]
            # 统一转为 list 保证 JSON 序列化的确定性
            if isinstance(val, (tuple, list, np.ndarray)):
                norm_select[k] = np.asarray(val).tolist()
            else:
                norm_select[k] = val

        # 2. 规范化 wavelength
        norm_wave = {}
        if "range" in wavelength:
            norm_wave["range"] = list(wavelength["range"])
            
        if "resample" in wavelength:
            resamp = wavelength["resample"]
            if isinstance(resamp, dict):
                norm_wave["resample"] = {k: resamp[k] for k in sorted(resamp.keys())}
            elif isinstance(resamp, (list, np.ndarray)):
                # 如果是自定义波长列，提取字节序列进行哈希，避免 JSON 过长
                arr_bytes = np.asarray(resamp).tobytes()
                norm_wave["resample_hash"] = hashlib.md5(arr_bytes).hexdigest()

        # 3. 生成 JSON 字符串并 Hash
        fingerprint = {
            "selection": norm_select,
            "wavelength": norm_wave
        }
        fingerprint_str = json.dumps(fingerprint, sort_keys=True)
        return hashlib.md5(fingerprint_str.encode('utf-8')).hexdigest()

    def _create_symlink(self, target_path, symlink_path, overwrite):
        """创建操作系统软链接，作为预设网格的别名"""
        try:
            if os.path.exists(symlink_path) or os.path.islink(symlink_path):
                if overwrite:
                    os.remove(symlink_path)
                else:
                    return
            os.symlink(os.path.abspath(target_path), symlink_path)
            print(f"Symlink created: {os.path.basename(symlink_path)} -> {os.path.basename(target_path)}")
        except OSError as e:
            print(f"Warning: Failed to create symlink ({e}).")

    def _valid_wavelength(self, wavelength: dict | None = None):
        if wavelength is None:
            wave_range = (np.min(self.grid.wave), np.max(self.grid.wave))
            flag_resample = False
            method = 'custom'
            new_wave = None
            return wave_range, flag_resample, method, new_wave

        wave_range = wavelength.get("range")
        method = wavelength.get("method")
        step = wavelength.get("step")
        grid = wavelength.get("grid")

        # grid cannot coexist with other wavelength-generation options
        if grid is not None:
            if wave_range is not None:
                raise ValueError("'grid' cannot be used together with 'range'.")
            if method is not None or step is not None:
                raise ValueError("'grid' cannot be used together with 'method' or 'step'.")

            wave_range = (np.min(grid), np.max(grid))
            flag_resample = True
            return wave_range, flag_resample, method, grid

        # method and step must appear together
        if (method is None) != (step is None):
            raise ValueError("'method' and 'step' must be specified together.")

        if wave_range is None:
            wave_range = (np.min(self.grid.wave), np.max(self.grid.wave))
            
        # validate method
        if method is None:
            flag_resample = False
            method = 'custom'
            new_wave = None
            return wave_range, flag_resample, method, new_wave
        elif method == 'linear':
            min_w, max_w = wave_range
            new_wave = np.linspace(min_w, max_w, num=int((max_w - min_w) / step) + 1, endpoint=True)
            flag_resample = True
            return wave_range, flag_resample, method, new_wave
        elif method == 'log':
            min_w, max_w = wave_range
            log_min_w, log_max_w = np.log(min_w), np.log(max_w)
            log_new_wave = np.linspace(log_min_w, log_max_w, num=int((log_max_w - log_min_w) / step) + 1, endpoint=True)
            new_wave = np.exp(log_new_wave)
            flag_resample = True
            return wave_range, flag_resample, method, new_wave
        else:
            raise ValueError(f"Unsupported resampling method '{method}'.")

    def derive(
        self,
        select: dict | None = None,
        wavelength: dict | None = None,
        *,
        alias: str | None = None,
        cache_dir: str | Path | None = None,
        overwrite: bool = False,
        progress: bool = False,
    ) -> "SpecModel":
        """
        Create a derived spectral model from the current model.

        The ``derive`` method creates a new spectral grid by selecting a
        subset of the original parameter grid and/or modifying the wavelength
        sampling. It does not generate new physical models or interpolate
        along model parameter axes.

        Model parameters (e.g. ``teff``, ``logg``, ``feh``) can only be
        selected from existing grid points. The wavelength axis is treated
        differently and supports both wavelength selection and resampling.

        Parameters
        ----------
        select : dict, optional
            Selection rules for model parameters.

            The dictionary keys must correspond to existing parameter axes
            in the current model grid. Each parameter selection can be
            specified either as a range or as an explicit list of grid values.

            Examples
            --------
            Select a parameter range:

            >>> model.derive(
            ...     select={
            ...         "teff": (4500, 7000),
            ...         "logg": (2.0, 5.0)
            ...     }
            ... )

            Select specific existing grid points:

            >>> model.derive(
            ...     select={
            ...         "teff": numpy.ndarray([4500, 5000, 5500]),
            ...         "logg": (3.0, 4.0)
            ...     }
            ... )

            Notes
            -----
            - Parameters not included in ``select`` are preserved.
            - No interpolation is performed on model parameter axes.
            - Different spectral model families may contain different
            parameters. For example, a white dwarf model may only contain
            ``teff`` and ``logg``, while a stellar atmosphere model may
            additionally contain ``feh`` or ``alphaFe``.

        wavelength : dict, optional
            Operations applied to the wavelength axis.

            Supported operations include:

            ``range``
                Select a wavelength interval while preserving the original
                wavelength sampling.

                Example:

                >>> wavelength={
                ...     "range": (4000, 8000)
                ... }

            ``resample``
                Generate a regular wavelength grid using a predefined sampling
                scheme.

                Example (linear sampling):

                >>> wavelength={
                ...     "resample": {
                ...         "method": "linear",
                ...         "step": 1.0
                ...     }
                ... }

                Example (logarithmic sampling):

                >>> wavelength={
                ...     "resample": {
                ...         "method": "log",
                ...         "step": 1e-4
                ...     }
                ... }

                Example (custom sampling):

                >>> wavelength={
                        "resample": numpy.linspace(4000, 8000, 1000)
                    }

            Notes
            -----
            The wavelength operation order is:

            1. Crop the wavelength range (if specified).
            2. Generate the target wavelength grid.
            3. Interpolate spectra onto the new wavelength grid.

            Unlike model parameters, wavelength sampling can be changed.

        alias : str, optional
            User-defined name for the derived model.

            The alias is stored in the model metadata and can be used for
            identification in model registries or cache systems.

        cache_dir : str or pathlib.Path, optional
            Directory used to store the derived model.

            If provided, the derived model will be saved in this directory.
            Existing derived models may be reused if the same derivation
            configuration is found in the cache.

            If ``None``, the default package cache directory is used.

        overwrite : bool, default=False
            Whether to overwrite an existing cached derived model.

            If ``False``, an existing cached model with the same derivation
            hash may be reused.

            If ``True``, the derived model will be regenerated and the
            existing cached file will be replaced.

        Returns
        -------
        SpecModel
            A new derived spectral model.

            The returned model contains provenance information describing
            the parent model and all derivation operations.

        Raises
        ------
        ValueError
            If an invalid parameter name, parameter value, or wavelength
            operation is provided.

            Examples
            --------
            >>> model.derive(
            ...     select={"feh": (-1, 0)}
            ... )
            ValueError:
                Parameter 'feh' is not available.
                Available parameters: ['teff', 'logg']

        Notes
        -----
        The derived model keeps a complete record of its origin:

        Examples
        --------
        >>> derived.metadata["derived_from"]
        {
            "parent": "MARCS",
            "selection": {
                "teff": [4500, 7000],
                "logg": [2, 5]
            },
            "wavelength": {
                "range": [4000, 9000],
                "resample": {
                    "method": "log",
                    "step": 1e-4
                }
            }
        }

        This provenance information is used for:

        - reproducibility,
        - cache hashing,
        - model registry management,
        - future plugin extensions.
        """

        select = select or {}
        wavelength = wavelength or {}

        if cache_dir is None:
            active_cache_dir = Path(config.cache_PATH).expanduser()
        else:
            active_cache_dir = Path(cache_dir).expanduser()
        
        active_cache_dir.mkdir(parents=True, exist_ok=True)
        cache_hash = self._generate_cache_key(select, wavelength)
        model_name = self.grid.metadata['model_name']
        cache_filename = f"{model_name}_derived_{cache_hash}.h5"
        cache_filepath = active_cache_dir / cache_filename

        if alias:
            alias_PATH = Path(config.alias_PATH).expanduser()
            alias_PATH.mkdir(parents=True, exist_ok=True)
            alias_filepath = alias_PATH / f'{alias}.h5'
            if alias_filepath.exists() and not overwrite:
                raise AliasAlreadyExistsError(
                    f"Alias '{alias}' already exists. Pass overwrite=True to overwrite."
                )

        if cache_filepath.exists() and not overwrite:
            print(f"Cache hit! Loading derived grid from {cache_filepath}...")
            derived_model = self.__class__.load(cache_filepath)
            if alias:
                self._create_symlink(cache_filepath, alias_filepath, overwrite=overwrite)
            return derived_model

        # check the select par
        select_values = {}
        for param, val in select.items():
            if param not in self.grid.axis_names:
                raise ValueError(
                    f"Parameter '{param}' is not available. "
                    f"Available parameters: {list(self.grid.axis_names)}"
                )
            if isinstance(val, (list, np.ndarray)):
                valid_pts = set(self.grid.axes[param])
                invalid_pts = set(np.asarray(val)) - valid_pts
                if invalid_pts:
                    raise ValueError(f"Grid points {invalid_pts} for '{param}' do not exist.")
                select_values[param] = np.asarray(val)
            elif isinstance(val, tuple):
                if len(val) != 2:
                    raise ValueError(f"Grid point range for '{param}' must be a tuple of (min, max).")
                min_pt, max_pt = val
                valid_pts = set(self.grid.axes[param])
                arg = (self.grid.axes[param] >= min_pt) & (self.grid.axes[param] <= max_pt)
                select_value = self.grid.axes[param][arg]
                if len(select_value) == 0:
                    raise ValueError(f"Grid point range for {param} ({min_pt}, {max_pt}) does not match any existing grid points.")
                select_values[param] = select_value
            else:
                raise ValueError(f"Invalid selection value for '{param}': {val}. Must be a tuple, list, or ndarray.")

        wave_range, flag_resample, method, new_wave = self._valid_wavelength(wavelength)

        slice_indices = []
        new_axes = {}

        # 2.1 提取各个恒星参数轴的索引
        for param in self.grid.axis_names:
            axis_vals = self.grid.axes[param]
            if param in select_values:
                val = select_values[param]
                mask = np.any([np.isclose(axis_vals, v) for v in val], axis=0)
            else:
                mask = np.ones_like(axis_vals, dtype=bool)
            keep_idx = np.where(mask)[0]
            slice_indices.append(keep_idx)
            new_axes[param] = axis_vals[keep_idx]

        # 2.2 提取波长轴的索引
        wave_vals = self.grid.wave
        wave_mask = (wave_vals >= wave_range[0]) & (wave_vals <= wave_range[1])
        # if wave_range:
        #     wave_mask = (wave_vals >= wave_range[0]) & (wave_vals <= wave_range[1])
        # else:
        #     wave_mask = np.ones_like(wave_vals, dtype=bool)
        wave_keep_idx = np.where(wave_mask)[0]

        flux_slice_tuple = tuple(slice_indices + [wave_keep_idx])
        mask_slice_tuple = tuple(slice_indices)
        
        cropped_flux = self.grid.flux_tensor[np.ix_(*flux_slice_tuple)]
        cropped_mask = self.grid.valid_mask[np.ix_(*mask_slice_tuple)]
        cropped_wave = wave_vals[wave_keep_idx]
        grid_pars = self.grid.grid_parameters
        cropped_grid_pars = {param: grid_pars[param][np.ix_(*mask_slice_tuple)] for param in grid_pars}

        # ==========================================
        # 3. 波长重采样阶段 (只对切片后的小流量矩阵计算)
        # ==========================================
        # if target_wave is not None:
        if new_wave is not None:
            nflux_tensor = np.full(cropped_mask.shape + (len(new_wave),), np.nan, dtype=cropped_flux.dtype)
            iterator = np.ndindex(cropped_mask.shape)
            if progress:
                iterator = tqdm(
                    iterator,
                    total=cropped_mask.size,
                    desc="Resampling spectra"
                )
            for idx in iterator:
                if cropped_mask[idx]:
                    flux64 = 10**(cropped_flux[idx].astype(np.float64))
                    nflux64 = pyrebin.rebin_padvalue(cropped_wave, flux64, new_wave)
                    lognflux = np.log10(nflux64).astype(cropped_flux.dtype)
                    nflux_tensor[idx] = lognflux
        else:
            nflux_tensor = cropped_flux
            # final_wave = cropped_wave
            # final_flux = cropped_flux
            # space_type = self.grid.metadata.get('wave_sampling', 'unknown')

        # ==========================================
        # 4. 元数据血缘追踪与对象重装配
        # ==========================================
        new_metadata = self.grid.metadata.copy()
        new_metadata['is_derived'] = True
        # new_metadata['parent_hash'] = cache_key
        new_metadata['wave_sampling'] = method
        
        # 实例化新的底层网格
        new_grid = SpecGrid(
            wave=new_wave,
            axes=new_axes,
            axis_names=self.grid.axis_names,
            flux_tensor=nflux_tensor,
            valid_mask=cropped_mask,
            grid_parameters=cropped_grid_pars,
            metadata=new_metadata
        )
        
        new_model = self.__class__(new_grid)

        # ==========================================
        # 5. 数据落盘与别名绑定阶段
        # ==========================================
        # if cache_dir is not None:
        print(f"Caching derived grid to {cache_filepath}...")
        new_model.grid.to_hdf5(cache_filepath) # 直接利用底层 Grid 的 I/O 接口
        
        if alias:
            self._create_symlink(cache_filepath, alias_filepath, overwrite)

        return new_model

    def get_flux(self, **kwargs):
        missing_params = set(self.grid.axis_names) - set(kwargs.keys())
        if missing_params:
            raise ValueError(
                f"Missing required grid parameters: {list(missing_params)}. "
                f"Required parameters for this model are: {list(self.grid.axis_names)}"
            )
        query_point = []
        for param in self.grid.axis_names:
            val = kwargs[param]
            axis_arr = self.grid.axes[param]
            min_val, max_val = np.min(axis_arr), np.max(axis_arr)
            
            if val < min_val or val > max_val:
                raise ValueError(
                    f"Parameter '{param}'={val} is outside of the grid range [{min_val}, {max_val}]."
                )
            query_point.append(val)
        
        if not hasattr(self, '_interpolator'):
            points = tuple(self.grid.axes[param] for param in self.grid.axis_names)
            self._interpolator = RegularGridInterpolator(
                points, 
                self.grid.flux_tensor, 
                method='linear', 
                bounds_error=True, 
                fill_value=np.nan
            )
        
        interpolated_flux = self._interpolator([query_point])[0]
        if np.any(np.isnan(interpolated_flux)):
            raise ValueError(
                f"The requested parameters {kwargs} fall into a physical hole (invalid model region) in the grid."
            )

        return 10**interpolated_flux
