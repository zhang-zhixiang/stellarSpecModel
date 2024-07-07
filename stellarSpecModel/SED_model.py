import io
import contextlib
import pyphot
import spectool
import numpy as np
from extinction import apply
from extinction import fitzpatrick99
from astropy import constants as cs
from . import stellarSpecModel
from .phot_util import flux_to_mag as f2m
from .phot_util import mag_to_flux as m2f
from .phot_util import filtername2pyphotname
from .phot_util import fluxes_to_mags, mags_to_fluxes


class SEDModel:
    def __init__(self, bands=None, teff=5700, logg=4.5, feh=0.0, 
                 R=1.0, distance=10.0, Av=0.0,
                 specmodel=stellarSpecModel.BTCond_Model()):
        """a class to generate stellar model SED

        Args:
            bands (list, optional): the band names of the SED, a example ['SDSSu', 'JV', 'SDSSr', '2MASSJ', 'W2']. Defaults to None.
            teff (float, optional): temperature of the SED. Defaults to 5700.
            logg (float, optional): surface gravity. Defaults to 4.5.
            feh (float, optional): Metallicity. Defaults to 0.0.
            R (float, optional): radius of the stellar, unit is R_sun. Defaults to 1.0.
            distance (float, optional): distance of the stellar, unit is pc. Defaults to 10.0.
            Av (float, optional): Extinction Coefficient. Defaults to 0.0.
            specmodel (stellarSpecModel, optional): stellarSpecModel used to generate the stellar spectrum. Defaults to stellarSpecModel.BTCond_Model().

        Raises:
            ValueError: if specmodel is not an instance of StellarSpecModel, raise ValueError
        """
        if isinstance(specmodel, stellarSpecModel.StellarSpecModel):
            self.stellar_model = specmodel
        else:
            raise ValueError('specmodel must be an instance of StellarSpecModel')
        self.bands = []
        self.filters = {}
        self.pyphot_lib = pyphot.get_library()
        self._rat_rsun_pc = cs.R_sun.to('pc').value
        self.Rv = 3.1
        self.teff = teff
        self.logg = logg
        self.feh = feh
        self.distance = distance
        self.Av = Av
        self.rad = R
        if bands is not None:
            for band in bands:
                self.add_band(band)

    def set_SED_pars(self, teff, logg, feh, R, distance, Av=0.0):
        self.teff = teff
        self.logg = logg
        self.feh = feh
        self.distance = distance
        self.rad = R
        self.Av = Av

    def set_teff(self, teff):
        self.teff = teff

    def set_logg(self, logg):
        self.logg = logg

    def set_feh(self, feh):
        self.feh = feh

    def set_distance(self, distance):
        self.distance = distance

    def set_radius(self, R):
        self.rad = R

    def set_Av(self, Av):
        self.Av = Av
        
    def add_band(self, band):
        std_filtername = filtername2pyphotname(band)
        self.bands.append(std_filtername)
        tfilter = self.pyphot_lib[std_filtername]
        waves = tfilter.wavelength.to('AA').value
        trans = tfilter.transmit
        eff_wave = tfilter.leff.to('AA').value
        width = tfilter.width.to('AA').value
        self.filters[std_filtername] = [waves, trans, eff_wave, width]

    def add_bands(self, bands):
        for band in bands:
            self.add_band(band)

    def get_SED_spec(self):
        waves = self.stellar_model.wavelength
        fluxes = self.stellar_model.get_flux(self.teff, self.feh, self.logg)
        rat = (self.rad / self.distance * self._rat_rsun_pc) ** 2
        fluxes *= rat
        ext = fitzpatrick99(waves, self.Av, self.Rv)
        nfluxes = apply(ext, fluxes)
        return waves, nfluxes

    def get_SED(self):
        waves, fluxes = self.get_SED_spec()
        waves_out, fluxes_out = [], []
        for band in self.bands:
            wave_filter, transmit, eff_wave, width = self.filters[band]
            flux_interp = spectool.pyrebin.rebin_padvalue(waves, fluxes, wave_filter)
            widths = np.diff(wave_filter)
            widths = np.append(widths, widths[-1])
            flux_int = np.sum(flux_interp * transmit * widths) / np.sum(transmit * widths)
            waves_out.append(eff_wave)
            fluxes_out.append(flux_int)
        return np.array(waves_out), np.array(fluxes_out)

    def get_SED_mags(self):
        waves, fluxes = self.get_SED()
        mags = [f2m(flux, 0.0, filtname)[0] for filtname, flux in zip(self.bands, fluxes)]
        mags = np.array(mags)
        return mags

    def plot(self, ax=None, show=False):
        wave_spec, flux_spec = self.get_SED_spec()
        wave_sed, fluxe_sed = self.get_SED()
        sedmin, sedmax = np.min(fluxe_sed), np.max(fluxe_sed)
        logsedmin, logsedmax = np.log10(sedmin), np.log10(sedmax)
        deltay = (logsedmax - logsedmin) / 10
        ymin = 10 ** (logsedmin - deltay)
        ymax = 10 ** (logsedmax + deltay)
        wmin = 1200
        wmax = 30 * 10000
        arg = np.where((wave_spec > wmin) & (wave_spec < wmax))
        wave_spec = wave_spec[arg]
        flux_spec = flux_spec[arg]
        if ax is None:
            import matplotlib.pyplot as plt
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111)
            ax.plot(wave_spec, flux_spec, label='Spec')
            ax.errorbar(wave_sed, fluxe_sed, fmt='o', color='red', ms=5, label='Phot')
            ax.set_xlabel('wavelength (AA)')
            ax.set_ylabel('flux (erg/s/cm2/AA)')
            ax.legend()
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_ylim(ymin, ymax)
            if show is True:
                plt.show()
            else:
                return ax

    def _display_pars(self):
        params = [
            ("Effective Temperature (Teff)", self.teff, "K"),
            ("Surface Gravity (logg)", self.logg, "dex"),
            ("Metallicity (feh)", self.feh, "[Fe/H]"),
            ("Distance", self.distance, "pc"),
            ("Radius", self.rad, "R_sun"),
            ("Extinction Coefficient (Av)", self.Av, ""),
        ]
        
        max_length = max(len(name[0]) for name in params)
        
        print("SED Model Parameters".center(max_length + 30, "="))
        for name, value, unit in params:
            if unit:
                print(f"{name.ljust(max_length)} : {value} {unit}")
            else:
                print(f"{name.ljust(max_length)} : {value}")

        # 分块输出 bands
        print("Bands".ljust(max_length), ": ", end="")
        bands_per_line = 5
        for i in range(0, len(self.bands), bands_per_line):
            band_chunk = self.bands[i:i+bands_per_line]
            if i > 0:
                print(" " * (max_length + 3), end="")
            print(", ".join(band_chunk))
        print()
        # print("=" * (max_length + 30))

    def _display_band_data(self):
        headers = ["Band", "Wavelength", "Width", "Flux", "Mag"]
        units = ["", "AA", "AA", "erg/s/cm2/AA", ""]
        row_format = "{:<15} {:<15} {:<15} {:<15} {:<15}"
        maxlength = len(row_format.format(*headers))
        print("SED Model Band Data".center(maxlength, "="))
        print(row_format.format(*headers))
        print(row_format.format(*units))
        print("=" * maxlength)
        waves, fluxes = self.get_SED()
        mags = [f2m(flux, 0.0, filtname)[0] for filtname, flux in zip(self.bands, fluxes)]
        bands = self.bands
        waves_phot = [self.filters[band][2] for band in bands]
        widths_phot = [self.filters[band][3] for band in bands]
        for band, wave, width, flux, mag in zip(bands, waves_phot, widths_phot, fluxes, mags):
            # print('type flux:', type(flux))
            # print('flux =', flux)
            flux_str = f"{flux:.2e}"
            mag_str = f"{mag:.2f}"
            wave_str = f"{wave:.1f}"
            width_str = f"{width:.1f}"
            print(row_format.format(band, wave_str, width_str, flux_str, mag_str))
        print("=" * maxlength)

    def __str__(self):
        output_capture = io.StringIO()
        with contextlib.redirect_stdout(output_capture):
            self._display_pars()
            self._display_band_data()
        captured_output = output_capture.getvalue()
        return captured_output

    def __repr__(self):
        return self.__str__()


class ObservedSEDModel(SEDModel):
    def __init__(self, bands=None, teff=5700, logg=4.5, feh=0.0, 
                 R=1.0, distance=10.0, Av=0.0,
                 specmodel=stellarSpecModel.BTCond_Model(),
                 observed_fluxes=None, observed_errors=None,
                 observed_mags=None, observed_mag_errors=None):
        """a class to generate and compare stellar model SED with observed data

        Args:
            bands (list, optional): the band names of the SED, a example ['SDSSu', 'JV', 'SDSSr', '2MASSJ', 'W2']. Defaults to None.
            teff (float, optional): temperature of the SED. Defaults to 5700.
            logg (float, optional): surface gravity. Defaults to 4.5.
            feh (float, optional): Metallicity. Defaults to 0.0.
            R (float, optional): radius of the stellar, unit is R_sun. Defaults to 1.0.
            distance (float, optional): distance of the stellar, unit is pc. Defaults to 10.0.
            Av (float, optional): Extinction Coefficient. Defaults to 0.0.
            specmodel (stellarSpecModel, optional): stellarSpecModel used to generate the stellar spectrum. Defaults to stellarSpecModel.BTCond_Model().
            observed_fluxes (list, optional): observed fluxes for the bands. Defaults to None.
            observed_errors (list, optional): errors in the observed fluxes. Defaults to None.
            observed_mags (list, optional): observed magnitudes for the bands. Defaults to None.
            observed_mag_errors (list, optional): errors in the observed magnitudes. Defaults to None.

        Raises:
            ValueError: if specmodel is not an instance of StellarSpecModel, raise ValueError
        """
        super().__init__(bands, teff, logg, feh, R, distance, Av, specmodel)
        self.obs_mags = []
        self.obs_mag_errs = []
        self.obs_fluxes = []
        self.obs_flux_errs = []
        self.add_data(bands)

    def _add_mere_data(self, bands=None, obs_mags=None, obs_mag_errs=None, obs_fluxes=None, obs_flux_errs=None):
        obs_mags, obs_mag_errs, obs_fluxes, obs_flux_errs = self._complete_obsdata(bands, obs_mags, obs_mag_errs, obs_fluxes, obs_flux_errs)
        for mag, mag_err, flux, flux_err in zip(obs_mags, obs_mag_errs, obs_fluxes, obs_flux_errs):
            self.obs_mags.append(mag)
            self.obs_mag_errs.append(mag_err)
            self.obs_fluxes.append(flux)
            self.obs_flux_errs.append(flux_err)

    def add_data(self, bands, obs_mags=None, obs_mag_errs=None, obs_fluxes=None, obs_flux_errs=None):
        """add data to the ObservedSEDModel"""
        if isinstance(bands, str):
            self.add_band(bands)
        else:
            self.add_bands(bands)
        self._add_mere_data(bands, obs_mags, obs_mag_errs, obs_fluxes, obs_flux_errs)

    def _complete_obsdata(self, bands, obs_mags=None, obs_mag_errs=None, obs_fluxes=None, obs_flux_errs=None):
        if obs_mags is None and obs_fluxes is None:
            raise ValueError("Either obs_mags or obs_fluxes must be provided.")
        if isinstance(bands, str):
            bands = [bands, ]
            if obs_mags is not None:
                obs_mags = [obs_mags, ]
                obs_mag_errs = [obs_mag_errs, ]
            else:
                obs_fluxes = [obs_fluxes, ]
                obs_flux_errs = [obs_flux_errs, ]
        if obs_mags is not None:
            obs_fluxes, obs_flux_errs = [], []
            for mag, err, band in zip(obs_mags, obs_mag_errs, bands):
                flux, flux_err = m2f(mag, err, band)
                obs_fluxes.append(flux)
                obs_flux_errs.append(flux_err)
        else:
            obs_mags, obs_mag_errs = [], []
            for flux, err, band in zip(obs_fluxes, obs_flux_errs, bands):
                mag, mag_err = f2m(flux, err, band)
                obs_mags.append(mag)
                obs_mag_errs.append(mag_err)
        return obs_mags, obs_mag_errs, obs_fluxes, obs_flux_errs

    def get_chisq(self):
        """calculate the chi-squared value of the observed data and the model
        """
        fluxes_model = self.get_SED()[1]
        chisq = np.sum((self.obs_fluxes - fluxes_model)**2/self.obs_flux_errs**2)
        return chisq

    def _display_observations(self):
        headers = ["Band", "Wavelength", "Observed Flux", "Error", "Observed Magnitude", "Mag Error"]
        units = ["", "AA", "erg/s/cm2/AA", "erg/s/cm2/AA", "", ""]
        row_format = "{:<15} {:<15} {:<15} {:<15} {:<15} {:<15}"
        maxlength = len(row_format.format(*headers))
        print("Observed Data".center(maxlength, "="))
        print(row_format.format(*headers))
        print(row_format.format(*units))
        print("=" * maxlength)
        
        bands = self.bands
        waves_phot = [self.filters[band][2] for band in bands]
        for band, wave, obs_flux, error, obs_mag, mag_err in zip(bands, waves_phot, self.obs_fluxes, self.obs_flux_errs, self.obs_mags, self.obs_mag_errs):
            obs_flux_str = f"{obs_flux:.2e}"
            error_str = f"{error:.2e}"
            obs_mag_str =f"{obs_mag:.2f}"
            obs_mag_err_str = f"{mag_err:.2f}"
            print(row_format.format(band, wave, obs_flux_str, error_str, obs_mag_str, obs_mag_err_str))
        # print("=" * maxlength)

    def __str__(self):
        output_capture = io.StringIO()
        with contextlib.redirect_stdout(output_capture):
            self._display_pars()
            self._display_band_data()
            self._display_observations()
        captured_output = output_capture.getvalue()
        return captured_output

    def __repr__(self):
        return self.__str__()