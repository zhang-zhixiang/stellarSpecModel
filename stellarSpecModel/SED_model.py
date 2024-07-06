from astropy import constants as cs
from extinction import fitzpatrick99
from extinction import apply
import pyphot
import spectool
import numpy as np
from . import stellarSpecModel
from .phot_util import flux_to_mag as f2m
from .phot_util import filtername2pyphotname


class SEDModel:
    def __init__(self, bands=None, specmodel=stellarSpecModel.BTCond_Model()):
        if isinstance(specmodel, stellarSpecModel.StellarSpecModel):
            self.stellar_model = specmodel
        else:
            raise ValueError('specmodel must be an instance of StellarSpecModel')
        self.bands = []
        self.filters = {}
        self.pyphot_lib = pyphot.get_library()
        self._rat_rsun_pc = cs.R_sun.to('pc').value
        self.Rv = 3.1
        self.teff = 5700
        self.logg = 4.5
        self.feh = 0.0
        self.distance = 10.0
        self.Av = 0.0
        self.rad = 1.0
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

    def get_SED_spec(self, teff, logg, feh, R, distance, Av=0.0):
        waves = self.stellar_model.wavelength
        fluxes = self.stellar_model.get_flux(teff, feh, logg)
        rat = (R / distance * self._rat_rsun_pc) ** 2
        fluxes *= rat
        ext = fitzpatrick99(waves, Av, self.Rv)
        nfluxes = apply(ext, fluxes)
        return waves, nfluxes

    def get_SED(self, teff, logg, feh, R, distance, Av=0.0):
        waves, fluxes = self.get_SED_spec(teff, logg, feh, R, distance, Av)
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

    def get_SED_mags(self, teff, logg, feh, R, distance, Av=0.0):
        waves, fluxes = self.get_SED(teff, logg, feh, R, distance, Av)
        mags = [f2m(flux, 0.0, filtname)[0] for filtname, flux in zip(self.bands, fluxes)]
        mags = np.array(mags)
        return mags

    def plot(self, teff, logg, feh, R, distance, Av=0.0, ax=None, show=False):
        wave_spec, fluxe_spec = self.get_SED_spec(teff, logg, feh, R, distance, Av)
        wave_sed, fluxe_sed = self.get_SED(teff, logg, feh, R, distance, Av)
        if ax is None:
            import matplotlib.pyplot as plt
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111)
            ax.plot(wave_spec, fluxe_spec, label='spec')
            ax.errorbar(wave_sed, fluxe_sed, fmt='o', color='red', ms=5, label='sed')
            ax.set_xlabel('wavelength (AA)')
            ax.set_ylabel('flux (erg/s/cm2/AA)')
            ax.legend()
            ax.set_xscale('log')
            ax.set_yscale('log')
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
            ("Extinction Coefficient (Av)", self.av, ""),
            ("Bands", ', '.join(self.bands), "")
        ]
        max_length = max([len(name[0]) for name in params])
        print("SED Model Parameters".center(max_length + 30, "="))
        for name, value, unit in params:
            if unit:
                print(f"{name.ljust(max_length)} : {value} {unit}")
            else:
                print(f"{name.ljust(max_length)} : {value}")
        print("=" * (max_length + 30))

    def _display_band_data(self):
        headers = ["Band", "Wavelength", "Width", "Flux", "Mag"]
        units = ["", "AA", "AA", "erg/s/cm2/AA", ""]
        row_format = "{:<10} {:<15} {:<15} {:<15} {:<15}"
        maxlength = len(row_format.format(*headers)) + 10
        print("SED Model Band Data".center(maxlength, "="))
        print(row_format.format(*headers))
        print(row_format.format(*units))
        print("=" * maxlength)
        fluxes = self.get_SED_mags(self.teff, self.logg, self.feh, self.rad, self.distance, self.Av)
        mags = [f2m(flux, 0.0, filtname)[0] for filtname, flux in zip(self.bands, fluxes)]
        bands = self.bands
        waves = [self.filters[i][2] for i in len(bands)]
        widths = [self.filters[i][3] for i in len(bands)]
        for band, wave, width, flux, mag in zip(bands, waves, widths, fluxes, mags):
            flux_str = f"{flux:.2e}"
            mag_str = f"{mag:.2f}"
            print(row_format.format(band, wave, width, flux_str, mag_str))
        print("=" * maxlength)

    def __str__(self):
        self._display_pars()
        self._display_band_data()
