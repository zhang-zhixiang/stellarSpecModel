import pyphot
import numpy as np
from astropy import constants as cs
import matplotlib.pyplot as plt
from extinction import apply
from extinction import fitzpatrick99
import spectool
from . import stellarSpecModel
from .phot_util import fluxes_to_mags as f2ms
from .phot_util import mags_to_fluxes as m2fs
from .phot_util import filtername2pyphotname as f2p
from .phot_util import load_local_filter


class BinarySEDModel:
    def __init__(self, teff1=None, feh1=None, logg1=None, R1=None, 
                 D=None, Av=0.0, teff2=None, feh2=None, logg2=None, R2=None, 
                 syserr=None, specmodel=stellarSpecModel.BTCond_Model()):
        self.teff1 = teff1
        self.feh1 = feh1
        self.logg1 = logg1
        self.R1 = R1
        self.D = D
        self.Av = Av
        self.teff2 = teff2
        self.feh2 = feh2
        self.logg2 = logg2
        self.R2 = R2
        self.stellar_model = specmodel
        self.syserr = syserr

        self.filters = {}
        self._bands = []
        self._eff_waves_SED = []
        self._widths_band = []
        self._pyphot_lib = pyphot.get_library()
        self._rat_rsun_pc = cs.R_sun.to('pc').value
        self._Rv = 3.1
        self._obs_mags = []
        self._obs_magerrs = []
        self._obs_fluxes = []
        self._obs_fluxerrs = []

    def set_pars(self, teff1=None, feh1=None, logg1=None, R1=None, 
                 D=None, Av=None, teff2=None, feh2=None, logg2=None, R2=None,
                 syserr=None):
        self.teff1 = teff1 if teff1 is not None else self.teff1
        self.feh1 = feh1 if feh1 is not None else self.feh1
        self.logg1 = logg1 if logg1 is not None else self.logg1
        self.R1 = R1 if R1 is not None else self.R1
        self.D = D if D is not None else self.D
        self.Av = Av if Av is not None else self.Av
        self.teff2 = teff2 if teff2 is not None else self.teff2
        if feh2 is not None:
            self.feh2 = feh2
        elif feh1 is not None:
            self.feh2 = feh1
        self.logg2 = logg2 if logg2 is not None else self.logg2
        self.R2 = R2 if R2 is not None else self.R2
        self.syserr = syserr if syserr is not None else self.syserr

    def _load_filter(self, bandname):
        try:
            tfilter = self._pyphot_lib[bandname]
            waves = tfilter.wavelength.to('AA').value
            trans = tfilter.transmit
            eff_wave = tfilter.leff.to('AA').value
            width = tfilter.width.to('AA').value
            return waves, trans, eff_wave, width
        except Exception:
            return load_local_filter(bandname)

    def add_data(self, bands, obs_mags=None, obs_magerrs=None, 
                 obs_fluxes=None, obs_fluxerrs=None, 
                 obs_nufnus=None, obs_nufnuerrs=None):
        """add photometric data points to the model. You can add only bands but no observed data or bands with observed data. The observed data can be in magnitudes or fluxes. If you add only bands, you can get the model SED, but you cannot obtain chisq or likelihood from the model.

        If you set obs_mags, obs_fluxes will be ignored.

        Args:
            bands (list): band names of the photometric data points.
            obs_mags (list or None, optional): observed mags of the photometric data. Defaults to None.
            obs_magerrs (list or None, optional): obs mag errs. Defaults to None.
            obs_fluxes (list or None, optional): obs fluxes, in units of erg/s/cm^2/AA. Defaults to None.
            obs_fluxerrs (list or None, optional): obs flux errs. Defaults to None.
            obs_nufnus (list or None, optional): observed nufnus, in units of erg/s/cm^2. Defaults to None.
            obs_nufnuerrs (list or None, optional): obs nufnu errs. Defaults to None.
        """
        std_bands = [f2p(b) for b in bands]
        eff_waves = []
        for bandname in std_bands:
            waves, trans, eff_wave, width = self._load_filter(bandname)
            self.filters[bandname] = [waves, trans]
            eff_waves.append(eff_wave)
            self._widths_band.append(width)
        if obs_mags is not None:
            obs_fluxes, obs_fluxerrs = m2fs(obs_mags, obs_magerrs, std_bands)
        elif obs_fluxes is not None:
            obs_mags, obs_magerrs = f2ms(obs_fluxes, obs_fluxerrs, std_bands)
        elif obs_nufnus is not None:
            obs_fluxes = np.array(obs_nufnus) / np.array(eff_waves)
            obs_fluxerrs = np.array(obs_nufnuerrs) / np.array(eff_waves)
            obs_mags, obs_magerrs = f2ms(obs_fluxes, obs_fluxerrs, std_bands)
        else:
            obs_mags = [np.nan,] * len(bands)
            obs_magerrs = [np.nan,] * len(bands)
            obs_fluxes = [np.nan,] * len(bands)
            obs_fluxerrs = [np.nan,] * len(bands)
        self._bands += list(std_bands)
        self._obs_mags += list(obs_mags)
        self._obs_magerrs += list(obs_magerrs)
        self._obs_fluxes += list(obs_fluxes)
        self._obs_fluxerrs += list(obs_fluxerrs)
        self._eff_waves_SED += list(eff_waves)

    def get_SED_spec1(self):
        teff = self.teff1
        feh = self.feh1
        if feh is None:
            feh = 0.0
        logg = self.logg1
        if logg is None:
            logg = 4.4
        waves = self.stellar_model.wavelength
        fluxes = self.stellar_model.get_flux(teff, feh, logg)
        rat = (self.R1 / self.D * self._rat_rsun_pc) ** 2
        fluxes *= rat
        ext = fitzpatrick99(waves, self.Av, self._Rv)
        nfluxes = apply(ext, fluxes)
        return waves, nfluxes

    def get_SED_spec2(self):
        teff = self.teff2
        feh = self.feh2
        if feh is None:
            if self.feh1 is not None:
                feh = self.feh1
            else:
                feh = 0.0
        logg = self.logg2
        if logg is None:
            logg = 4.4
        waves = self.stellar_model.wavelength
        fluxes = self.stellar_model.get_flux(teff, feh, logg)
        rat = (self.R2 / self.D * self._rat_rsun_pc) ** 2
        fluxes *= rat
        ext = fitzpatrick99(waves, self.Av, self._Rv)
        nfluxes = apply(ext, fluxes)
        return waves, nfluxes

    def get_SED_spec(self):
        wave_spec1, spec1 = self.get_SED_spec1()
        wave_spec2, spec2 = self.get_SED_spec2()
        return wave_spec1, spec1 + spec2

    def _SED_from_spec(self, wave_spec, spec):
        SED_outs = []
        for ind, band in enumerate(self._bands):
            wave_filter, trans = self.filters[band]
            flux_interp = spectool.pyrebin.rebin_padvalue(wave_spec, spec, wave_filter)
            widths = np.diff(wave_filter)
            widths = np.append(widths, widths[-1])
            flux_int = np.sum(flux_interp * trans * widths) / np.sum(trans * widths)
            SED_outs.append(flux_int)
        return SED_outs

    def get_SED1(self):
        wave_spec, spec = self.get_SED_spec1()
        wave_SED = np.array(self._eff_waves_SED)
        SED_outs = np.array(self._SED_from_spec(wave_spec, spec))
        return wave_SED, SED_outs

    def get_SED2(self):
        wave_spec, spec = self.get_SED_spec2()
        wave_SED = np.array(self._eff_waves_SED)
        SED_outs = np.array(self._SED_from_spec(wave_spec, spec))
        return wave_SED, SED_outs

    def get_SED(self):
        wave_spec1, spec1 = self.get_SED_spec1()
        wave_spec2, spec2 = self.get_SED_spec2()
        spec = spec1 + spec2
        wave_SED = np.array(self._eff_waves_SED)
        SED_outs = np.array(self._SED_from_spec(wave_spec1, spec))
        return wave_SED, SED_outs

    def get_SED_mags1(self):
        wave_SED, SED_outs = self.get_SED1()
        errs = [0, ] * len(self._bands)
        mags, mag_errs = f2ms(SED_outs, errs, self._bands)
        return self._bands, mags

    def get_SED_mags2(self):
        wave_SED, SED_outs = self.get_SED2()
        errs = [0, ] * len(self._bands)
        mags, mag_errs = f2ms(SED_outs, errs, self._bands)
        return self._bands, mags

    def get_SED_mags(self):
        wave_SED, SED_outs = self.get_SED()
        errs = [0, ] * len(self._bands)
        mags, mag_errs = f2ms(SED_outs, errs, self._bands)
        return self._bands, mags

    def get_chisq(self):
        wave_SED, SED_model = self.get_SED()
        fluxes_obs = np.array(self._obs_fluxes)
        flux_errs_obs = np.array(self._obs_fluxerrs)
        chisq = np.sum((fluxes_obs - SED_model) ** 2 / flux_errs_obs ** 2)
        return chisq

    def get_chisq_syserr(self):
        wave_SED, SED_model = self.get_SED()
        fluxes_obs = np.array(self._obs_fluxes)
        flux_errs_obs = np.array(self._obs_fluxerrs)
        arg = np.isnan(flux_errs_obs)
        flux_errs_obs[arg] = 0
        if self.syserr is not None:
            syserrs = fluxes_obs * self.syserr
            sigma = np.sqrt(flux_errs_obs ** 2 + syserrs ** 2)
        else:
            sigma = flux_errs_obs
        chisq = np.sum((fluxes_obs - SED_model) ** 2 / sigma ** 2)
        return chisq

    def get_lnlike(self):
        wave_SED, SED_model = self.get_SED()
        fluxes_obs = np.array(self._obs_fluxes)
        flux_errs_obs = np.array(self._obs_fluxerrs)
        arg = np.isnan(flux_errs_obs)
        flux_errs_obs[arg] = 0
        if self.syserr is not None:
            syserrs = fluxes_obs * self.syserr
            sigma = np.sqrt(flux_errs_obs ** 2 + syserrs ** 2)
        else:
            sigma = flux_errs_obs
        lnlike = -0.5 * np.sum((fluxes_obs - SED_model)**2 / sigma**2 + np.log(sigma**2))
        return lnlike

    def plot(self, ax=None, show=False):
        wave_spec1, flux_spec1 = self.get_SED_spec1()
        wave_spec2, flux_spec2 = self.get_SED_spec2()
        wave_spec, flux_spec = self.get_SED_spec()

        wave_SED1, SED1_model = self.get_SED1()
        _, SED2_model = self.get_SED2()
        _, SED_model = self.get_SED()

        fluxes_obs = np.array(self._obs_fluxes)
        fluxerrs_obs = np.array(self._obs_fluxerrs)

        sedmin, sedmax = np.min(SED_model), np.max(SED_model)
        logsedmin, logsedmax = np.log10(sedmin), np.log10(sedmax)
        deltay = (logsedmax - logsedmin) / 10
        ymin = 10 ** (logsedmin - deltay)
        ymax = 10 ** (logsedmax + deltay)

        wmin, wmax = np.min(wave_SED1), np.max(wave_SED1)
        logwmin, logwavemax = np.log10(wmin), np.log10(wmax)
        deltax = (logwavemax - logwmin) / 10
        xmin = 10 ** (logwmin - deltax)
        xmax = 10 ** (logwavemax + deltax)

        if ax is None:
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111)
        
        ax.plot(wave_spec1, flux_spec1, label='Spec1')
        ax.plot(wave_spec2, flux_spec2, label='Spec2')
        ax.plot(wave_spec, flux_spec, label='Spec1 + Spec2')

        # ax.errorbar(wave_SED1, SED1_model, fmt='o', ms=5, label='SED model 1')
        # ax.errorbar(wave_SED1, SED2_model, fmt='o', ms=5, label='SED model 2')
        ax.errorbar(wave_SED1, SED_model, fmt='o', ms=5, label='Combined SED model')

        # print(fluxes_obs)
        arg = np.isfinite(fluxes_obs)
        wave_SED1 = wave_SED1[arg]
        fluxes_obs = fluxes_obs[arg]
        fluxerrs_obs = fluxerrs_obs[arg]
        if len(wave_SED1) > 0:
            ax.errorbar(wave_SED1, fluxes_obs, yerr=fluxerrs_obs, fmt='s', label='Observed data', markersize=5, color='k')

        ax.set_xlabel('wavelength (AA)')
        ax.set_ylabel('flux (erg/s/cm2/AA)')

        ax.legend()
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

        if show:
            plt.show()
            return None
        else:
            return ax
