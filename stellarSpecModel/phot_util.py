import os
import pyphot
import numpy as np
from astropy import constants as cs
from astropy import units as u


_absdir = os.path.dirname(os.path.abspath(__file__))
_ftablename = os.path.join(_absdir, 'filter_name_table.txt')
_dic_filternames = dict([i.split() for i in open(_ftablename)])
_lib = pyphot.get_library()


def load_local_filter(bandname):
    abs_finfo = os.path.join(_absdir, 'filter_data', 'filter_info.txt')
    dic_local_filters = dict([i.split() for i in open(abs_finfo)])
    if bandname in dic_local_filters:
        fname = dic_local_filters[bandname]
        abs_fname = os.path.join(_absdir, 'filter_data', fname)
        data = np.loadtxt(abs_fname)
        waves = data[:, 0]
        trans = data[:, 1]
        eff_wave = np.sum(waves*trans) / np.sum(trans)
        widths = np.diff(waves)
        widths = np.append(widths, widths[-1])
        sumflux = np.sum(trans*widths)
        width = sumflux / np.mean(trans)
        return waves, trans, eff_wave, width
    else:
        raise ValueError('filter name {} is not supported'.format(bandname))


def filtername2pyphotname(filtername):
    """convert a input filter name to pyphot filter name. For example, input 'SDSS:r' or 'SDSSr' will return 'SDSS_r'

    Args:
        filtername (str): the inputing filter name
    """
    if filtername not in _dic_filternames:
        raise ValueError('filter name {} is not supported'.format(filtername))
    return _dic_filternames[filtername]


def get_effective_wavelength(filtername):
    """get the effective wavelength of a filter

    Args:
        filtername (str): the inputing filter name
    """
    pyphot_name = filtername2pyphotname(filtername)
    return _lib[pyphot_name].leff.to('AA').value


def mag_to_flux_AB(mag, mag_err):
    """convert a input magnitude to flux

    Args:
        mag (float): the inputing magnitude
        mag_err (float): the inputing magnitude error
    """
    flux = 10**(-0.4 * (mag + 48.6))
    flux_err = flux * 0.4 * np.log(10) * mag_err
    return flux, flux_err


def get_filter_width(filtername):
    """get the filter width
    unit is AA
    """
    pyphot_name = filtername2pyphotname(filtername)
    return _lib[pyphot_name].width.to('AA').value


def convert_f_nu_to_f_lambda(fnu, wave):
    """convert a input f_nu to f_lambda, 
       unit of fnu: erg/s/cm^2/Hz
       unit of wave: AA
       output unit: erg/s/cm^2/AA
    """
    flambda = fnu * cs.c.to(u.AA/u.s).value / wave**2
    return flambda


def convert_f_lambda_to_f_nu(flambda, wave):
    """convert a input f_lambda to f_nu, 
       unit of flambda: erg/s/cm^2/AA
       unit of wave: AA
       output unit: erg/s/cm^2/Hz
    """
    fnu = flambda * wave**2 / cs.c.to(u.AA/u.s).value
    return fnu


def get_zero_flux(filtername):
    pyphot_name = filtername2pyphotname(filtername)
    return _lib[pyphot_name].Vega_zero_flux.to('erg/s/cm**2/AA').magnitude


def mag_to_flux(mag, mag_err, filtername):
    """convert a input magnitude to flux

    Args:
        mag (float): the inputing magnitude
        filtername (str): the inputing filter name
    """
    pyphot_name = filtername2pyphotname(filtername)
    if 'PS1_' in pyphot_name or 'SDSS_' in pyphot_name or 'GALEX_' in pyphot_name:
        flux, flux_err = mag_to_flux_AB(mag, mag_err)
        wave_eff= get_effective_wavelength(pyphot_name)
        flambda = convert_f_nu_to_f_lambda(flux, wave_eff)
        flambda_err = flux_err / flux * flambda
        return flambda, flambda_err
    else:
        f0 = get_zero_flux(pyphot_name)
        flux = 10 ** (-0.4 * mag) * f0
        flux_err = 0.4 * flux * np.log(10) * mag_err
        return flux, flux_err


def mags_to_fluxes(mags, mag_errs, filternames):
    if isinstance(filternames, str):
        return mag_to_flux(mags, mag_errs, filternames)
    else:
        fluxes, flux_errs = [], []
        for mag, mag_err, filtername in zip(mags, mag_errs, filternames):
            flux, flux_err = mag_to_flux(mag, mag_err, filtername)
            fluxes.append(flux)
            flux_errs.append(flux_err)
        return np.array(fluxes), np.array(flux_errs)


def flux_to_mag(flux, flux_err, filtername):
    """convert a input flux to magnitude
       unit of flux: erg/s/cm**2/AA
    """
    if filtername not in _lib.content:
        return np.nan, np.nan
    pyphot_name = filtername2pyphotname(filtername)
    if 'PS1_' in pyphot_name or 'SDSS_' in pyphot_name or 'GALEX_' in pyphot_name:
        wave_eff= get_effective_wavelength(pyphot_name)
        fnu = convert_f_lambda_to_f_nu(flux, wave_eff)
        mag = -48.6 - 2.5 * np.log10(fnu)
        mag_err = 2.5 * flux_err / (np.log(10) * flux)
        return mag, mag_err
    else:
        f0 = get_zero_flux(pyphot_name)
        mag = -2.5 * np.log10(flux / f0)
        mag_err = 2.5 * flux_err / (np.log(10) * flux)
        return mag, mag_err


def fluxes_to_mags(fluxes, flux_errs, filternames):
    if isinstance(filternames, str):
        return flux_to_mag(fluxes, flux_errs, filternames)
    else:
        mags, mag_errs = [], []
        for flux, flux_err, filtername in zip(fluxes, flux_errs, filternames):
            mag, mag_err = flux_to_mag(flux, flux_err, filtername)
            mags.append(mag)
            mag_errs.append(mag_err)
        return np.array(mags), np.array(mag_errs)