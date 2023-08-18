from stellarSpecModel import MARCS_Model, BTCond_Model
import matplotlib.pyplot as plt


def main():
    btcond_model = BTCond_Model()
    marcs_model = MARCS_Model()
    print('MARCS wavelength unit =', marcs_model.wavelength_units)
    print('MARCS flux units =', marcs_model.flux_units)
    print('min and max Teff of MARCS grid =', marcs_model.min_teff, marcs_model.max_teff)
    print('min and max FeH of MARCS grid =', marcs_model.min_feh, marcs_model.max_feh)
    print('min and max logg of MARCS grid =', marcs_model.min_logg, marcs_model.max_logg)

    print('BTCond wavelength unit =', btcond_model.wavelength_units)
    print('BTCond flux units =', btcond_model.flux_units)
    print('min and max Teff of BTCond grid =', btcond_model.min_teff, btcond_model.max_teff)
    print('min and max FeH of BTCond grid =', btcond_model.min_feh, btcond_model.max_feh)
    print('min and max logg of BTCond grid =', btcond_model.min_logg, btcond_model.max_logg)

    teff = 5000
    feh = 0
    logg = 4.5
    marcs_flux = marcs_model.get_flux(teff, feh, logg)
    btcond_flux = btcond_model.get_flux(teff, feh, logg)
    plt.plot(marcs_model.wavelength, marcs_flux, label='MARCS')
    plt.plot(btcond_model.wavelength, btcond_flux, label='BTCond')
    plt.legend()
    plt.xlabel('Wavelength ({})'.format(marcs_model.wavelength_units))
    plt.ylabel('Flux ({})'.format(marcs_model.flux_units))
    plt.show()


if __name__ == '__main__':
    main()
