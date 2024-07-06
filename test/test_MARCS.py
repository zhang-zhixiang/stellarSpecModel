from stellarSpecModel import MARCS_Model
import matplotlib.pyplot as plt
import spectool


def test_MARCS_Model():
    model = MARCS_Model()
    print('teff grid =', model.teff_grid)
    print('logg grid =', model.logg_grid)
    print('feh grid =', model.feh_grid)
    teff = 3750
    logg = 4.5
    feh = 0.0
    print('teff =', teff, 'logg =', logg, 'feh =', feh)
    wave = model.wavelength
    flux = model.get_flux(teff, feh, logg)
    nflux = spectool.spec_filter.gaussian_filter(wave, flux, 120)
    plt.plot(wave, flux, label='Teff = %d, logg = %.1f, [Fe/H] = %.1f' % (teff, logg, feh))
    plt.plot(wave, nflux, color='red')
    plt.legend()
    plt.xlabel('Wavelength (Angstrom)')
    plt.ylabel('Flux')
    plt.title('MARCS Model')
    plt.show()


if __name__ == '__main__':
    test_MARCS_Model()