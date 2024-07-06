import matplotlib.pyplot as plt
from stellarSpecModel import TlustyWDModel


def plot_model():
    model = TlustyWDModel()

    print(model.teff_grid)
    print(model.logg_grid)

    teff = 20000
    logg = 8.0
    wave = model.wavelength
    flux1 = model.get_flux(teff, logg)

    teff2 = 20000
    logg2 = 8.5
    flux2 = model.get_flux(teff2, logg2)

    teff3 = 16000
    logg3 = 8.0
    flux3 = model.get_flux(teff3, logg3)

    plt.plot(wave, flux1, label=f'Teff={teff}, logg={logg}')
    plt.plot(wave, flux2, label=f'Teff={teff2}, logg={logg2}')
    plt.plot(wave, flux3, label=f'Teff={teff3}, logg={logg3}')
    plt.legend()
    plt.xlabel('Wavelength (AA)')
    plt.ylabel('Flux')
    plt.show()


if __name__ == '__main__':
    plot_model()