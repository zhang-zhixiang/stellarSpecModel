from stellarSpecModel import TlustyModel
import numpy as np
import matplotlib.pyplot as plt


def main():
    model = TlustyModel()
    print('teff grid:', model.teff_grid)
    print('logg grid:', model.logg_grid)
    print('metal grid:', model.feh_grid)
    wave = model.wavelength
    teff = 16500
    logg = 2.3
    metal = 0.4
    flux = model.get_flux(teff, metal, logg)
    plt.plot(wave, flux)
    plt.show()


if __name__ == '__main__':
    main()