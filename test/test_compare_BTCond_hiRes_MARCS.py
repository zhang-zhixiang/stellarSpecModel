import stellarSpecModel
import matplotlib.pyplot as plt
import numpy as np


def main():
    btmodel = stellarSpecModel.BTCond_Model_hiRes()
    marcsmodel = stellarSpecModel.MARCS_Model_hiRes()

    print('BT model high')
    print(np.diff(btmodel.wavelength))
    print('MARCS model high')
    print(np.diff(marcsmodel.wavelength))

    btmodel_low = stellarSpecModel.BTCond_Model()
    marcsmodel_low = stellarSpecModel.MARCS_Model()

    print('BT model low')
    print(np.diff(btmodel_low.wavelength))
    print('MARCS model low')
    print(np.diff(marcsmodel_low.wavelength))

    plt.plot(btmodel.wavelength[1:], np.diff(btmodel.wavelength), label='BT-Cond high')

    teff = 5777
    logg = 4.44
    feh = 0.1

    flux_bt = btmodel.get_flux(teff, feh+0.3, logg)
    flux_marcs = marcsmodel.get_flux(teff, feh, logg)

    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)
    ax.plot(marcsmodel.wavelength, flux_marcs, label='MARCS')
    ax.plot(btmodel.wavelength, flux_bt, label='BT-Cond')
    ax.set_xlabel(f'Wavelength ({btmodel.wavelength_units})')
    ax.set_ylabel(f'Flux ({btmodel.flux_units})')
    ax.legend()
    plt.show()


if __name__ == '__main__':
    main()
