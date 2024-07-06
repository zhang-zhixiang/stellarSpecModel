import stellarSpecModel
import matplotlib.pyplot as plt


def test_stellarSpecModel():
    model = stellarSpecModel.BTCond_Model()
    model_hiRes = stellarSpecModel.BTCond_Model_hiRes()
    model_R7500 = stellarSpecModel.BTCond_Model_R7500()
    model_R1800 = stellarSpecModel.BTCond_Model_R1800()
    model_R500 = stellarSpecModel.BTCond_Model_R500()
    model_R100 = stellarSpecModel.BTCond_Model_R100()

    teff = 5700
    logg = 4.5
    feh = 0.0
    flux_R100 = model_R100.get_flux(teff, feh, logg)
    flux_R500 = model_R500.get_flux(teff, feh, logg)
    flux_R1800 = model_R1800.get_flux(teff, feh, logg)
    flux_R7500 = model_R7500.get_flux(teff, feh, logg)
    flux_hiRes = model_hiRes.get_flux(teff, feh, logg)
    flux = model.get_flux(teff, feh, logg)

    plt.plot(model_hiRes.wavelength, flux_hiRes, label='hiRes')
    plt.plot(model.wavelength, flux, label='BTCond')
    plt.plot(model_R7500.wavelength, flux_R7500, label='R=7500')
    plt.plot(model_R1800.wavelength, flux_R1800, label='R=1800')
    plt.plot(model_R500.wavelength, flux_R500, label='R=500')
    plt.plot(model_R100.wavelength, flux_R100, label='R=100')
    plt.legend()

    plt.xscale('log')
    plt.yscale('log')

    plt.show()


if __name__ == '__main__':
    test_stellarSpecModel()