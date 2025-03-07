import matplotlib.pyplot as plt
import pyphot
import numpy as np


def read_trmit():
    fname = '../stellarSpecModel/filter_data/WISE_WISE.W3.dat'
    data = np.loadtxt(fname)
    wave = data[:, 0]
    trmit = data[:, 1]
    return wave, trmit


def main():
    lib = pyphot.get_library()
    fw3 = lib["WISE_RSR_W3"]
    wave = fw3.wavelength.to('AA').value
    trans = fw3.transmit
    plt.plot(wave, trans/np.max(trans), label='pyphot')

    wave, transmit = read_trmit()
    plt.plot(wave, transmit/np.max(transmit), label='SVO')
    plt.legend()

    plt.show()


def read_trmit_w4():
    fname = '../stellarSpecModel/filter_data/WISE_WISE.W4.dat'
    data = np.loadtxt(fname)
    wave = data[:, 0]
    trmit = data[:, 1]
    return wave, trmit


def main2():
    lib = pyphot.get_library()
    fw4 = lib["WISE_RSR_W4"]
    wave = fw4.wavelength.to('AA').value
    trans = fw4.transmit
    plt.plot(wave, trans/np.max(trans), label='pyphot')

    wave, transmit = read_trmit_w4()
    plt.plot(wave, transmit/np.max(transmit), label='AAVSO')
    plt.legend()

    plt.show()


if __name__ == '__main__':
    main()
    main2()