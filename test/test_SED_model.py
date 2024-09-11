import stellarSpecModel


def test_SED_model():
    teff = 5700
    logg = 4.5
    feh = 0.0
    rad = 1.0
    dist = 3000
    av = 0.4
    bands = [
        'SDSSu',
        'SDSSg',
        'SDSSr',
        'SDSSi',
        'SDSSz',
        '2MASSJ',
        '2MASSH',
        '2MASSKs',
        'W1',
        'W2',
        'W3',
        'W4',
    ]
    model = stellarSpecModel.SEDModel(bands)
    model.set_SED_pars(teff, logg, feh, rad, dist, av)
    print(model)
    model.plot(show=True)


def read_data(fname):
    lst = open(fname).readlines()
    lst = [i.split() for i in lst if i[0] != '#']
    filternames = [i[0] for i in lst]
    mags = [float(i[1]) for i in lst]
    mag_errs = [float(i[2]) for i in lst]
    return filternames, mags, mag_errs


def test_SED_observed_model():
    fname = 'mags.dat'
    filternames, mags, mag_errs = read_data(fname)
    teff = 13600
    feh = 0.05
    logg = 2.0
    rad = 15.7
    dist = 3726
    Av = 2.94
    model = stellarSpecModel.SED_model.ObservedSEDModel(bands=filternames,
        teff=teff, logg=logg, feh=feh, R=rad, distance=dist, Av=Av,
        specmodel=stellarSpecModel.BTCond_Model_R100(),
        observed_mags=mags, observed_mag_errors=mag_errs
        )
    print(model)
    model.plot(show=True)


if __name__ == '__main__':
    test_SED_model()
    test_SED_observed_model()