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


if __name__ == '__main__':
    test_SED_model()