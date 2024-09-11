from stellarSpecModel import BinarySEDModel
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def read_SEDinfo():
    csvname = 'SEDinfo.csv'
    df = pd.read_csv(csvname)
    nufnu = df['nufnu'].values
    nufnuerr = df['nufnuerr'].values
    # arg = np.isfinite(nufnuerr)
    arg_nan = np.isnan(nufnuerr)
    nufnuerr[arg_nan] = 0.1 * nufnu[arg_nan]
    df['nufnuerr'] = nufnuerr
    arg = np.isfinite(nufnu)
    df = df[arg]
    return df


def test():
    teff1 = 6171
    teff2 = 5122
    R1 = 1.8
    R2 = 2.93
    D = 169.89
    Av = 0.0
    feh = -0.26
    df = pd.read_csv('SEDinfo.csv')
    filternames = df['filtname'].values
    nufnus = df['nufnu'].values
    nufnu_errs = df['nufnuerr'].values
    binary_model = BinarySEDModel(teff1=teff1, teff2=teff2, R1=R1, R2=R2, D=D, Av=Av, feh1=feh)
    binary_model.add_data(filternames, obs_nufnus=nufnus, obs_nufnuerrs=nufnu_errs)
    binary_model.plot(show=True)

    syserr = 0.1
    binary_model.set_pars(syserr=syserr)
    print('binarymodel.syserr =', binary_model.syserr)
    print('binarymodel.obs_fluxes =', binary_model._obs_fluxes)
    print('binarymodel.obs_fluxerrs =', binary_model._obs_fluxerrs)

    chisq = binary_model.get_chisq()
    chisq_sys = binary_model.get_chisq_syserr()
    lnlike = binary_model.get_lnlike()

    print('chisq = ', chisq)
    print('chisq_sys = ', chisq_sys)
    print('lnlike = ', lnlike)

    syserrs = np.logspace(-3, -0.4, 50)
    chisq_syss = []
    lnlikes = []
    for syserr in syserrs:
        binary_model.set_pars(syserr=syserr)
        chisq_syss.append(binary_model.get_chisq_syserr())
        lnlikes.append(binary_model.get_lnlike())
    plt.plot(syserrs, chisq_syss)
    plt.xlabel('syserr')
    plt.ylabel('chisq_sys')
    plt.show()

    plt.plot(syserrs, lnlikes)
    plt.xlabel('syserr')
    plt.ylabel('lnlike')
    plt.show()


if __name__ == '__main__':
    test()
