# stellarSpecModel

This is a simple package to interpolate the stellar spectral grid. Users provide stellar parameters (Teff, FeH, logg), the package will return the corresponding stellar spectrum.

## Installation

You can install the package by using pip:
```bash
pip install git+https://github.com/zzxihep/stellarSpecModel.git
```
Some system may raise an error `error: externally-managed-environment`, you can try to add `--break-system-packages` to the command above.
```bash
pip install git+https://github.com/zzxihep/stellarSpecModel.git --break-system-packages
```
or clone the repository and install it:
```bash
git clone https://github.com/zzxihep/stellarSpecModel.git
cd stellarSpecModel
python setup.py install --user
```

## Dependencies

The package depends on the following packages:
- numpy
- scipy
- astropy
- [selenium](https://github.com/SeleniumHQ/selenium) (optional, for downloading the grid)

The package need to download the stellar spectral grid from website, so you need to install [selenium](https://github.com/SeleniumHQ/selenium) to download the grids. Otherwise, you can download the grids manually and put them in the `~/.stellarSpecModel/grid_data/` folder. The URL of the grid files are:
- [https://www.jianguoyun.com/p/DZmcNoUQ2ZfcCBjW-5cFIAA](https://www.jianguoyun.com/p/DZmcNoUQ2ZfcCBjW-5cFIAA)
- [https://www.jianguoyun.com/p/Def5fgsQ2ZfcCBiYmpgFIAA](https://www.jianguoyun.com/p/Def5fgsQ2ZfcCBiYmpgFIAA)

see more details in config.py

## Usage

The package currently provides two classes: `MARCS_Model` and `BTCond_Model` for the MARCS and BT-Cond grid respectively. The usage of the two classes are the same. MARCS grid has a wavelength coverage of 3800-9700 Angstrom with high sampling, while BT-Cond grid has a wavelength coverage of 300-3E5 Angstrom with low sampling. You can also use your own grid by inheriting the `StellarSpecModel` class.

Here we show an example of using the package:
```python
from stellarSpecModel import MARCS_Model, BTCond_Model
import matplotlib.pyplot as plt

btcond_model = BTCond_Model()
marcs_model = MARCS_Model()
print('MARCS wavelength unit =', marcs_model.wavelength_units)
print('MARCS flux units =', marcs_model.flux_units)
print('min and max Teff of MARCS grid =', marcs_model.min_teff, marcs_model.max_teff)
print('min and max FeH of MARCS grid =', marcs_model.min_feh, marcs_model.max_feh)
print('min and max logg of MARCS grid =', marcs_model.min_logg, marcs_model.max_logg)
print('Teff grid of MARCS =', marcs_model.teff_grid)
print('FeH grid of MARCS =', marcs_model.feh_grid)
print('logg grid of MARCS =', marcs_model.logg_grid)
print('\n')

print('BTCond wavelength unit =', btcond_model.wavelength_units)
print('BTCond flux units =', btcond_model.flux_units)
print('min and max Teff of BTCond grid =', btcond_model.min_teff, btcond_model.max_teff)
print('min and max FeH of BTCond grid =', btcond_model.min_feh, btcond_model.max_feh)
print('min and max logg of BTCond grid =', btcond_model.min_logg, btcond_model.max_logg)
print('Teff grid of BTCond =', btcond_model.teff_grid)
print('FeH grid of BTCond =', btcond_model.feh_grid)
print('logg grid of BTCond =', btcond_model.logg_grid)

teff = 5700
feh = 0
logg = 4.5
marcs_flux = marcs_model.get_flux(teff, feh, logg)
btcond_flux = btcond_model.get_flux(teff, feh, logg)
plt.plot(marcs_model.wavelength, marcs_flux, label='MARCS')
plt.plot(btcond_model.wavelength, btcond_flux, label='BTCond')
plt.legend()
plt.xlabel('Wavelength ({})'.format(marcs_model.wavelength_units))
plt.ylabel('Flux ({})'.format(marcs_model.flux_units))
plt.xlim(2750, 10000)
plt.ylim(0, 1.1 * max(marcs_flux))
plt.show()
```
The output of the above code is:
```bash
MARCS wavelength unit = Angstrom
MARCS flux units = erg / (Angstrom cm2 s)
min and max Teff of MARCS grid = 2700.0 8000.0
min and max FeH of MARCS grid = -4.0 1.0
min and max logg of MARCS grid = 2.0 5.0
Teff grid of MARCS = [2700. 2800. 2900. 3000. 3100. 3200. 3300. 3400. 3500. 3600. 3700. 3800.
 3900. 4000. 4250. 4500. 4750. 5000. 5250. 5500. 5750. 6000. 6250. 6500.
 6750. 7000. 7250. 7500. 7750. 8000.]
FeH grid of MARCS = [-4.   -3.   -2.5  -2.   -1.5  -1.   -0.75 -0.5  -0.25  0.    0.25  0.5
  0.75  1.  ]
logg grid of MARCS = [2.  2.5 3.  3.5 4.  4.5 5. ]


BTCond wavelength unit = Angstrom
BTCond flux units = erg / (Angstrom cm2 s)
min and max Teff of BTCond grid = 3500.0 8000.0
min and max FeH of BTCond grid = -4.0 0.5
min and max logg of BTCond grid = 1.0 5.0
Teff grid of BTCond = [3500. 3600. 3700. 3800. 3900. 4000. 4100. 4200. 4300. 4400. 4500. 4600.
 4700. 4800. 4900. 5000. 5100. 5200. 5300. 5400. 5500. 5600. 5700. 5800.
 5900. 6000. 6100. 6200. 6300. 6400. 6500. 6600. 6700. 6800. 6900. 7000.
 7200. 7400. 7600. 7800. 8000.]
FeH grid of BTCond = [-4.  -3.5 -3.  -2.5 -2.  -1.5 -1.  -0.5  0.   0.3  0.5]
logg grid of BTCond = [1.  1.5 2.  2.5 3.  3.5 4.  4.5 5. ]
```
![example](https://github.com/zzxihep/stellarSpecModel/blob/main/example.png)

## Note

The flux of the grid spectrum corresponds to the flux received by an observer located at the surface of a star. Therefore, assuming that the radius of a star is equal to $1 R_\odot$, and the observer's position is at a distance of 100 parsecs, the received flux will be `nflux = flux * (cs.R_sun / (100*u.pc)).to('').value**2`.