# StellarSpecModel Documentation

## Overview

`StellarSpecModel` is a Python package to interpolate the stellar spectral grid. Users provide stellar parameters (Teff, FeH, logg), the package will return the corresponding stellar spectrum.

This packagge also designed for generating and analyzing theoretical stellar spectral energy distributions (SEDs). The package includes functionality for both single and binary star systems, incorporating extinction models and the ability to handle photometric data in various filter bands.

### Key Modules:

1. **SED_model**: This module is responsible for generating theoretical SEDs for a single star, given its parameters like effective temperature (`Teff`), metallicity (`FeH`), surface gravity (`logg`), radius of the star, distance of the star to observer, and extinction parameter. The model generates the flux across a range of wavelengths.

2. **binary_SED_model**: This module extends the capabilities of the `SED_model` to handle binary star systems. It allows for the creation of composite SEDs from two stars, each with its own set of parameters. The model also includes the ability to apply extinction and calculate the resulting observed fluxes for each component of the binary system.

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

## Detailed Descriptions of New Modules

### 1. `SED_model`

This module allows you to compute the SED for a single star based on its fundamental properties. The input parameters include the star's effective temperature (`Teff`), metallicity (`FeH`), surface gravity (`logg`), radius, and distance from the observer. The output is a wavelength-flux relationship that can be adjusted for extinction using the Fitzpatrick extinction law.

#### Key Methods:
- **get_SED_spec**: Generates the flux of a star as a function of wavelength.
- **get_SED1 / get_SED2**: Get the SED for the first or second star (in the case of binary systems).
- **get_SED**: Returns the combined SED of two stars (binary system).
- **plot**: Plots the modelled SED against the observed data, if available, in both linear and logarithmic scales.

#### Example Usage:

```python
from stellarSpecModel import SED_model

# Initialize the model with parameters
sed_model = SED_model(teff1=6000, feh1=-0.5, logg1=4.0, R1=1.0, D=100)

# Get the SED for the first star
wave_SED1, SED1 = sed_model.get_SED1()

# Plot the SED
sed_model.plot(show=True)
```

### 2. `binary_SED_model`

This module handles the modeling of binary star systems by combining the individual SEDs of two stars. It computes the combined SED from both stars, considering parameters like their effective temperatures, metallicities, gravities, radii, and distances. Additionally, it can compute the chi-squared value and likelihood of the model given observed photometric data.

#### Key Methods:
- **get_SED_spec1 / get_SED_spec2**: Return the SED for each individual star in the binary system.
- **get_SED**: Computes the combined SED for both stars in the system.
- **add_data**: Allows the addition of photometric data (in magnitudes or fluxes) for comparing the model to observed data.
- **get_chisq**: Calculates the chi-squared statistic for the model fit to observed data.
- **get_lnlike**: Computes the log-likelihood for the observed data given the model.
- **plot**: Plots the combined SED model with any available observed data.

#### Example Usage:

```python
from stellarSpecModel import binary_SED_model

# Initialize the binary model with two stars
binary_model = binary_SED_model.BinarySEDModel(teff1=6000, feh1=-0.5, logg1=4.0, R1=1.0, 
                                               teff2=5000, feh2=-1.0, logg2=4.5, R2=1.0, D=100)

# Add photometric data (optional)
binary_model.add_data(bands=["V", "J"], obs_mags=[12.3, 13.1], obs_magerrs=[0.02, 0.02])

# Get the combined SED for the binary system
wave_SED, combined_SED = binary_model.get_SED()

# Plot the binary model SED
binary_model.plot(show=True)
```

## Requirements

To run `StellarSpecModel`, the following packages are required:

- `numpy`
- `pyphot`
- `matplotlib`
- `astropy`
- `extinction`
- [`spectool`](https://github.com/zhang-zhixiang/spectool)
- [`selenium`](https://github.com/SeleniumHQ/selenium) (optional, for downloading the grids)

## Installation

You can install `StellarSpecModel` directly from the GitHub repository:

```bash
pip install git+https://github.com/zhang-zhixiang/stellarSpecModel.git
```