# pIMF
pIMF (python Initial Mass Function, pronounced pie-em-eff) is a small library designed for interacting with intial mass functions of stars.

The current scope is generating and integrating over different IMFs.

This README is split into different sections:
1. [Definitions and Choices](#definitions-and-choices)
2. [Available functional forms](#available-functional-forms)
3. [Installation](#installation)
4. [Brief tutorial](#tutorial)
5. [Contributing](#contributing)
6. [(no) Testing](#testing)

## Definitions and Choices
There are several ways one can define an IMF. I have tried to be consistent as defining my IMFs as the number of stars that form in a given mass bin: $$\mathrm{d}N = \xi(M)\mathrm{d}M$$.

This means that the `integrate` method calculates the total **number** of stars that form in a given mass range, and `integrate_product` calculates the total **mass**.

Note the definition also uses $\mathrm{d}M$, not $\mathrm{d}\log M$.

These differences effectively just introduce/remove factors of $M$ so it's not hard to switch between them, but it is important to make sure you get it right.

## Available functional forms
We currently have implemented several functional forms:
- Power law, like that of [Salpeter (1955)](https://ui.adsabs.harvard.edu/abs/1955ApJ...121..161S/abstract)
- Broken power law, similar to [Kroupa (2001)](https://ui.adsabs.harvard.edu/abs/2001MNRAS.322..231K/abstract).
- Power law with an exponential cutoff, suggested by [Wise+ (2012)](https://ui.adsabs.harvard.edu/abs/2012ApJ...745...50W/abstract)
- Lognormal distribution, as in [Larson (1973)](https://ui.adsabs.harvard.edu/abs/1973MNRAS.161..133L/abstract)
- [Chabrier (2003)](https://ui.adsabs.harvard.edu/abs/2003PASP..115..763C/abstract), which is defined piecewise as a lognormal distribution below a critical mass, and as a power law above.
- The [Maschberger (2013)](https://ui.adsabs.harvard.edu/abs/2013MNRAS.429.1725M/abstract) L3 form

Check out [Contributing](#contributing) if you want anything not listed above!

## Installation
I will hopefully publish to [pypi](https://pypi.org) and [conda](https://docs.conda.io/en/latest/) once I have a setup.py file in version 1.0.

### From source
Download the source code, e.g. by cloning this repo.

Run `pip install .`

For everything needed to run tests, build documentation, etc. install the development version (note the quotes):

`pip install '.[dev]'`

### Requirements
For normal use, only `numpy` and `scipy` are required. The `pyproject.toml` file includes a full list.

## Tutorial
It is hopefully pretty easy to get started with pIMF. The below code shows you how to plot a few different literature IMFs to compare them. See also our [Read the Docs](https://pimf.readthedocs.io/en/latest/)

```python
from pimf import PowerLawIMF, ChabrierIMF, BrokenPowerLawIMF

import numpy as np
import matplotlib.pyplot as plt

s95 = PowerLawIMF(
    normalisation="mass",  # What should the IMF be normalised to?
    normalisation_value=1  # I will define this to be one solar mass
)
k01 = BrokenPowerLawIMF(normalisation="mass")  # Normalisation value is 1 by default
c03 = ChabrierIMF(normalisation="mass")

M = np.geomspace(0.1, 100)
plt.loglog(
    M,
    s95(M),  # The __call__ method will evaluate our IMF at that mass.
    label="Salpeter (1955)"
    )
plt.loglog(M, k01(M), label="Kroupa (2001)")
plt.loglog(M, c03(M), label="Chabrier (2003)")
plt.xlabel("Mass")
plt.ylabel("Number of stars")
plt.legend()
plt.show()
```

## Contributing
If there are any features/[functional forms](#available-functional-forms) currently missing that you would like to see, feel free to [create an issue](https://github.com/E-W-Jones/pIMF/issues/new), or even [create your own pull request](https://github.com/E-W-Jones/pIMF/pulls)!

## Testing
There is currently no testing infrastructure, however I do have some random functions/plots to check things look right. [This is a known issue](https://github.com/E-W-Jones/pIMF/issues/2).