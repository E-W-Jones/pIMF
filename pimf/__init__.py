"""
pIMF
====

pIMF (python Initial Mass Functions, pronounced pie-em-eff) is a small library for manipulating stellar Initial Mass Functions.

Quick Start
-----------
To find the percentage of stars between 1 and 100 solar masses:

```python
>>> from pimf import PowerLawIMF
>>> s55 = PowerLawIMF()
>>> s55.integrate(1, 100) / s55.integrate(0.1, 100)
>>> 0.04458320760384313
```

All of the IMFs in {py:class}`pimf.initialmassfunction` are also exposed here:
```{autodoc2-summary}
   :renderer: myst
   ~pimf.InitialMassFunction
   ~pimf.PowerLawIMF
   ~pimf.ChabrierIMF
   ~pimf.BrokenPowerLawIMF
   ~pimf.L3IMF
   ~pimf.LognormalIMF
   ~pimf.ExponentialCutoffPowerLawIMF
```
"""

from . import initialmassfunction
from .initialmassfunction import (
    InitialMassFunction,
    PowerLawIMF,
    ChabrierIMF,
    BrokenPowerLawIMF,
    L3IMF,
    LognormalIMF,
    ExponentialCutoffPowerLawIMF
    )

__version__ = "0.1.0"

__all__ = initialmassfunction.__all__
