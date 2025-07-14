# pIMF Home

## Contents
```{toctree}
:maxdepth: 1
self
tutorial
subclassing
theory
apidocs/index
```

## Available IMFs
There are several IMFs available right out of the box:
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

<!-- 
.. pIMF documentation master file, created by
   sphinx-quickstart on Sat Jul  5 15:16:33 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pIMF documentation
==================

Add your content using ``reStructuredText`` syntax. See the
`reStructuredText <https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html>`_
documentation for details.

Contents
--------
.. toctree::
   self
   api

.. automodule:: pimf -->