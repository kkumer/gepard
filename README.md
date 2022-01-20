[![Build Status](https://app.travis-ci.com/kkumer/gepard.svg?token=LxGQXDrTERP6sTZEjfVh&branch=ps)](https://app.travis-ci.com/kkumer/gepard)
[![codecov](https://codecov.io/gh/kkumer/gepard/branch/ps/graph/badge.svg?token=YWFZALEK33)](https://codecov.io/gh/kkumer/gepard)

## Synopsis

**Gepard** - tool for studying the 3D quark and gluon distributions in the nucleon


   * Modelling Generalized Parton Distributions (GPD) and Compton form factors (CFF).
   * Perturbative NLO QCD evolution of GPDs
   * Calculation of deeply virtual Compton scattering (DVCS) and deeply virtual meson production (DVMP) observables to NLO accuracy.
   * Fitting parametrized models to the experimental data.


## Installation

You need Python version 3.7 at the least.
Just clone the git repository and then

```
cd gepard
pip install -e .
```

For Python packages which are required to run Gepard code, see
the file [requirements.txt](requirements.txt).

For the public release (real soon now), we intend to create
the public PyPI package for easy installation.

## Using

See the documentation,
available at  [https://gepard.phy.hr/](https://gepard.phy.hr/).
Also there are lots of examples in the [tests](tests/) subdirectory.

## Citing

If you use Gepard for your scientific publications, you might want to cite the
[arXiv:hep-ph/0703179](https://arxiv.org/abs/hep-ph/0703179) paper where theoretical
framework that is implemented in the Gepard code is described in some detail.


## Copyright and license

Copyright (C) 2022 Krešimir Kumerički

This program is free software: you can redistribute it and/or modify it under the terms of the 
GNU Affero General Public License v3.0 or later. See the file [LICENCE](LICENSE) for details.
