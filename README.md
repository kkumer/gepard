[![Tests](https://github.com/kkumer/gepard/actions/workflows/python-package.yml/badge.svg)](https://github.com/kkumer/gepard/actions/workflows/python-package.yml)
[![codecov](https://codecov.io/gh/kkumer/gepard/branch/ps/graph/badge.svg)](https://codecov.io/gh/kkumer/gepard)

![Gepard logo](docs/source/media/Gepard_logo.png)

## Synopsis

**Gepard** - tool for studying the 3D quark and gluon distributions in the nucleon


   * Modelling Generalized Parton Distributions (GPD) and Compton form factors (CFF).
   * Perturbative NLO QCD evolution of GPDs
   * Calculation of deeply virtual Compton scattering (DVCS) and deeply virtual meson production (DVMP) observables to NLO accuracy.
   * Fitting parametrized models to the experimental data.


## Installation

Gepard is a Python package available from standard PyPI repository so you
should be able to install the latest version together with the required
dependencies by simply issuing

```
pip install gepard
```

Of course, you should have Python first (version 3.7 at the least), 
as well as `pip`, which all major Linux distributions have as standard package.
See the documentation for alternative instalation from sources, which is also easy.

Gepard is developed on Linux, but is tested also on Windows and MacOS. 
Please complain if you have problems.


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
