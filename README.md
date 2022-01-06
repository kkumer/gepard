[![Build Status](https://app.travis-ci.com/kkumer/gepard.svg?token=LxGQXDrTERP6sTZEjfVh&branch=ps)](https://app.travis-ci.com/kkumer/gepard)
[![codecov](https://codecov.io/gh/kkumer/gepard/branch/ps/graph/badge.svg?token=YWFZALEK33)](https://codecov.io/gh/kkumer/gepard)

## Synopsis

**Gepard** - tool for studying the 3D quark and gluon distributions in the nucleon


Modelling Generalized Parton Distributions (GPD) and Compton form factors (CFF).
Perturbative QCD evolution and calculation of deeply virtual
Compton scattering (DVCS) and deeply virtual meson production (DVMP) to NLO accuracy. 
Fitting parametrized models to experimental data.


## Installation

Just clone the git repository and then

```
cd gepard
pip install -e .
```

For Python packages which are required to run Gepard code, see
the file [requirements.txt](requirements.txt).

For public release, we intend to create public PyPI and conda packages.

## Using

In Python, or, better, Jupyter, import Gepard like this

```
import gepard as g
```

See documentation, 
available at  [https://calculon.phy.hr/gpd/docs/](https://calculon.phy.hr/gpd/docs/).
and in `docs` subdir.
Also there are lot of examples in the [tests](tests/) subdirectory.

## Citing

If you use Gepard for your scientific publications, you might cite the
[paper](https://arxiv.org/abs/hep-ph/0703179) where most of the formalism
that is implemented in the Gepard code is described.


## License

GNU [AGPLv3](https://www.gnu.org/licenses/why-affero-gpl.en.html).
