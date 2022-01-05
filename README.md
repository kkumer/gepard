[![Build Status](https://app.travis-ci.com/kkumer/gepard.svg?token=LxGQXDrTERP6sTZEjfVh&branch=ps)](https://app.travis-ci.com/kkumer/gepard)
[![codecov](https://codecov.io/gh/kkumer/gepard/branch/ps/graph/badge.svg?token=YWFZALEK33)](https://codecov.io/gh/kkumer/gepard)

## Synopsis

**Gepard** - package for working with generalized parton distributions (GPDs)

(**This git branch is a complete rewrite of the code, and not yet ready for production use!**)

Modelling GPDs in momentum fraction or conformal moment space, perturbative QCD evolution and calculation of conformal Compton form factors (CFFs) up to NNLO accuracy. Modelling CFFs using dispersion relations. Calculation of DVCS and DVMP observables. Fitting to experimental data (both least-squares and neural nets).


## Installation

Just clone the git repository and then

```
cd gepard
pip install -e .
```


## Using

In Python, or, better, Jupyter,

```
import gepard as g
```

For least-squares fitting you need Minuit python package iminuit version 2.

See `tests` subdir for examples of use.

## Documentation

Available at  [https://calculon.phy.hr/gpd/docs/](https://calculon.phy.hr/gpd/docs/).
and in `docs` subdir.

## Developing

Before pushing anything to github master `pytest -q --runslow` *must* run without errors,
or commit message should indicate that code is not clean.

We aim for 100% documented, type-hinted code that passes `flake8`, 
`pydocstyle` (Google doc conventions) and `mypy` linters, but this
is not presently strictly required for commits.


## License

GPL, I guess, or LGPL, once gepard becomes public.
