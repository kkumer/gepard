## Synopsis

**Gepard** - package for working with generalized parton distributions (GPDs)

(**This git branch is a complete rewrite of the code, and not yet ready for production use!**)

Modelling GPDs in momentum fraction or conformal moment space, perturbative QCD evolution and calculation of conformal Compton form factors (CFFs) up to NNLO accuracy. Modelling CFFs using dispersion relations. Calculation of DVCS and DVMP observables. Fitting to experimental data (both least-squares and neural nets).


## Installation

Just clone the git repository. Note that presently you need Fortran extension `pygepard.so`
in the package directory. 
You can obtain it from the author for the Linux `x86_64` architecture or compile it 
yourself using `master` branch.

There is no Python `setup.py` yet, so you need to add the package directory
to your Python path, see below.

## Using

```
import sys
sys.path.append('<path to package dir>')
import gepard as g
```

For least-squares fitting you need Minuit python package. Tested with iminuit-1 (version 2 
will *not* work).

See `tests` subdir for examples of use.

## Documentation

Available at  [https://calculon.phy.hr/gpd/docs/](https://calculon.phy.hr/gpd/docs/).
and in `docs` subdir.

## Developing

Before pushing anything to github master `pytest` *must* run without errors,
or commit message should indicate that code is not clean.
We aim for 100% documented, type-hinted code that passes `flake8`, 
`pydocstyle` (Google doc conventions) and `mypy` linters, but this
is not presently strictly required for commits.


## License

GPL, once gepard becomes public.
