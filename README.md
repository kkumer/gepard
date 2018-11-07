## Synopsis

**Gepard** - package for working with generalized parton distributions (GPDs)

Modelling GPDs in momentum fraction or conformal moment space, perturbative QCD evolution and calculation of conformal Compton form factors (CFFs) up to NNLO accuracy. Modelling CFFs using dispersion relations. Calculation of DVCS and DVMP observables. Fitting to experimental data (both least-squares and neural nets).


## Installation

First make the pygepard extension library.
**Prerequisites**: C and Fortran compilers. 
```sh
make pygepards
```
Then go into pype directory and compile auxilliary library
```sh
make optModel.so
```

## Using

Note that python code is still version 2!
Start (i)python v2 and 
run either `pype.py` which is generic or one of more specific example scripts in `ex` subdir.


## License

GPL, once gepard becomes public.
