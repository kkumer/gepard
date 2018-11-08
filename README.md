## Synopsis

**Gepard** - package for working with generalized parton distributions (GPDs)

Modelling GPDs in momentum fraction or conformal moment space, perturbative QCD evolution and calculation of conformal Compton form factors (CFFs) up to NNLO accuracy. Modelling CFFs using dispersion relations. Calculation of DVCS and DVMP observables. Fitting to experimental data (both least-squares and neural nets).


## Installation

First make multiple copies of the pygepard extension library.
**Prerequisites**: C and Fortran compilers. 
```sh
make pygepards
```
Then go into pype directory and compile auxilliary library
```sh
cd pype
make optModel.so
```

## Using

Note that python code is still version 2!
Start (i)python v2 and 
run either `pype.py` which is generic or one of more specific example scripts in `ex` subdir.

For least-squares fitting you need Minuit python package. Tested with iminuit-1.3.3. (Hardwired limit on number of paramters maybe needs to be increased for some models. It's 50 I think.)

## Testing

If you are changing gepard code and want to be sure that you haven't broken something important, tests of many functions are available in subdir `test`. Fast suite of most important tests:
```
nosetests2 --rednose -vA "not newfeature and not long and not extendedtesting"
```


## License

GPL, once gepard becomes public.
