## Synopsis

**Gepard** - package for working with generalized parton distributions (GPDs)

This git branch contains the last version of old hybrid Python-Fortran code. 
You can use it (it's not easy to compile and run) if you need some features
which are not yet transfer to new pure Python version in master branch.
Important features which are only in this branch are:
   * NNLO evolution (partial) and NNLO hard scattering coefficients
   * Neural network fitting using PyBrain
   

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

Start (i)python and 
run either `pype.py` which is generic or one of more specific example scripts in `ex` subdir.

For least-squares fitting you need Minuit python package. Tested with iminuit-1.3.3. 

## Testing

If you are changing gepard code and want to be sure that you haven't broken something important, tests of many functions are available in subdir `test`. Fast suite of most important tests:
```
nosetests --rednose -vA "not newfeature and not long and not extendedtesting"
```


## License

GPL, once gepard becomes public.
