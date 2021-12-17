"""Base classes for models of "soft" hadronic structure functions.

    * Elastic electromagnetic form factors
    * Generalized parton distribution functions (GPDs)
    * Compton form factors (for DVCS)
"""
from __future__ import annotations

from . import theory


class Model(theory.Theory):
    """Base class for all models.

    Instance of this class specifies structure of relevant hadrons.
    Methods provided are typically GPDs, CFFs, elastic FFs,
    DVMP transition TFFs, DAs, etc.
    Main subclasses are:

     - ParameterModel which depends on real parameters (some of which can be
       provided by minimization routine in the fitting procedure).
     - NeuralModel where structure functions are represented as neural nets
       (not yet implemented)
     - GridModel where values of structure functions are represented as
       grids of numbers, which may be interpolated
       (not yet implemented)

    These are then further subclassed to model actual structure functions.

    """
    def __init__(self, **kwargs) -> None:
        """Init Model object."""
        super().__init__(**kwargs)


class ParameterModel(Model):
    """Base class for all parametrized models.

    Attributes:
        parameters: dict {'par0': float, ...}. Actual value of parameter.
        parameters_fix: dict {'par0': bool, ...}. Is parameter value
            fixed? Considered False if non-existent.
        parameters_limits: dict {'par0: (float, float), ...}. Range for
            fitting. Considered (-inf, inf) if non-existent.
    """
    parameters: dict = {}
    parameters_fix: dict = {}
    parameters_limits: dict = {}

    def __init__(self, **kwargs) -> None:
        """Init ParameterModel object."""
        # print('ParameterModel init done')
        super().__init__(**kwargs)


    def add_parameters(self, newpars: dict):
        """Append newpars to parameters."""
        # Create parameters dict if it doesn't exist yet
        try:
            self.parameters
        except AttributeError:
            self.parameters = {}
        self.parameters.update(newpars)


    def release_parameters(self, *pars: str):
        """Release parameters for fitting.

        Args:
            *pars: Names of parameters to be released

        Notes:
            If allowed parameter ranges have to be changed from default, user needs
            to modify parameters dictionary directly.

        Todo:
            Note that relasing and fixing parameters for model instance is sensible
            only before creating Fitter instance! Afterwards, one has to fix and release
            params both for Fitter and for model.

        """
        for par in pars:
            if par not in self.parameters:
                raise ValueError('Parameter {} is not defined in model {}'.format(
                        par, self))
            self.parameters_fix[par] = False

    def fix_parameters(self, *pars: str):
        """Fix parameters so they are not fitting variables.

        Args:
            *pars: Names of parameters to be fixed. If first name is 'ALL'
                then fix all parameters.

        Todo:
            Note that relasing and fixing parameters for model instance is sensible
            only before creating Fitter instance! Afterwards, one has to fix and release
            params both for Fitter and for model.

        """
        if pars[0] == 'ALL':
            # fix 'em all
            for par in self.parameters:
                self.parameters_fix[par] = True
        else:
            for par in pars:
                if par not in self.parameters:
                    raise ValueError('Parameter {} is not defined in model {}'.format(
                            par, self))
                self.parameters_fix[par] = True

    def free_parameters(self) -> List[str]:
        """Return list of names of free fitting parameters."""
        return [p for p in self.parameters if p not in self.parameters_fix
                or not self.parameters_fix[p]]

    def print_parameters_errors(self, pvalues=False, ndof=0):
        """Print fitting parameters and their errors."""
        tolerance2 = 1   # sqrt(2*Npoints)/1.65^2 according to GJV/MRST
        for p in self.free_parameters():
            val = self.parameters[p]
            # err = sqrt(tolerance2)*sqrt(self.covariance[p,p])
            err = sqrt(tolerance2)*self.parameters_errors[p]
            if pvalues:
                pval = 2*(1.-scipy.stats.t.cdf(abs(val/err), ndof))
                print('%5s = %8.3f +- %5.3f  (p = %.3g)' % (p, val, err, pval))
            else:
                print('%5s = %8.3f +- %5.3f' % (p, val, err))
