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
        parameters_fixed: dict {'par0': bool, ...}. Is parameter value
            fixed? Considered False if non-existent. Once Fitter object
            is created, user should manipulate this dict only via Fitter methods.
        parameters_limits: dict {'par0: (float, float), ...}. Allowed range
            for fitting. Considered (-inf, inf) if non-existent. Once Fitter object
            is created, user should manipulate this dict only via Fitter methods.

    """
    parameters: dict = {}
    parameters_fixed: dict = {}
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


    def _release_parameters(self, *pars: str):
        """Release parameters for fitting.

        Args:
            *pars: Names of parameters to be released

        Notes:
            User should relase and fix parameters using methods of Fitter instance.
            This then calls this private ParameterModel method.

        """
        for par in pars:
            if par not in self.parameters:
                raise ValueError('Parameter {} is not defined in model {}'.format(
                        par, self))
            self.parameters_fixed[par] = False

    def _fix_parameters(self, *pars: str):
        """Fix parameters so they are not fitting variables.

        Args:
            *pars: Names of parameters to be fixed. If first name is 'ALL'
                then fix all parameters.

        Notes:
            User should relase and fix parameters using methods of Fitter instance.
            This then calls this private ParameterModel method.

        """
        if pars[0] == 'ALL':
            # fix 'em all
            for par in self.parameters:
                self.parameters_fixed[par] = True
        else:
            for par in pars:
                if par not in self.parameters:
                    raise ValueError('Parameter {} is not defined in model {}'.format(
                            par, self))
                self.parameters_fixed[par] = True

    def free_parameters(self) -> List[str]:
        """Return list of names of free fitting parameters."""
        return [p for p in self.parameters if p not in self.parameters_fixed
                or not self.parameters_fixed[p]]

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
