"""Base classes for models of "soft" hadronic structure functions.

    * Elastic electromagnetic form factors
    * Generalized parton distribution functions (GPDs)
    * Compton form factors (for DVCS)

Actual implementations are in other files. Here only generic
treatment of generalized parametrized model is coded.

Todo:
    * Model defined by grid of numbers
    * Flavored models

"""
from __future__ import annotations

from math import sqrt
from typing import List

from . import theory
from .constants import tolerance2


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
    def __init__(self, **kwargs) -> None:
        # If subclases don't create these dicts, this is the
        # last chance. They have to exist.
        for d in ['parameters', 'parameters_fixed', 'parameters_limits']:
            if not hasattr(self, d):
                setattr(self, d, {})
        # By default, all params are fixed
        self._fix_parameters('ALL')
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

    def add_parameters_limits(self, newlimits: dict):
        """Append newlimits to parameters_limits."""
        # Create parameters_limits dict if it doesn't exist yet
        try:
            self.parameters_limits
        except AttributeError:
            self.parameters_limits = {}
        self.parameters_limits.update(newlimits)

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

    def print_parameters(self):
        """Print fitting parameters and their errors."""
        for p in self.free_parameters():
            val = self.parameters[p]
            err = sqrt(tolerance2)*self.parameters_errors[p]
            print('{:5s} = {:8.3f} +- {:5.3f}'.format(p, val, err))
