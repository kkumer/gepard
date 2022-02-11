"""Implementation of fitting algorithms."""

from __future__ import annotations

from typing import List
from iminuit import Minuit  # type: ignore

from . import data, theory


class Fitter(object):
    """Superclass for fitting procedures/algorithms."""

    def __init__(self, **kwargs) -> None:
        for key in kwargs:
            setattr(self, key, kwargs[key])


class MinuitFitter(Fitter):
    """Fits using iminuit Python frontend to MINUIT2 C++ library.

    Args:
        fitpoints (data.DataSet): points to fit to
        theory (theory.Theory): theory/model to be fitted

    """
    def __init__(self, fitpoints: data.DataSet,
                 theory: theory.Theory, **kwargs) -> None:
        self.fitpoints = fitpoints
        self.theory = theory
        init_vals = [v for v in self.theory.parameters.values()]
        names = [k for k in self.theory.parameters.keys()]

        def fcn(p):
            """Cost function for minimization - chi-square."""
            # Update model parameters.
            # Note: This relies on Python>3.7 where dicts are ordered!!
            self.theory.parameters.update(zip(self.theory.parameters.keys(), p))
            chisq = self.theory.chisq(self.fitpoints)
            return chisq

        fcn.errordef = Minuit.LEAST_SQUARES

        self.minuit = Minuit(fcn, init_vals, name=names)
        for k, v in self.theory.parameters_fixed.items():
            self.minuit.fixed[k] = v
        for k, v in self.theory.parameters_limits.items():
            self.minuit.limits[k] = v
        Fitter.__init__(self, **kwargs)

    def covsync(self):
        """Synchronize covariance and theory parameter errors with iminuit."""
        self.theory.parameters_errors = self.minuit.errors.to_dict()
        # iminuit gives covariance table for all parameters, here
        # we take only free ones:
        try:
            self.theory.covariance = {(p1, p2): self.minuit.covariance[p1, p2]
                                      for p1 in self.theory.free_parameters()
                                      for p2 in self.theory.free_parameters()}
        except AttributeError:
            print("Something's problematic. No covariances available.")

    def fit(self):
        """Perform simple fit.

        For better control, use iminuit's functions.

        """
        self.minuit.migrad()
        self.covsync()

    # The following methods keep status of parameters (fixed, limits)
    # in sync between Theory object and minuit object

    def fix_parameters(self, *args):  # noqa: D402
        """fix_parameters('p1', 'p2', ...)."""
        if args[0] == 'ALL':
            # fix 'em all
            self.theory._fix_parameters('ALL')
            for par in self.theory.model.parameters:
                self.minuit.fixed[par] = True
        else:
            self.theory._fix_parameters(*args)
            for par in args:
                self.minuit.fixed[par] = True

    def release_parameters(self, *args):  # noqa: D402
        """release_parameters('p1', 'p2', ...)."""
        self.theory._release_parameters(*args)
        for par in args:
            self.minuit.fixed[par] = False

    def limit_parameters(self, dct):  # noqa: D402
        """limit_parameters({'p1': (lo, hi), ...}."""
        self.theory.parameters_limits.update(dct)
        for k, v in dct.items():
            self.minuit.limits[k] = v

    def free_parameters(self) -> List[str]:
        """Return list of names of free fitting parameters."""
        self.theory.free_parameters()

    def print_parameters(self):
        """Values and errors for free parameters."""
        self.theory.print_parameters()
