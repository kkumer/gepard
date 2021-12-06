"""Implementation of fitting algorithms."""


from iminuit import Minuit

import gepard.data
import gepard.theory


class Fitter(object):
    """Superclass for fitting procedures/algorithms."""

    def __init__(self, **kwargs) -> None:
        for key in kwargs:
            setattr(self, key, kwargs[key])


class FitterMinuit(Fitter):
    """Fits using iminuit Python frontend to MINUIT2 C++ library."""

    def __init__(self, fitpoints: gepard.data.DataSet,
                 theory: gepard.theory.Theory, **kwargs) -> None:
        """Set what is fitted to what and how."""
        self.fitpoints = fitpoints
        self.theory = theory
        init_vals = [v for v in self.theory.m.parameters.values()]
        names = [k for k in self.theory.m.parameters.keys()]

        def fcn(p):
            """Cost function for minimization - chi-square."""
            # Update model parameters.
            # FIXME: This relies on Python>3.7 where dicts are ordered!!
            #  It's OK, one just needs to specify Python version dependency
            self.theory.m.parameters.update(zip(self.theory.m.parameters.keys(), p))
            chisq = self.theory.chisq(self.fitpoints)
            return chisq

        fcn.errordef = Minuit.LEAST_SQUARES

        self.minuit = Minuit(fcn, init_vals, name=names)
        for k, v in self.theory.m.parameters_limits.items():
            self.minuit.limits[k] = v
        Fitter.__init__(self, **kwargs)

    def fit(self):
        """Start fitting."""
        self.minuit.migrad()
        self.theory.model.parameters_errors = self.minuit.errors.to_dict()
        self.theory.model.covariance = self.minuit.covariance.to_table()

    def fix_parameters(self, *args):
        """fix_parameters('p1', 'p2', ...)."""

        if args[0] == 'ALL':
            # fix 'em all
            self.theory.model.fix_parameters('ALL')
            for par in self.theory.model.parameters:
                self.minuit.fixed[par] = True
        else:
            self.theory.model.fix_parameters(*args)
            for par in args:
                self.minuit.fixed[par] = True

    def release_parameters(self, *args):
        """release_parameters('p1', 'p2', ...)."""

        self.theory.model.release_parameters(*args)
        for par in args:
            self.minuit.fixed[par] = False

    def limit_parameters(self, dct):
        """limit_parameters({'p1': (lo, hi), ...}."""

        self.theory.model.parameters_limits.update(dct)
        for k, v in dct.items():
            self.minuit.limits[k] = v

    def print_parameters(self):
        for par in self.theory.model.parameters:
            if not self.theory.model.parameters['fix_'+par]:
                print('%5s = %4.2f +- %4.2f' % (par,
                        self.minuit.values[par], self.minuit.errors[par]))
