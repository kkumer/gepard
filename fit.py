"""Classes and other stuff for fitting."""

import sys

#FIXME: this is needed only for FitterMinuit. It should
# not raise exception on systems without pyminuit installed
# after another Fitter is implemented, say NN
try: # if you have ROOT you might want minuit2
    from minuit2 import Minuit2 as Minuit
except:
    from minuit import Minuit


class Fitter(object):
    """Superclass for fitting procedures/algorithms."""


class FitterMinuit(Fitter):
    """Fits using pyminuit."""



    def __init__(self, fitpoints, theory, **kwargs):
        self.fitpoints = fitpoints
        self.theory = theory

        # FIXME: ugly hack because Minuit counts the arguments of fcn so 'self'
        #        is not allowed
        # fcnargs = "NS, alS, ..."
        # pardict = "'NS': NS, 'alS': alS, ..."
        fcnargs = ", ".join(theory.model.parameter_names) 
        pardict = ", ".join(map(lambda x: "'%s': %s" % x, 
                          zip(theory.model.parameter_names, theory.model.parameter_names)))
        exec(
"""
def fcn(%s):
    theory.model.parameters.update({%s})
    chisq = 0.
    for pt in fitpoints:
        chisq = chisq + (
                (getattr(theory, pt.yaxis)(pt) - pt.val)**2 / pt.err**2 )
    return chisq
""" % (fcnargs, pardict), locals(),locals())
        self.minuit = Minuit(fcn, **theory.model.parameters)
        for key in kwargs:
            setattr(self.minuit, key, kwargs[key])


    def fit(self):
        self.minuit.migrad()
        print "ncalls = ", self.minuit.ncalls
        self.theory.model.print_chisq(self.fitpoints, self.theory)
        self.theory.model.print_parameters()
        return self.theory.model


