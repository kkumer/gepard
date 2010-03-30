"""Classes and other stuff for fitting."""

import sys

#FIXME: this is needed only for FitterMinuit
try: # if you have ROOT you might want minuit2
    from minuit2 import Minuit2 as Minuit
except:
    from minuit import Minuit


class Fitter(object):
    """Superclass for fitting procedures/algorithms."""


class FitterMinuit(Fitter):
    """Fits using pyminuit."""



    def __init__(self, fitpoints, approach, model, **kwargs):
        self.fitpoints = fitpoints
        self.approach = approach
        self.model = model

        # FIXME: ugly hack because Minuit counts the arguments of fcn so 'self'
        #        is not allowed
        # fcnargs = "NS, alS, ..."
        # pardict = "'NS': NS, 'alS': alS, ..."
        fcnargs = ", ".join(model.parameter_names) 
        pardict = ", ".join(map(lambda x: "'%s': %s" % x, 
                          zip(model.parameter_names, model.parameter_names)))
        exec(
"""
def fcn(%s):
    pars = {%s}
    chisq = 0.
    for pt in fitpoints:
        chisq = chisq + (
                (getattr(approach, pt.yaxis)(pt, pars) - pt.val)**2 / pt.err**2 )
    return chisq
""" % (fcnargs, pardict), locals(),locals())
        self.m = Minuit(fcn, **model.parameter_dict)
        for key in kwargs:
            setattr(self.m, key, kwargs[key])


    def fit(self):
        self.m.migrad()
        print "ncalls = ", self.m.ncalls
        print self.model.prettyprint(self.fitpoints, self.approach)
        return self.model


