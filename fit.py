"""Classes and other stuff for fitting."""

import sys

#FIXME: this is needed only for FitterMinuit
try: # if you have ROOT you might want minuit2
    from minuit2 import Minuit2 as Minuit
except:
    from minuit import Minuit

import utils 
from scipy.special import gammainc


class Fitter(object):
    """Superclass for fitting procedures/algorithms."""


class FitterMinuit(Fitter):
    """Fits using pyminuit."""



    def __init__(self, fitpoints, approach, model, **kwargs):
        self.fitpoints = fitpoints
        self.approach = approach

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

        # DMGLO
        self.m = Minuit(fcn, **model.parameter_dict)
        self.m.tol = 0.01
        #self.m.strategy = 2
        self.m.printMode = 0
        #self.m.maxcalls = 350


    def printres(self, printsigmas=0):
        """Print out the fitting result (chi-squares)."""
        nfreepars=utils.npars(self.m)
        dof = len(self.fitpoints) - nfreepars
        sigmas = [(getattr(self.approach, pt.yaxis)(pt, self.m.values) - pt.val) / pt.err for
                    pt in self.fitpoints]
        chi = sum(s*s for s in sigmas)  # equal to m.fval if minuit fit is done
        fitprob = (1.-gammainc(dof/2., chi/2.)) # probability of this chi-sq
        fitres = (chi, dof, fitprob)
        print 'P(chi-square, d.o.f) = P(%1.2f, %2d) = %5.4f' % fitres
        if printsigmas:
            print sigmas
            utils.prettyprint(sigmas)

    def fit(self):
        self.m.migrad()
        print "ncalls = ", self.m.ncalls


