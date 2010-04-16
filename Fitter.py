"""Classes for fitting."""

import sys

import numpy as np

#FIXME: this is needed only for FitterMinuit. It should
# not raise exception on systems without pyminuit installed
# after another Fitter is implemented, say NN
try: # if you have ROOT you might want minuit2
    from minuit2 import Minuit2 as Minuit
except:
    from minuit import Minuit

from pybrain.tools.shortcuts import buildNetwork
import brain
import trans  # output layer transformation for FitterBrain


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
        self.theory.print_chisq(self.fitpoints)
        self.theory.model.print_parameters()
        return self.theory


class FitterBrain(Fitter):
    """Fits using PyBrain neural net library."""

    def __init__(self, fitpoints, theory, **kwargs):
        self.fitpoints = fitpoints
        self.theory = theory
        self.verbose = 0

    def artificialData(self, datapoints, trainsize=13):
        """Create artificial data replica.
        
        Replica is created by randomly picking value around mean value taken from
        original data, using normal Gaussian distribution with width equal to
        uncertainty of the original data point. Resulting set of artificial
        datapoints is then shuffled and divided into two SupervisedDataSet
        instances: training and testing, which are returned.
        Input datapoints can be DataSet instance or just a list of DataPoints
        instances.

        Keyword arguments:
        trainsize -- size of subset used for training (rest is for testing)

           
        """
        training = brain.SupervisedDataSetTransformed(2, 1)  # FIXME: get size from data
        testing = brain.SupervisedDataSetTransformed(2, 1)
        i = 0
        trans.map2pt.clear()
        for pt in np.random.permutation(datapoints):
            xs = [pt.xB, pt.t, pt.Q2]
            # FIXME: This abs() below is for HERMES->BKM. Should be done using info
            # from .dat Rounding the number, to make matching of trans.map2pt work
            # regardless of computer rounding behaviour
            y = [pt.val + round(np.random.normal(0, pt.err, 1)[0], 5)]
            # Numerical derivative of transformation function w.r.t. net output:
            deriv = (self.theory.predict(pt, parameters={'outputvalue':1.2}) -
                     self.theory.predict(pt, parameters={'outputvalue':1.0})) / 0.2
            trans.map2pt[y[0]] = (self.theory, pt, deriv)
            # FIXME: trainsize should be specified by percentage and not by value
            if i < trainsize:
                training.addSample(xs[:-1], y) # we don't use Q2 for training
            else:
                testing.addSample(xs[:-1], y)
            i += 1
        return training, testing

    def makenet(self, datapoints):
        """Create trained net and return tuple (net, error)."""

        dstrain, dstest = self.artificialData(datapoints)

        net = buildNetwork(2, 7, 1)

        t = brain.RPropMinusTrainerTransformed(net, learningrate = 0.9, lrdecay = 0.98, 
                momentum = 0.0, batchlearning = True, verbose = False)

        #t = brain.BackpropTrainerTransformed(net, learningrate = 0.7, lrdecay = 0.98,
        #        momentum = 0.1 , batchlearning = False, verbose = False)

        # Train in batches of batchlen epochs and repeat nbatch times
        nbatch = 20
        batchlen = 5
        memerr = 1.  # large initial error, certain to be bettered
        for k in range(nbatch):
            t.trainOnDataset(dstrain, batchlen)
            trainerr, testerr = (t.testOnData(dstrain), t.testOnData(dstest))
            if testerr < memerr:
                memerr = testerr
                memnet = net
                if self.verbose:
                    print "Epoch: %6i   ---->    Error: %8.3g  TestError: %8.3g" % (
                            t.epoch, trainerr, testerr)
        return net, memerr
    
    def fit(self, nnets=12):
        """Create and train nnets (default: 12) neural networks."""
        for n in range(nnets):
            net, memerr = self.makenet(self.fitpoints)
            print "Net No. %2i  --->  TestError: %8.3g" % (n, memerr)
            self.theory.model.nets.append(net)
        return self.theory

