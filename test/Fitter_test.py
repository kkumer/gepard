
import copy
from nose.tools import *
import numpy as np

import utils, Model, Approach, Fitter
from results import DMGLO1  #use some testpars here?

data = utils.loaddata('data/ep2epgamma', approach=Approach.hotfixedBMK)  

m = Model.ModelDR()
m.parameters.update(DMGLO1)
m.parameters['limit_Mv'] = (0.9, 1.1)  # for compatibility with old g.
m.release_parameters('bS', 'Mv')
t = Approach.hotfixedBMK(m, optimization = False)

# Optimized model
mopt = Model.ModelDR(optimization = True)
mopt.parameters.update(DMGLO1)
mopt.ndparameters = np.array([mopt.parameters[name] for name in mopt.parameter_names])
mopt.parameters['limit_Mv'] = (0.9, 1.1)  # for compatibility with old g.
mopt.release_parameters('bS', 'Mv')
topt = Approach.hotfixedBMK(mopt, optimization = False)

# testing data set
testpoints = [data[31][12]] + [data[8][1]] + [data[29][2]] + [data[30][3]]
# testing data set for fits
fitpoints = data[31][12:14] + data[8][1:3] + data[30][2:4]


def test_fit():
    """Testing set of fitting observables."""
    chisq = 0.
    for pt in testpoints:
        chisq = chisq + (
                (getattr(t, pt.yaxis)(pt) - pt.val)**2 / pt.err**2 )
    assert_almost_equal(chisq, 8.8479922732922276)

def test_fit2():
    """Testing actual fitting by FitterMinuit."""
    f = Fitter.FitterMinuit(fitpoints, t)
    f.fit()
    chisq = t.chisq(fitpoints)[0]
    assert_almost_equal(chisq, 8.4891170857950087)

test_fit2.long = 1
test_fit2.batch = 1

def test_fit2opt():
    """Testing actual fitting by FitterMinuit."""
    fopt = Fitter.FitterMinuit(fitpoints, topt)
    fopt.fit()
    chisq = topt.chisq(fitpoints)[0]
    assert_almost_equal(chisq, 8.4891170857950087, 5)

test_fit2opt.long = 1
test_fit2opt.batch = 1
test_fit2opt.optimization = 1

def test_fit_neural():
    """Testing neural network fitting by FitterBrain."""
    np.random.seed(68)
    mNN = Model.ModelNN()
    tNN = Approach.hotfixedBMK(mNN)
    fNN = Fitter.FitterBrain(2*fitpoints, tNN, nnets=1, nbatch=6)
    fNN.fit()
    chisq = tNN.chisq(2*fitpoints)[0]
    assert_almost_equal(chisq, 22.754514666801569)

test_fit_neural.long = 1
