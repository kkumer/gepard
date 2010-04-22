"""
This is collection of tests. 
Just run 'nosetest' in the pype directory.

"""

import copy
from nose.tools import *
import numpy as np

import utils, Model, Approach, Fitter
from results import DMGLO1  #use some testpars here?

data = utils.loaddata('data/ep2epgamma')  

m = Model.ModelDR()
m.parameters.update(DMGLO1)
m.release_parameters('bS', 'Mv')
t = Approach.hotfixedBMK(m, optimization = False)

# testing data point
pt0 = copy.deepcopy(data[31][12])  # was data[1][0]
pt0.in1polarization = 1
pt0.in1charge = -1
pt0.to_conventions(Approach.BMK)
pt0.prepare(Approach.BMK)
# testing data set
testpoints = [data[31][12]] + [data[8][1]] + [data[29][2]] + [data[30][3]]
# testing data set for fits
fitpoints = data[31][12:14] + data[8][1:3] + data[30][2:4]

def test_CFF():
    """Calculate CFF H."""
    assert_almost_equal(m.ImH(pt0), 17.67971592396648)
    assert_almost_equal(m.ReH(pt0), -2.4699741916859592)

test_CFF.one = 1

def test_Xunp():
    """Calculate basic cross section Xunp."""
    assert_almost_equal(t.Xunp(pt0, vars={'phi':1.}), 1.8934179005138172)
    ar = t.Xunp(pt0, vars={'phi':np.array([0.3, 1.1])})
    assert isinstance(ar, np.ndarray)
    assert_almost_equal(ar[0], 1.4431498922585635)
    assert_almost_equal(ar[1], 1.9830803062066602)

def test_Xunp2():
    """Any kinematic variable can be in vars."""
    assert_almost_equal(t.Xunp(pt0, 
        vars={'phi':1., 'xB':0.07}), 3.0240991086297577)

def test_Xunp3():
    """ndarray of Q2 could be in vars."""
    ar = t.Xunp(pt0, vars={'phi':1, 'Q2':np.array([2.3, 2.5])})
    assert isinstance(ar, np.ndarray)
    assert_almost_equal(ar[0], 1.9818420283023106)
    assert_almost_equal(ar[1], 2.0973314901673334)

def test_Xunp4():
    """New feature: ndarray of any kinematical variable could be in vars."""
    ar = t.Xunp(pt0, vars={'phi':1, 'xB':np.array([0.07, 0.11])})
    assert isinstance(ar, np.ndarray)
    assert_almost_equal(ar[1], 0.68881435298587579)

test_Xunp4.newfeature = 1

def test_fit():
    """Testing set of fitting observables."""
    chisq = 0.
    for pt in testpoints:
        chisq = chisq + (
                (getattr(t, pt.yaxis)(pt) - pt.val)**2 / pt.err**2 )
    assert_almost_equal(chisq, 5.115294374919535)


def test_fit2():
    """Testing actual fitting by FitterMinuit."""
    f = Fitter.FitterMinuit(fitpoints, t)
    f.fit()
    chisq = t.chisq(fitpoints)[0]
    assert_almost_equal(chisq, 6.7638690027046771)

test_fit2.long = 1

def test_fit_neural():
    """Testing neural network fitting by FitterBrain."""
    np.random.seed(68)
    mNN = Model.ModelNN()
    tNN = Approach.hotfixedBMK(mNN)
    fNN = Fitter.FitterBrain(2*fitpoints, tNN, nnets=1, nbatch=6)
    fNN.fit()
    chisq = tNN.chisq(2*fitpoints)[0]
    assert_almost_equal(chisq, 20.894673139809264)

test_fit_neural.long = 1
