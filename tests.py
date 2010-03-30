"""
This is collection of tests. 
Just run 'nosetest' in the pype directory.

"""

import copy
from nose.tools import *
import numpy as np

import utils, models, Approach, fit
from results import DMGLO1  #use some testpars here?

ff = models.ModelDR()
ff.parameter_dict.update(DMGLO1)
ff.parameter_dict.update({'fix_bS':False, 'fix_Mv':False})
b = Approach.hotfixedBMK(ff, optimization = False)

data = utils.loaddata('data/ep2epgamma', approach=b)  
# testing data point
pt0 = copy.deepcopy(data[31][12])  # was data[1][0]
pt0.in1polarization = 1
pt0.in1charge = -1
pt0.to_conventions(b)
pt0.prepare(b)
# testing data set
testpoints = [data[31][12]] + [data[8][1]] + [data[29][2]] + [data[30][3]]
# testing data set for fits
fitpoints = data[31][12:14] + data[8][1:3] + data[30][2:4]

def test_CFF():
    assert_almost_equal(ff.ImH(pt0), 17.695607175490565)
    assert_almost_equal(ff.ReH(pt0), -2.529020432735325)
    # Old non-generic model
    #assert_almost_equal(ff.ImH(pt0), 17.681527454585797)
    #assert_almost_equal(ff.ReH(pt0), -2.471928504474433)

def test_Xunp():
    assert_almost_equal(b.Xunp(pt0, DMGLO1, vars={'phi':1.}), 1.8934179005138172)
    ar = b.Xunp(pt0, DMGLO1, vars={'phi':np.array([0.3, 1.1])})
    assert isinstance(ar, np.ndarray)
    assert_almost_equal(ar[0], 1.4431498922585635)
    assert_almost_equal(ar[1], 1.9830803062066602)

def test_Xunp2():
    """New feature: any kinematic variable could be in vars."""
    assert_almost_equal(b.Xunp(pt0, DMGLO1, 
        vars={'phi':1., 'xB':0.07}), 3.0240991086297577)

def test_Xunp3():
    """New feature: ndarray of Q2 could be in vars."""
    ar = b.Xunp(pt0, DMGLO1, vars={'phi':1, 'Q2':np.array([2.3, 2.5])})
    assert isinstance(ar, np.ndarray)
    assert_almost_equal(ar[0], 1.9818420283023106)
    assert_almost_equal(ar[1], 2.0973314901673334)

def test_Xunp4():
    """New feature: ndarray of any kinematical variable could be in vars."""
    ar = b.Xunp(pt0, DMGLO1, vars={'phi':1, 'xB':np.array([0.07, 0.11])})
    assert isinstance(ar, np.ndarray)
    assert_almost_equal(ar[1], 0.68881435298587579)

test_Xunp4.newfeature = 1

def test_fit():
    """Testing set of fitting observables."""
    chisq = 0.
    for pt in testpoints:
        chisq = chisq + (
                (getattr(b, pt.yaxis)(pt, DMGLO1) - pt.val)**2 / pt.err**2 )
    assert_almost_equal(chisq, 5.115294374919535)


def test_fit2():
    """Testing actual fitting by FitterMinuit."""
    f = fit.FitterMinuit(fitpoints, b, ff)
    f.fit()
    assert_almost_equal(f.m.fval, 6.7638634368267949, 4)

test_fit2.long = 1
