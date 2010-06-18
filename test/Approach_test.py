
import copy
from nose.tools import *
import numpy as np

import utils, Model, Approach
from results import DMGLO1  #use some testpars here?

data = utils.loaddata('data/ep2epgamma')  

m = Model.ModelDR()
m.parameters.update(DMGLO1)
t = Approach.hotfixedBMK(m, optimization = False)

# testing data point
pt0 = copy.deepcopy(data[31][12])  # was data[1][0]
pt0.in1polarization = 1
pt0.in1charge = -1
pt0.to_conventions(Approach.BMK)
pt0.prepare(Approach.BMK)


def test_CFF():
    """Calculate CFF H."""
    assert_almost_equal(m.ImH(pt0), 17.67971592396648)
    assert_almost_equal(m.ReH(pt0), -2.4699741916859592)

test_CFF.one = 1

def test_Xunp():
    """Calculate basic cross section Xunp."""
    assert_almost_equal(t.Xunp(pt0, vars={'phi':1.}),1.8872296859809268)
    ar = t.Xunp(pt0, vars={'phi':np.array([0.3, 1.1])})
    assert isinstance(ar, np.ndarray)
    assert_almost_equal(ar[0], 1.4412915895026261)
    assert_almost_equal(ar[1], 1.976296757099798)

def test_Xunp2():
    """Any kinematic variable can be in vars."""
    assert_almost_equal(t.Xunp(pt0, 
        vars={'phi':1., 'xB':0.07}), 3.0167858349753569)

def test_Xunp3():
    """ndarray of Q2 could be in vars."""
    ar = t.Xunp(pt0, vars={'phi':1, 'Q2':np.array([2.3, 2.5])})
    assert isinstance(ar, np.ndarray)
    assert_almost_equal(ar[0], 1.9769620102137611)
    assert_almost_equal(ar[1], 2.0929037204130285)

def test_Xunp4():
    """New feature: ndarray of any kinematical variable could be in vars."""
    ar = t.Xunp(pt0, vars={'phi':1, 'xB':np.array([0.07, 0.11])})
    assert isinstance(ar, np.ndarray)
    assert_almost_equal(ar[1], 0.68881435298587579)

test_Xunp4.newfeature = 1

