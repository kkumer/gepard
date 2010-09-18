
import copy
from nose.tools import *
import numpy as np

import utils, Model, Approach
from results import DMGLO1, DMepsGLO1  #use some testpars here?

data = utils.loaddata('data/ep2epgamma', approach=Approach.hotfixedBMK)  

m = Model.ModelDR()
m.parameters.update(DMGLO1)
t = Approach.hotfixedBMK(m, optimization = False)

mBM10 = Model.ModelDR()
mBM10.parameters.update(DMepsGLO1)
tBM10 = Approach.BM10(mBM10, optimization = False)

# testing data point for hotfixedBMK
pt0 = copy.deepcopy(data[31][12])  # was data[1][0]
pt0.in1polarization = 1
pt0.in1charge = -1

# testing data point for BM10
pt1 = copy.deepcopy(data[33][-1])
pt1.in2polarization = 1
pt1.phi = 2.*np.pi/5.
pt1.prepare(Approach.BM10)

def test_CFF():
    """Calculate CFF H."""
    assert_almost_equal(m.ImH(pt0), 17.67971592396648)
    assert_almost_equal(m.ReH(pt0), -2.4699741916859592)

test_CFF.one = 1

def test_Xunp():
    """Calculate basic cross section Xunp."""
    assert_almost_equal(t.Xunp(pt0, vars={'phi':1.}), 1.8872666478756337)
    ar = t.Xunp(pt0, vars={'phi':np.array([0.3, 1.1])})
    assert isinstance(ar, np.ndarray)
    assert_almost_equal(ar[0], 1.4413231946120821)
    assert_almost_equal(ar[1], 1.9763350136864286)

def test_Xunp2():
    """Any kinematic variable can be in vars."""
    assert_almost_equal(t.Xunp(pt0, 
        vars={'phi':1., 'xB':0.07}),  3.0168274215074025)

def test_Xunp3():
    """ndarray of Q2 could be in vars."""
    ar = t.Xunp(pt0, vars={'phi':1, 'Q2':np.array([2.3, 2.5])})
    assert isinstance(ar, np.ndarray)
    assert_almost_equal(ar[0], 1.9769930014185824)
    assert_almost_equal(ar[1], 2.0929323473733308)

def test_Xunp4():
    """New feature: ndarray of any kinematical variable could be in vars."""
    ar = t.Xunp(pt0, vars={'phi':1, 'xB':np.array([0.07, 0.11])})
    assert isinstance(ar, np.ndarray)
    assert_almost_equal(ar[1], 0.68881435298587579)

test_Xunp4.newfeature = 1

def test_XunpBM10():
    """Calculate unpolarized cross section Xunp in BM10 Approach."""
    assert_almost_equal(tBM10.PreFacSigma(pt1)*tBM10.TBH2unp(pt1), 0.01502937358336803)
    assert_almost_equal(tBM10.PreFacSigma(pt1)*tBM10.TDVCS2unp(pt1), 0.012565093106990456)
    assert_almost_equal(tBM10.PreFacSigma(pt1)*tBM10.TINTunp(pt1), 0.0011255158978939425)
    assert_almost_equal(tBM10.Xunp(pt1), 0.028719982588252427)

def test_XLP():
    """Calculate long. polarized cross section XLP in BM10 Approach."""
    assert_almost_equal(tBM10.PreFacSigma(pt1)*tBM10.TBH2LP(pt1), 0.009495908777414035)
    assert_almost_equal(tBM10.PreFacSigma(pt1)*tBM10.TDVCS2LP(pt1), -0.0032470111398419628)
    assert_almost_equal(tBM10.PreFacSigma(pt1)*tBM10.TINTLP(pt1), 0.0085102074298275109)
    assert_almost_equal(tBM10.XLP(pt1), 0.014759105067399584)
