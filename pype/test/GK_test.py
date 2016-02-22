"""Testing Goloskokov-Kroll model GK12 """

from nose.tools import *
import Model, Approach, Data

mGK = Model.GK12()
thGK = Approach.BM10(mGK)

gpdargs = (0.21, 0.13, -0.2, 7.)

pt = Data.DummyPoint()
pt.xi = 0.1/(2-0.1)
pt.t = -0.28
pt.Q2 = 2.

def test_GKGPDHuval():
    """Test GK model: H^u_val"""
    res = thGK.m.Huval(*gpdargs)
    assert_almost_equal(res, 2.249834779299078)

test_GKGPDHuval.extendedtesting = 1

def test_GKGPDHudsea():
    """Test GK model: H^u_sea"""
    res = thGK.m.Hudsea(*gpdargs)
    assert_almost_equal(res, 0.136954792369421)

test_GKGPDHudsea.extendedtesting = 1

def test_GKImH():
    """Test GK model: Im(CFFH)"""
    res = thGK.m.ImH(pt)
    assert_almost_equal(res, 13.1057763522748)

test_GKImH.extendedtesting = 1

def test_GKReH():
    """Test GK model: Re(CFFH)"""
    res = thGK.m.ReH(pt)
    assert_almost_equal(res, -1.94670099930845, 5)

#test_GKReH.extendedtesting = 1

def test_GKImHt():
    """Test GK model: Im(CFFHt)"""
    res = thGK.m.ImHt(pt)
    assert_almost_equal(res, 2.3675474323)

test_GKImHt.extendedtesting = 1

def test_GKReHt():
    """Test GK model: Re(CFFHt)"""
    res = thGK.m.ReHt(pt)
    assert_almost_equal(res, 1.0193853587666957)

test_GKReHt.extendedtesting = 1

def test_GKImE():
    """Test GK model: Im(CFFE)"""
    res = thGK.m.ImE(pt)
    assert_almost_equal(res, -3.3357207038)

test_GKImE.extendedtesting = 1

def test_GKReE():
    """Test GK model: Re(CFFE)"""
    res = thGK.m.ReE(pt)
    assert_almost_equal(res, -5.1576181607)

test_GKReE.extendedtesting = 1

def test_GKImEt():
    """Test GK model: Im(CFFEt)"""
    res = thGK.m.ImEt(pt)
    assert_almost_equal(res, 1.960464892)

test_GKImEt.extendedtesting = 1

def test_GKReEt():
    """Test GK model: Re(CFFEt)"""
    res = thGK.m.ReEt(pt)
    assert_almost_equal(res, 63.030901624)

test_GKReEt.extendedtesting = 1
