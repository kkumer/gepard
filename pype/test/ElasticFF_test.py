
from nose.tools import *

import Model, Data

md = Model.ModelDR()  # dipole FFs
m = Model.ModelDRKelly()  # Kelly's FFs

# proton
pt = Data.DummyPoint(init = {'in2particle':'p', 'Q2':4., 't':0., 'xB':0.1})
pt1 = Data.DummyPoint(init = {'in2particle':'p', 'Q2':4., 't':-1., 'xB':0.1})

# neutron
npt = Data.DummyPoint(init = {'in2particle':'n', 'Q2':4., 't':0., 'xB':0.1})
npt1 = Data.DummyPoint(init = {'in2particle':'n', 'Q2':4., 't':-1., 'xB':0.1})

mup = 2.7928473
mun = -1.913043

# Proton::

def test_DipoleFFp_zero():
    assert_almost_equal(md.F1(pt)/10., 1./10., 3)
    assert_almost_equal(md.F2(pt)/10., (mup-1.)/10., 2)

def test_DipoleFFp_one():
    assert_almost_equal(md.F1(pt1), 0.240568, 5)
    assert_almost_equal(md.F2(pt1), 0.241579, 5)

def test_KellyFFp_zero():
    assert_almost_equal(m.F1(pt), 1.)
    assert_almost_equal(m.F2(pt), mup-1., 6)

def test_KellyFFp_one():
    assert_almost_equal(m.F1(pt1), 0.238729, 5)
    assert_almost_equal(m.F2(pt1), 0.260398, 5)

# Neutron:

def test_KellyFFn_zero():
    assert_almost_equal(m.F1(npt), 0.)
    assert_almost_equal(m.F2(npt), mun, 6)

def test_KellyFFn_none():
    assert_almost_equal(m.F1(npt1)*10., -0.0441597*10., 2)
    assert_almost_equal(m.F2(npt1), -0.306797, 2)

