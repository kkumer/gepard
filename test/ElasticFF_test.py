
from nose.tools import *

import Model

md = Model.ModelDR()  # dipole FFs
m = Model.ModelDRKelly()  # Kelly's FFs

mup = 2.7928473
mun = -1.913043

# Proton::

def test_DipoleFFp_zero():
    assert_almost_equal(md.F1(0.)/10., 1./10., 3)
    assert_almost_equal(md.F2(0.)/10., (mup-1.)/10., 2)

def test_DipoleFFp_one():
    assert_almost_equal(md.F1(-1.), 0.240568, 5)
    assert_almost_equal(md.F2(-1.), 0.241579, 5)

def test_KellyFFp_zero():
    assert_almost_equal(m.F1(0.), 1.)
    assert_almost_equal(m.F2(0.), mup-1., 6)

def test_KellyFFp_one():
    assert_almost_equal(m.F1(-1.), 0.238729, 5)
    assert_almost_equal(m.F2(-1.), 0.260398, 5)

# Neutron:

def test_KellyFFn_zero():
    assert_almost_equal(m.F1n(0.), 0.)
    assert_almost_equal(m.F2n(0.), mun, 6)

def test_KellyFFn_none():
    assert_almost_equal(m.F1n(-1.)*10., -0.0441597*10., 2)
    assert_almost_equal(m.F2n(-1.), -0.306797, 2)

