"""Tests for electromagnetic elastic form factors."""

import gepard as g
from pytest import approx, raises

md = g.eff.DipoleEFF()  # dipole FFs
m = g.eff.KellyEFF()  # Kelly's FFs

# proton
pt = g.DataPoint(in2particle='p', Q2=4., t=0., xB=0.1)
pt1 = g.DataPoint(in2particle='p', Q2=4., t=-1., xB=0.1)

# neutron
npt = g.DataPoint(in2particle='n', Q2=4., t=0., xB=0.1)
npt1 = g.DataPoint(in2particle='n', Q2=4., t=-1., xB=0.1)

mup = 2.7928473
mun = -1.913043

# Proton::

def test_DipoleFFp_zero():
    assert md.F1(pt)/10. == approx(1./10., 3)
    assert md.F2(pt)/10. == approx((mup-1.)/10., rel=1.e-2)

def test_DipoleFFp_one():
    assert md.F1(pt1) == approx(0.240568, rel=1.e-5)
    assert md.F2(pt1) == approx(0.241579, rel=1.e-5)

def test_KellyFFp_zero():
    assert m.F1(pt) == 1.
    assert m.F2(pt) == approx(mup-1., rel=1.e-6)

def test_KellyFFp_one():
    assert m.F1(pt1) == approx(0.238729, rel=1.e-5)
    assert m.F2(pt1) == approx(0.260398, rel=1.e-5)

# Neutron:
def test_DipoleFFn_exception():
    with raises(Exception,
                match='Neutron dipole elastic FFs are not implemented yet! Use Kelly.'):
        md.F1(npt) == 0.
    with raises(Exception,
                match='Neutron dipole elastic FFs are not implemented yet! Use Kelly.'):
        md.F2(npt) == 0.

def test_KellyFFn_zero():
    assert m.F1(npt) == 0.
    assert m.F2(npt) == approx(mun, rel=1.e-6)

def test_KellyFFn_none():
    assert m.F1(npt1)*10. == approx(-0.0441597*10., rel=1.e-2)
    assert m.F2(npt1) == approx(-0.306797, rel=1.e-2)
