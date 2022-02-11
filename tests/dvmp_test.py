"""Tests for DVEM."""

import gepard as g
from gepard.fits import par_KM10b
from pytest import approx, fixture, mark

par_dvmp = {'ns':  0.152, 'al0s': 1.158, 'alps': 0.15, 'ms2': 0.446,
            'secs': -0.442, 'this': 0.089,  # 'ng': 0.5,  # provided by ns
            'al0g': 1.247, 'alpg': 0.15, 'mg2': 0.7, 'secg': -2.309, 'thig': 0.812}


class MyTheory(g.gpd.PWNormGPD, g.dvmp.MellinBarnesTFF, g.dvmp.DVMP):
    pass

@fixture
def th_dvmp():
    """Generic LO model."""
    th = MyTheory()
    th.parameters.update(par_dvmp)
    return th


@fixture
def th_nlo():
    """Generic NLO model."""
    th = MyTheory(p=1)
    th.parameters.update(par_dvmp)
    return th


@fixture
def th_KM10b():
    """KM10b model."""
    th = MyTheory()
    th.parameters.update(par_KM10b)
    return th


def test_dvmp_TFFs_LO(th_dvmp):
    """Calculate LO DVMP TFFs for rho0 production at input scale."""
    xB = 1e-4
    pt = g.data.DataPoint({'Q2': 4., 't': 0, 'xB': xB})
    tffs = th_dvmp.tff(pt.xi, pt.t, pt.Q2)
    reh, imh = tffs[0], tffs[1]
    # following agrees with DM to best than percent
    assert imh == approx(12395.53, rel=1e-7)
    assert reh == approx(4766.8993, rel=1e-7)


def test_gepardTFFsEvol(th_KM10b):
    """Calculate LO DVMP TFFs for rho0 production + evolution."""
    pt = g.data.DataPoint({'Q2': 6.6, 'W': 75., 't': -0.025})
    tffs = th_KM10b.tff(pt.xi, pt.t, pt.Q2)
    reh, imh = tffs[0], tffs[1]
    assert reh == approx(285.15512, rel=1e-2)
    assert imh == approx(511.39622404, rel=1e-2)


def test_gepardXrho0t(th_KM10b):
    """Calculate LO DVMP cross section d sigma / dt"""
    pt = g.data.DataPoint({'Q2': 6.6, 'W': 75., 't': -0.025,
                           'process': 'gammastarp2rho0p'})
    assert th_KM10b.XGAMMA(pt) == approx(1212.62165, rel=1.e-2)


def test_c1_NLO(th_dvmp):
    """Calculate NLO DVMP hard scattering coefficients"""
    xB = 1e-4
    pt = g.data.DataPoint({'Q2': 4., 't': 0, 'xB': xB})
    # to get agreement with these old numbers:
    th_dvmp.corr_c1dvmp_one = 0
    th_dvmp.corr_c1dvmp_sgn = -1
    # c1dvmp(model, sgntr, j, k)
    aux = g.c1dvmp.c1dvmp(th_dvmp, 1, (0.5+1.j), 2)
    # comparing to DM's DVEM-c1_forKreso.nb
    # quark part
    assert aux[0] == approx(30.3836+9.2856j,  rel=1e-5)
    # "pure singlet" part
    assert aux[1] == approx(-0.496221+3.85004j, rel=1e-5)
    # gluon part
    # DM's notebook value is a permil away ...
    assert aux[2] == approx(35.2725+34.0699j,  rel=1e-5)
    # ... from formula from "Towards DVMP" paper:
    # assert aux[2] == approx(35.2818+34.0699j, rel=1e-5)


@mark.slow
def test_dvmp_TFFs_NLO(th_nlo):
    """Calculate NLO DVMP TFFs for rho0 production at input scale."""
    xB = 1e-4
    pt = g.data.DataPoint({'Q2': 4., 't': 0, 'xB': xB})
    # to get agreement with these old numbers:
    th_nlo.corr_c1dvmp_one = 0
    th_nlo.corr_c1dvmp_sgn = -1
    tffs = th_nlo.tff(pt.xi, pt.t, pt.Q2)
    reh, imh = tffs[0], tffs[1]
    # following agrees with gepard-fortran ...
    assert reh == approx(5410.6143, rel=1e-5)
    assert imh == approx(22869.412, rel=1e-5)
    # which for unknown reason regressed from situation:
    # following agrees with DM to best than percent
    # assert reh == approx(5328.3678, rel=1e-3)
    # assert imh == approx(22676.063, rel=1e-3)
