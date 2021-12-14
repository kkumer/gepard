"""Tests for DVEM."""

import gepard as g
from pytest import approx, fixture, mark

par_dvmp = {'ns':  0.152, 'al0s': 1.158, 'alps': 0.15, 'ms2': 0.446,
            'secs': -0.442, 'this': 0.089,  # 'ng': 0.5,  # provided by ns
            'al0g': 1.247, 'alpg': 0.15, 'mg2': 0.7, 'secg': -2.309, 'thig': 0.812}

par_KM10b = {'tMv': 0.8, 'rS': 1.0, 'rpi': 4.0201, 'alv': 0.43, 'Nsea': 0.0,
             'Nv': 1.35, 'rv': 0.8081, 'Mpi': 1.5369, 'alS': 1.13, 'alpS': 0.15,
             'C': 5.4259, 'tNv': 0.6, 'bS': 2.0, 'bv': 0.7706, 'Mv': 0.8,
             'tbv': 1.0, 'alpv': 0.85, 'MC': 1.3305, 'MS': 0.707, 'trv': 3.2931,
             'EAL0G': 1.1, 'ESECS': 0.0, 'EDELM2S': 0.0, 'EPS': 2.0, 'ETHIS': 0.0,
             'ESECG': 0.0, 'EPG': 2.0, 'EDELM2G': 0.0, 'PS': 2.0, 'EALPG': 0.15,
             'EKAPG': 0.0, 'ESKEWG': 0.0, 'M02S': 0.49754317018981614,
             'EALPS': 0.15, 'EKAPS': 0.0, 'DELB': 0.0, 'ESKEWS': 0.0, 'SKEWS': 0.0,
             'ETHIG': 0.0, 'EM02G': 0.7, 'EAL0S': 1.0, 'DELM2S': 0.0, 'EM02S': 1.0,
             'SKEWG': 0.0, 'PG': 2.0, 'DELM2G': 0.0, 'ns': 0.15203911208796006,
             'al0s': 1.1575060246398083, 'alps': 0.15, 'al0g': 1.247316701070471,
             'secs': -0.4600511871918772, 'this': 0.09351798951979662,
             'alpg': 0.15, 'mg2': 0.7, 'secg': -2.5151319493485427,
             'ms2': 0.49754317018981614,
             'thig': 0.8915757559175185, 'kaps': 0.0, 'kapg': 0.0}

class MyTheory(g.gpd.PWNormGPD, g.cff.MellinBarnesCFF, g.theory.DVMP):
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
    """Calculate LO DVMP TFFs for rho production at input scale."""
    xB = 1e-4
    pt = g.data.DataPoint({'Q2': 4., 't': 0, 'xB': xB})
    tffs = th_dvmp.tff(pt.xi, pt.t, pt.Q2)
    reh, imh = tffs[0], tffs[1]
    # following agrees with DM to best than percent
    assert imh == approx(12395.53, rel=1e-7)
    assert reh == approx(4766.8993, rel=1e-7)


def test_gepardTFFsEvol(th_KM10b):
    """Calculate LO DVMP TFFs for rho production + evolution."""
    pt = g.data.DataPoint({'Q2': 6.6, 'W': 75., 't': -0.025})
    tffs = th_KM10b.tff(pt.xi, pt.t, pt.Q2)
    reh, imh = tffs[0], tffs[1]
    assert reh == approx(285.15512, rel=1e-2)
    assert imh == approx(511.39622404, rel=1e-2)


def test_gepardXrhot(th_KM10b):
    """Calculate LO DVMP cross section d sigma / dt"""
    pt = g.data.DataPoint({'Q2': 6.6, 'W': 75., 't': -0.025,
                           'process': 'gammastarp2rho0p'})
    #th = g.theory.DVMP(m)
    assert th_KM10b.X(pt) == approx(1212.62165, rel=1.e-2)


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
    """Calculate NLO DVMP TFFs for rho production at input scale."""
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
