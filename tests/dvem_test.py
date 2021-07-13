"""Tests for DVEM."""

import gepard as g
from pytest import approx, mark

par_dvmp = {'ns':  0.152, 'al0s': 1.158, 'alps': 0.15, 'ms': 0.446,
            'secs': -0.442, 'this': 0.089,  # 'ng': 0.5,  # provided by ns
            'al0g': 1.247, 'alpg': 0.15, 'mg': 0.7, 'secg': -2.309, 'thig': 0.812}


def test_dvem_TFFs_LO():
    """Calculate LO DVMP TFFs for rho production at input scale."""
    xB = 1e-4
    pt = g.data.DataPoint({'Q2': 4., 't': 0, 'xB': xB})
    # generic LO model from big DVMP draft
    fit_gpd = g.model.Fit()
    m = g.model.MellinBarnesModel(gpds=fit_gpd)
    m.parameters.update(par_dvmp)
    th = g.theory.BMK(m)
    tffs = th.m.tff(pt.xi, pt.t, pt.Q2)
    reh, imh = tffs[0], tffs[1]
    # following agrees with DM to best than percent
    assert imh == approx(12395.53, rel=1e-3)
    assert reh == approx(4766.8993, rel=1e-3)


@mark.skip('KM10b not yet implemented')
def test_gepardTFFsEvol():
    """Calculate LO DVMP TFFs for rho production + evolution."""
    pt = g.data.DataPoint({'Q2': 6.6, 'W': 75., 't': -0.025})
    fit_gpd = g.model.Fit()
    m = g.model.MellinBarnesModel(gpds=fit_gpd)
    m.parameters.update(par_KM10b)
    th = g.theory.BMK(m)
    tffs = th.m.tff(pt.xi, pt.t, pt.Q2)
    reh, imh = tffs[0], tffs[1]
    assert reh == approx(285.15512, rel=1e-2)
    assert imh == approx(511.39622404, rel=1e-2)


@mark.skip('KM10b not yet implemented')
def test_gepardXrhot():
    """Calculate LO DVMP cross section d sigma / dt"""
    pt = g.data.DataPoint({'Q2': 6.6, 'W': 75., 't': -0.025,
                           'process': 'gammastarp2rho0p'})
    fit_gpd = g.model.Fit()
    m = g.model.MellinBarnesModel(gpds=fit_gpd)
    m.parameters.update(par_KM10b)
    th = g.theory.BMK(m)
    assert th.X(pt) == approx(1212.62165, rel=1.e-2)


def test_c1_NLO():
    """Calculate NLO DVMP hard scattering coefficients"""
    xB = 1e-4
    pt = g.data.DataPoint({'Q2': 4., 't': 0, 'xB': xB})
    # generic LO model from big DVMP draft
    fit_gpd = g.model.Fit()
    m = g.model.MellinBarnesModel(gpds=fit_gpd)
    m.parameters.update(par_dvmp)
    # c1dvmp(model, sgntr, j, k)
    aux = g.c1dvmp.c1dvmp(m, 1, (0.5+1.j), 2)
    # comparing to DM's DVEM-c1_forKreso.nb
    # quark part
    assert aux[0] == approx(30.3836+9.2856j,  rel=1e-5)
    # "pure singlet" part
    assert aux[1] == approx(-0.496221+3.85004j, rel=1e-5)
    # gluon part
    assert aux[2] == approx(35.2725+34.0699j,  rel=1e-5)


@mark.skip('NLO fails - needs careful check')
def test_dvem_TFFs_NLO():
    """Calculate NLO DVMP TFFs for rho production at input scale."""
    xB = 1e-4
    pt = g.data.DataPoint({'Q2': 4., 't': 0, 'xB': xB})
    fit_gpd = g.model.Fit(p=1)
    m = g.model.MellinBarnesModel(gpds=fit_gpd)
    m.parameters.update(par_dvmp)
    th = g.theory.BMK(m)
    re, im =  (th.m.ReHrho(pt), th.m.ImHrho(pt))
    tffs = th.m.tff(pt.xi, pt.t, pt.Q2)
    reh, imh = tffs[0], tffs[1]
    # following agrees with gepard-fortran ...
    # assert reh == approx(5410.6143, rel=1e-5)
    # assert imh == approx(22869.412, rel=1e-5)
    # which for unknown reason regressed from situation:
    # following agrees with DM to best than percent
    assert reh == approx(5328.3678, rel=1e-3)
    assert imh == approx(22676.063, rel=1e-3)
