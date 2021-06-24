"""Tests for DVEM."""

import gepard as g
from pytest import approx, mark

par_dvmp = {'ns':  0.152, 'al0s': 1.158, 'alps': 0.15, 'ms': 0.446,
            'secs': -0.442, 'this': 0.089,  # 'ng': 0.5,  # provided by ns
            'al0g': 1.247, 'alpg': 0.15, 'mg': 0.7, 'secg': -2.309, 'thig': 0.812}


def test_gepardTFFs():
    """Calculate LO DVMP TFFs for rho production at input scale."""
    xB = 1e-4
    pt = g.data.DataPoint({'Q2': 4., 't': 0, 'xi': xB/2., 'xB': xB})
    g.utils.fill_kinematics(pt)
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

