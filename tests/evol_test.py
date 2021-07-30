"""Testing evolution code.

These tests are not necessary. They have been used
for development of evolution operator. Could be marked
"extended testing" or something.
"""

import gepard as g
import numpy as np
from pytest import approx

par_fit = {'ns':  0.152039, 'al0s': 1.15751, 'alps': 0.15, 'ms': 0.478391,
           'secs': -0.15152, 'this': 0.,  # 'ng': 0.4,  # provided by ns
           'al0g': 1.24732, 'alpg': 0.15, 'mg': 0.7, 'secg': -0.81217, 'thig': 0.,
           'kaps': 0.7, 'kapg': -0.2}

pt0_fit = g.data.DataPoint({'xi': 0.01, 'Q2': 4., 't': -0.2})  # noevol
pt_fit = g.data.DataPoint({'xi': 0.01, 'Q2': 8., 't': -0.2})  # evol


def test_evol():
    """Test LO evolution operator."""
    fit_gpd = g.model.Fit()
    m_fit = g.model.MellinBarnesModel(gpds=fit_gpd)
    m_fit.parameters.update(par_fit)
    # leading PW
    assert g.evolution.lambdaf(fit_gpd, fit_gpd.jpoints)[:, 0] == approx(
           np.array([-22.79064075+0.01967868j, 3.56546819+0.00069647j]))
    # nl PW
    assert g.evolution.lambdaf(fit_gpd, fit_gpd.jpoints+2)[:, 0] == approx(
           np.array([13.07805098+0.00080784j, 5.78554055+0.00036216j]))


def test_rnnlof():
    """Test LO singlet eigen projectors."""
    fit_gpd = g.model.Fit()
    m_fit = g.model.MellinBarnesModel(gpds=fit_gpd)
    m_fit.parameters.update(par_fit)
    # leading PW
    pr, r1proj = g.evolution.rnnlof(fit_gpd, fit_gpd.jpoints)
    assert pr[0, 0, :, :] == approx(
           np.array([[0.07529713+5.20367327e-05j, 0.14772791+8.44139525e-05j],
                     [0.47132239+2.98800251e-05j, 0.92470287-5.20367327e-05j]]))
