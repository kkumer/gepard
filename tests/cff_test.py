"""Testing code for CFFs."""

import sys

from pytest import approx, mark

import gepard as g

sys.path.append('/home/kkumer/g')
sys.path.append('/home/kkumer/g/gepard')

par_test = {'ns': 2./3. - 0.4, 'al0s': 1.1, 'alps': 0.25, 'ms': 1.1,
            'ng': 0.4, 'al0g': 1.2, 'alpg': 0.25, 'mg': 1.2}

par_fit = {'ns':  0.152039, 'al0s': 1.15751, 'alps': 0.15, 'ms': 0.478391,
           'secs': -0.15152, 'this': 0.,  # 'ng': 0.4,  # provided by ns
           'al0g': 1.24732, 'alpg': 0.15, 'mg': 0.7, 'secg': -0.81217, 'thig': 0.}

pt_test = g.data.DummyPoint({'W': 82., 'Q2': 1., 't': 0.})
pt_test.xi = pt_test.Q2 / (2.0 * pt_test.W * pt_test.W + pt_test.Q2)

pt0_fit = g.data.DummyPoint({'xi': 0.01, 'Q2': 4., 't': -0.2})  # noevol
pt_fit = g.data.DummyPoint({'xi': 0.01, 'Q2': 8., 't': -0.2})  # evol


def test_cff_H_test_noevol():
    """Test testing (ReH, ImH) (no evolution)."""
    test_gpd = g.model.Test()
    m_test = g.model.MellinBarnesModel(gpds=test_gpd)
    m_test.parameters.update(par_test)
    assert m_test.cff(pt_test.xi, pt_test.t, pt_test.Q2)[:2] == approx(
            [9839.566, 61614.9])


def test_cff_H_nlso3_noevol():
    """Test nl-so3  (ReH, ImH) (no evolution)."""
    fit_gpd = g.model.Fit()
    m_fit = g.model.MellinBarnesModel(gpds=fit_gpd)
    m_fit.parameters.update(par_fit)
    assert m_fit.cff(pt0_fit.xi, pt0_fit.t, pt0_fit.Q2)[:2] == approx(
            [18.9540, 60.1622], abs=0.001)


def test_cff_H_nlso3():
    """Test nl-so3 (ReH, ImH) (LO evolved to multiple Q2)."""
    fit_gpd = g.model.Fit()
    m_fit = g.model.MellinBarnesModel(gpds=fit_gpd)
    m_fit.parameters.update(par_fit)
    # Q2 can be changed during calls:
    assert m_fit.cff(pt0_fit.xi, pt0_fit.t, pt0_fit.Q2)[:2] == approx(
            [18.9540, 60.1622], abs=0.001)
    assert m_fit.cff(pt_fit.xi, pt_fit.t, pt_fit.Q2)[:2] == approx(
            (26.9984, 66.5255), abs=0.001)


# def test_cff_E_fit_noevol():
    # """Test fitting (ReE, ImE) (no evolution)."""
    # assert m_fit.cff(pt_fit.xi, pt_fit.t, pt_fit.Q2)[:2] == approx(
            # (13.2678, 42.11355), abs=0.001)


# def test_cff_E_fit():
    # """Test fitting (ReE, ImE) (LO evolved to Q2 = 8 GeV^2)."""
    # assert m_fit.cff(pt_fit.xi, pt_fit.t, pt_fit.Q2)[:2] == approx(
            # (13.5297, 43.6024), abs=0.001)
