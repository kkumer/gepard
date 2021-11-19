"""Testing code for observables."""

import sys

import gepard as g
from pytest import approx, mark

sys.path.append('/home/kkumer/g')
sys.path.append('/home/kkumer/g/gepard')

par_test = {'ns': 2./3. - 0.4, 'al0s': 1.1, 'alps': 0.25, 'ms2': 1.1**2,
            'ng': 0.4, 'al0g': 1.2, 'alpg': 0.25, 'mg2': 1.2**2}

par_fit = {'ns':  0.152039, 'al0s': 1.15751, 'alps': 0.15, 'ms2': 0.478391,
           'secs': -0.15152, 'this': 0.,  # 'ng': 0.4,  # provided by ns
           'al0g': 1.24732, 'alpg': 0.15, 'mg2': 0.7, 'secg': -0.81217, 'thig': 0.}

HERAtestpts = g.data.dset[39] + g.data.dset[45]
pt6 = g.data.dset[36][0]


def test_F2_NLO():
    """Test NLO DIS F2 evaluation."""
    fit_gpd = g.model.Fit(p=1)
    m = g.model.MellinBarnesModel(gpds=fit_gpd)
    th = g.theory.BMK(model=m)
    th.m.parameters.update({'ns': 0.15, 'al0s': 1., 'alps': 0.15, 'ms2': 1.,
                            'secs': 0., 'al0g': 1.1, 'alpg': 0.15, 'mg2': 0.7})
    f2 = th.F2(g.data.dset[201][0])
    assert f2 == approx(0.4220514008395502)


def test_XDVCSt_noevol():
    """Calculate LO DVCS partial cross section (no evolution)."""
    pt = g.data.DataPoint({'W': 82., 'Q2': 1., 't': 0.})
    pt.xi = pt.Q2 / (2.0 * pt.W * pt.W + pt.Q2)
    pt.xB = 2*pt.xi/(1.+pt.xi)
    pt.yaxis = 'X'
    test_gpd = g.model.Test()
    m = g.model.MellinBarnesModel(gpds=test_gpd)
    m.parameters.update(par_test)
    th = g.theory.hotfixedBMK(m)
    # Have to decrease precision due to slightly different formula
    assert th.predict(pt) == approx(5608.42804256, rel=1e-3)


def test_XDVCSt():
    """Calculate LO DVCS partial cross section (+ evolution)."""
    pt = g.data.DataPoint({'W': 82., 'Q2': 3., 't': -0.5})
    pt.xi = pt.Q2 / (2.0 * pt.W * pt.W + pt.Q2)
    pt.xB = 2*pt.xi/(1.+pt.xi)
    pt.yaxis = 'X'
    test_gpd = g.model.Test()
    m = g.model.MellinBarnesModel(gpds=test_gpd)
    m.parameters.update(par_test)
    th = g.theory.hotfixedBMK(m)
    # Comparison to old pure-Fortran 'src/test/test.F':
    # assert th._XDVCStApprox(pt) == approx(900.11627295345193)  # for t=0
    assert th._XDVCStApprox(pt) == approx(15.633422291049154)


def test_XDVCSt_NLO():
    """Calculate NLO DVCS partial cross section (no evolution)."""
    pt = g.data.DataPoint({'W': 82., 'Q2': 1., 't': 0.})
    pt.xi = pt.Q2 / (2.0 * pt.W * pt.W + pt.Q2)
    pt.xB = 2*pt.xi/(1.+pt.xi)
    pt.yaxis = 'X'
    test_gpd = g.model.Test(p=1)
    m = g.model.MellinBarnesModel(gpds=test_gpd)
    m.parameters.update(par_test)
    th = g.theory.hotfixedBMK(m)
    assert th._XDVCStEx(pt) == approx(821.0062045181508)


def test_XDVCSt_NLOevol():
    """Calculate NLO DVCS partial cross section (+ evolution)."""
    pt = g.data.DataPoint({'W': 82., 'Q2': 8., 't': 0.})
    pt.xi = pt.Q2 / (2.0 * pt.W * pt.W + pt.Q2)
    pt.xB = 2*pt.xi/(1.+pt.xi)
    pt.yaxis = 'X'
    test_gpd = g.model.Test(p=1)
    m = g.model.MellinBarnesModel(gpds=test_gpd)
    m.parameters.update(par_test)
    th = g.theory.hotfixedBMK(m)
    assert th._XDVCStEx(pt) == approx(90.76770897337423)


def test_XDVCS():
    """Calculate LO DVCS total cross section (+ evolution)."""
    pt = g.data.DataPoint({'W': 82., 'Q2': 3.})
    pt.xi = pt.Q2 / (2.0 * pt.W * pt.W + pt.Q2)
    pt.xB = 2*pt.xi/(1.+pt.xi)
    pt.yaxis = 'X'
    test_gpd = g.model.Test()
    m = g.model.MellinBarnesModel(gpds=test_gpd)
    m.parameters.update(par_test)
    th = g.theory.hotfixedBMK(m)
    # Comparison to old pure-Fortran 'src/test/test.F':
    # Have to decrease precision due to slightly different formula
    assert th.predict(pt) == approx(105.18391660404916, rel=1e-3)


def test_XDVCS_NLO():
    """Calculate NLO DVCS cross section (+ evolution)."""
    pt = g.data.DataPoint({'W': 55., 'Q2': 3.})
    pt.xi = pt.Q2 / (2.0 * pt.W * pt.W + pt.Q2)
    pt.xB = 2*pt.xi/(1.+pt.xi)
    pt.yaxis = 'X'
    test_gpd = g.model.Test(p=1)
    m = g.model.MellinBarnesModel(gpds=test_gpd)
    m.parameters.update(par_test)
    th = g.theory.hotfixedBMK(m)
    assert th.predict(pt) == approx(32.29728102, rel=1e-5)
    # To get complete agreement with Fortran take
    # (slower) tquadrature = quadSciPy10 and:
    # set pt.W=3.5  and use XDVCStApprox and
    # assert_almost_equal(aux, 6.8612469682766850, 5)


def test_predict():
    """Test theory predict."""
    fit_gpd = g.model.Fit()
    m = g.model.MellinBarnesModel(gpds=fit_gpd)
    th = g.theory.BMK(model=m)
    th.m.parameters.update(par_fit)
    assert th.predict(pt6) == approx(12.69069206084793)

# @mark.devel
# def test_chisq():
#     """Test chisq calculation."""
#     fit_gpd = g.model.Fit()
#     m = g.model.MellinBarnesModel(gpds=fit_gpd)
#     th = g.theory.BMK(model=m)
#     th.m.parameters.update(par_fit)
#     assert th.chisq_single([data[39][1]]) == approx(0.027310271896327565)
#     assert th.chisq_para([data[39][1]]) == approx(0.027310271896327565)


def test_chisq_Xt():
    """Test chisq Xt calculation."""
    fit_gpd = g.model.Fit()
    m = g.model.MellinBarnesModel(gpds=fit_gpd)
    th = g.theory.BMK(model=m)
    th.m.parameters.update(par_fit)
    assert th.chisq_single(g.data.dset[39][2:4]) == approx(1.012533085207716)
    # parallelization slows down fast testing
    # assert th.chisq_para(g.data.dset[39][2:4]) == approx(1.012533085207716)


def test_chisq_X():
    """Test chisq total X calculation."""
    fit_gpd = g.model.Fit()
    m = g.model.MellinBarnesModel(gpds=fit_gpd)
    th = g.theory.BMK(model=m)
    th.m.parameters.update(par_fit)
    assert th.chisq_single([g.data.dset[45][3]]) == approx(0.29405520100706245)
    # parallelization slows down fast testing
    # assert th.chisq_para([g.data.dset[45][3]]) == approx(0.29405520100706245)


# def test_chisq_aux():
#     """Test chisq total X calculation."""
#     fit_gpd = g.model.Fit()
#     m = g.model.MellinBarnesModel(gpds=fit_gpd)
#     th = g.theory.BMK(model=m)
#     th.m.parameters.update(par_fit)
#     # assert th.chisq_single(data[45][3:5]) == approx(0.29860340473733393)
#     # assert th.chisq_para(data[45][3:5]) == approx(0.29860340473733393)
#     # assert th.chisq_single(data[45][2:4]) == approx(0.4897683229894142)
#     # assert th.chisq_para(data[45][2:4]) == approx(0.4897683229894142)
#     # assert th.chisq_single(data[39][2:4]+data[45][2:4]) == approx(1.502301408197130)
#     # assert th.chisq_para(data[39][2:4]+data[45][2:4]) == approx(1.502301408197130)
#     # assert th.chisq_single(data[39][2:4]+data[45][2:5]) == approx(1.5068496119274017)
#     # assert th.chisq_para(data[39][2:4]+data[45][2:5]) == approx(1.5068496119274017)
#     # assert th.chisq_single(data[39][2:4]+data[45][2:]) == approx(2.23723298841219)
#     # assert th.chisq_para(data[39][2:4]+data[45][2:]) == approx(2.23723298841219)
#     # assert th.chisq_single(data[39][2:]+data[45][2:]) == approx(3.9535384269614986)
#     # assert th.chisq_para(data[39][2:]+data[45][2:]) == approx(3.9535384269614986)
#     assert th.chisq_single(data[39][2:]+data[45]) == approx(10.4732692425849)
#     assert th.chisq_para(data[39][2:]+data[45]) == approx(10.4732692425849)


def test_chisq_XtX():
    """Test chisq Xt+X calculation."""
    fit_gpd = g.model.Fit()
    m = g.model.MellinBarnesModel(gpds=fit_gpd)
    th = g.theory.BMK(model=m)
    th.m.parameters.update(par_fit)
    assert th.chisq_single(g.data.dset[39][2:]+g.data.dset[45]) == approx(10.47326924258)
    # assert th.chisq_para(data[39][2:]+data[45]) == approx(10.47326924258)
