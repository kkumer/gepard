"""Testing code for observables."""

import sys

import gepard as g
from pytest import approx, fixture, mark
from gepard.fits import th_KM15

sys.path.append('/home/kkumer/g')
sys.path.append('/home/kkumer/g/gepard')

par_test = {'ns': 2./3. - 0.4, 'al0s': 1.1, 'alps': 0.25, 'ms2': 1.1**2,
            'ng': 0.4, 'al0g': 1.2, 'alpg': 0.25, 'mg2': 1.2**2}

par_fit = {'ns':  0.152039, 'al0s': 1.15751, 'alps': 0.15, 'ms2': 0.478391,
           'secs': -0.15152, 'this': 0.,  # 'ng': 0.4,  # provided by ns
           'al0g': 1.24732, 'alpg': 0.15, 'mg2': 0.7, 'secg': -0.81217, 'thig': 0.}


class MyTest(g.gpd.PWNormGPD, g.cff.MellinBarnesCFF, g.dvcs.BMK):
    pass

class MyTestDIS(g.gpd.PWNormGPD, g.dis.DIS):
    pass

class MyTest2(g.gpd.TestGPD, g.cff.MellinBarnesCFF, g.dvcs.hotfixedBMK):
    pass

@fixture
def th_dis():
    th = MyTestDIS(p=1)
    th.parameters.update({'ns': 0.15, 'al0s': 1., 'alps': 0.15, 'ms2': 1.,
                            'secs': 0., 'al0g': 1.1, 'alpg': 0.15, 'mg2': 0.7})
    return th

@fixture
def thx():
    th = MyTest(p=0)
    th.parameters.update(par_fit)
    return th

@fixture
def th():
    th = MyTest(p=1)
    th.parameters.update(par_fit)
    return th

@fixture
def th_lo():
    th = MyTest2(p=0)
    th.parameters.update(par_test)
    return th

@fixture
def th_nlo():
    th = MyTest2(p=1)
    th.parameters.update(par_test)
    return th



def test_F2_NLO(th_dis):
    """Test NLO DIS F2 evaluation."""
    f2 = th_dis.DISF2(g.data.dset[201][0])
    assert f2 == approx(0.4220514008395502)


def test_XDVCSt_noevol(th_lo):
    """Calculate LO DVCS partial cross section (no evolution)."""
    pt = g.data.DataPoint({'W': 82., 'Q2': 1., 't': 0.})
    pt.xi = pt.Q2 / (2.0 * pt.W * pt.W + pt.Q2)
    pt.xB = 2*pt.xi/(1.+pt.xi)
    pt.observable = 'XGAMMA'
    # Have to decrease precision due to slightly different formula
    assert th_lo.predict(pt) == approx(5608.42804256, rel=1e-3)


def test_XDVCSt(th_lo):
    """Calculate LO DVCS partial cross section (+ evolution)."""
    pt = g.data.DataPoint({'W': 82., 'Q2': 3., 't': -0.5})
    pt.xi = pt.Q2 / (2.0 * pt.W * pt.W + pt.Q2)
    pt.xB = 2*pt.xi/(1.+pt.xi)
    pt.observable = 'XGAMMA'
    # Comparison to old pure-Fortran 'src/test/test.F':
    # assert th_lo._XGAMMA_DVCS_t_Approx(pt) == approx(900.11627295345193)  # for t=0
    assert th_lo._XGAMMA_DVCS_t_Approx(pt) == approx(15.633422291049154)


def test_XDVCSt_NLO(th_nlo):
    """Calculate NLO DVCS partial cross section (no evolution)."""
    pt = g.data.DataPoint({'W': 82., 'Q2': 1., 't': 0.})
    pt.xi = pt.Q2 / (2.0 * pt.W * pt.W + pt.Q2)
    pt.xB = 2*pt.xi/(1.+pt.xi)
    pt.observable = 'XGAMMA'
    assert th_nlo._XGAMMA_DVCS_t_Ex(pt) == approx(821.0062045181508)


def test_XDVCSt_NLOevol(th_nlo):
    """Calculate NLO DVCS partial cross section (+ evolution)."""
    pt = g.data.DataPoint({'W': 82., 'Q2': 8., 't': 0.})
    pt.xi = pt.Q2 / (2.0 * pt.W * pt.W + pt.Q2)
    pt.xB = 2*pt.xi/(1.+pt.xi)
    pt.observable = 'XGAMMA'
    assert th_nlo._XGAMMA_DVCS_t_Ex(pt) == approx(90.76770897337423)


def test_XDVCS(th_lo):
    """Calculate LO DVCS total cross section (+ evolution)."""
    pt = g.data.DataPoint({'W': 82., 'Q2': 3.})
    pt.xi = pt.Q2 / (2.0 * pt.W * pt.W + pt.Q2)
    pt.xB = 2*pt.xi/(1.+pt.xi)
    pt.observable = 'XGAMMA'
    # Comparison to old pure-Fortran 'src/test/test.F':
    # Have to decrease precision due to slightly different formula
    assert th_lo.predict(pt) == approx(105.18391660404916, rel=1e-3)


def test_XDVCS_NLO(th_nlo):
    """Calculate NLO DVCS cross section (+ evolution)."""
    pt = g.data.DataPoint({'W': 55., 'Q2': 3.})
    pt.xi = pt.Q2 / (2.0 * pt.W * pt.W + pt.Q2)
    pt.xB = 2*pt.xi/(1.+pt.xi)
    pt.observable = 'XGAMMA'
    assert th_nlo.predict(pt) == approx(32.29728102, rel=1e-5)
    # To get complete agreement with Fortran take
    # (slower) tquadrature = quadSciPy10 and:
    # set pt.W=3.5  and use XDVCStApprox and
    # assert_almost_equal(aux, 6.8612469682766850, 5)


def test_predict(thx):
    """Test theory predict."""
    pt6 = g.data.dset[36][0]
    assert thx.predict(pt6) == approx(12.69069206084793)

# @mark.devel
# def test_chisq():
#     """Test chisq calculation."""
#     assert th.chisq_single([data[39][1]]) == approx(0.027310271896327565)
#     assert th.chisq_para([data[39][1]]) == approx(0.027310271896327565)


def test_chisq_Xt(thx):
    """Test chisq Xt calculation."""
    assert thx.chisq_single(g.data.dset[39][2:4]) == approx(1.012533085207716)
    # parallelization slows down fast testing
    # assert th.chisq_para(g.data.dset[39][2:4]) == approx(1.012533085207716)


def test_chisq_X(thx):
    """Test chisq total X calculation."""
    assert thx.chisq_single([g.data.dset[45][3]]) == approx(0.29405520100706245)
    # parallelization slows down fast testing
    # assert th.chisq_para([g.data.dset[45][3]]) == approx(0.29405520100706245)


# def test_chisq_aux():
#     """Test chisq total X calculation."""
#     fit_gpd = g.gpd.PWNormGPD()
#     m = g.cff.MellinBarnesCFF(gpds=fit_gpd)
#     th = g.dvcs.BMK(model=m)
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


def test_chisq_XtX(thx):
    """Test chisq Xt+X calculation."""
    assert thx.chisq_single(g.data.dset[39][2:]+g.data.dset[45]) == approx(10.47326924258)
    # assert th.chisq_para(data[39][2:]+data[45]) == approx(10.47326924258)


@mark.skip('Uncertainties not yet implemented.')
def test_uncert():
    """Test uncertainty of prediction."""
    assert th_KM15.predict(g.dset[45][0], uncertainty=True) == (5.195327906304886, 0.2081213421895121)
