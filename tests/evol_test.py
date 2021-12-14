"""Testing evolution code.

These tests are not necessary. They have been used
for development of evolution operator. Could be marked
"extended testing" or something.
"""

import gepard as g
import numpy as np
from pytest import approx, fixture

par_test = {'ns': 2./3. - 0.4, 'al0s': 1.1, 'alps': 0.25, 'ms2': 1.1**2,
            'ng': 0.4, 'al0g': 1.2, 'alpg': 0.25, 'mg2': 1.2**2}

par_fit = {'ns':  0.152039, 'al0s': 1.15751, 'alps': 0.15, 'ms2': 0.478391,
           'secs': -0.15152, 'this': 0.,  # 'ng': 0.4,  # provided by ns
           'al0g': 1.24732, 'alpg': 0.15, 'mg2': 0.7, 'secg': -0.81217, 'thig': 0.,
           'kaps': 0.7, 'kapg': -0.2}

pt0_fit = g.data.DataPoint({'xi': 0.01, 'Q2': 4., 't': -0.2})  # noevol
pt_fit = g.data.DataPoint({'xi': 0.01, 'Q2': 8., 't': -0.2})  # evol

class MyTest(g.gpd.PWNormGPD, g.cff.MellinBarnesCFF):
    pass

class MyTest2(g.gpd.TestGPD, g.cff.MellinBarnesCFF):
    pass

@fixture
def th():
    th = MyTest()
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



def test_lambda(th):
    """Test LO singlet an. dim. eigenvalues."""
    # leading PW
    gam0 = g.adim.singlet_LO(th.jpoints+1, th.nf).transpose((2, 0, 1))
    assert g.evolution.lambdaf(gam0)[:, 0] == approx(
           np.array([-22.79064075+0.01967868j, 3.56546819+0.00069647j]))
    # nl PW
    gam0 = g.adim.singlet_LO(th.jpoints+3, th.nf).transpose((2, 0, 1))
    assert g.evolution.lambdaf(gam0)[:, 0] == approx(
           np.array([13.07805098+0.00080784j, 5.78554055+0.00036216j]))


def test_projectors_LO(th):
    """Test LO singlet eigen projectors."""
    gam0 = g.adim.singlet_LO(th.jpoints+1, th.nf).transpose((2, 0, 1))
    lam, pr = g.evolution.projectors(gam0)
    assert pr[0, 0, :, :] == approx(
           np.array([[0.07529713+5.20367327e-05j, 0.14772791+8.44139525e-05j],
                     [0.47132239+2.98800251e-05j, 0.92470287-5.20367327e-05j]]))


def test_rnlof(th_nlo):
    """Test projected NLO mu-indep part of evol.op."""
    # leading PW
    lam, pr, r1proj = g.evolution.rnlof(th_nlo, th_nlo.jpoints)
    assert pr[0, 0, :, :] == approx(
           np.array([[5.68471073518716855e-02+3.94946205837689893e-05j,
                      0.11226596856752971+6.51645223463349292e-05j],
                     [0.47757586947844916+3.45902444197994830e-05j,
                      0.94315289264812829-3.94946205837689961e-05j]]))
    assert r1proj[0, 0, 1, :, :] == approx(
           np.array([[0.25405852938085721-2.83493819853408962e-05j,
                      -3.02412560396950227e-02-1.54453206209269115e-05j],
                     [2.1343594746465175-1.56642538048529668e-03j,
                      -0.25405852938085760+2.83493819853415738e-05j]]))


def test_evolop_LO(th_lo):
    """Test LO evolution operator."""
    assert g.evolution.evolop(th_lo, th_lo.jpoints, 3.0,
                              'DVCS')[0, 0, :, :] == approx(
           np.array([[0.97360574083605833-4.26361786894205834e-05j,
                      0.12173639863278003-5.99655383745874504e-05j],
                     [0.51786249629889658-5.18175412363711162e-04j,
                      1.9346770097724764-1.15955006729918713e-03j]]))


def test_evolop_NLO(th_nlo):
    """Test NLO evolution operator."""
    # LO part (but with NLO alpha_strong)
    assert g.evolution.evolop(th_nlo, th_nlo.jpoints, 3.0,
                              'DVCS')[0, 0, :, :] == approx(
           np.array([[0.97506856774185890-6.13627841862173179e-05j,
                      0.16175930125082716-9.38626466680329819e-05j],
                     [0.68811853114565436-7.48865910115848387e-04j,
                      2.2521082168110360-1.65744690707721643e-03j]]))
    # NLO part
    assert g.evolution.evolop(th_nlo, th_nlo.jpoints, 3.0,
                              'DVCS')[0, 1, :, :] == approx(
           np.array([[1.3209647413503760-2.00187561001276184e-03j,
                      3.0958151106827598-4.13038076850684704e-03j],
                     [-2.7796825098769098+5.59500993233383315e-03j,
                      -8.6121807772420294+1.50792240021393256e-02j]]))
