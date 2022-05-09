"""Testing code for MINUIT fitting."""

import gepard as g
from pytest import approx, fixture, mark

# The following is not resulting in paralelization
# os.environ["OPENBLAS_MAIN_FREE"] = "1"

class FitTest(g.gpd.PWNormGPD, g.cff.MellinBarnesCFF, g.dvcs.BMK):
    pass

class FitTestDIS(g.gpd.PWNormGPD, g.dis.DIS):
    pass

class FitTestHot(g.gpd.PWNormGPD, g.cff.MellinBarnesCFF, g.eff.DipoleEFF, g.dvcs.hotfixedBMK):
    pass

@fixture
def th_dis_lo():
    th = FitTestDIS(p=0)
    th.m.parameters.update({'ns': 0.15, 'al0s': 1., 'alps': 0.15, 'ms2': 1.,
                            'secs': 0., 'al0g': 1.1, 'alpg': 0.15, 'mg2': 0.7})
    return th

@fixture
def th_dis_nlo():
    th = FitTestDIS(p=1)
    th.m.parameters.update({'ns': 0.15, 'al0s': 1., 'alps': 0.15, 'ms2': 1.,
                            'secs': 0., 'al0g': 1.1, 'alpg': 0.15, 'mg2': 0.7})
    return th

@fixture
def th():
    th = FitTest(p=0)
    th.m.parameters.update({'ns': 0.15203911208796006,
                            'al0s': 1.1575060246398083,
                            'alps': 0.15,
                            'ms2': 1.,
                            'secs': 0.,
                            'al0g': 1.247316701070471,
                            'alpg': 0.15,
                            'mg2': 0.7})
    return th

@fixture
def th_nlo():
    th = FitTestHot(p=1, scheme='msbar')
    th.m.parameters.update({'ns': 0.1678,
                            'al0s': 1.12835,
                            'alps': 0.15,
                            'ms2': 1.,
                            'secs': 0.,
                            'al0g': 1.099,
                            'alpg': 0.15,
                            'mg2': 0.7})
    return th

@mark.slow
def test_fit_DIS_LO(th_dis_lo):
    """Test LO fitting to HERA DIS F2 data."""
    DISpoints = g.dset[201].copy()
    for id in range(202, 213):
        DISpoints += g.dset[id]
    f = g.fitter.MinuitFitter(DISpoints, th_dis_lo)
    f.release_parameters('ns', 'al0s', 'al0g')
    f.minuit.migrad()
    assert f.minuit.fval == approx(49.731197679982)
    assert th_dis_lo.parameters['al0s'] == approx(1.1575179114289278)


@mark.slow
def test_fit_DIS_NLO(th_dis_nlo):
    """Test NLO fitting to HERA DIS F2 data."""
    DISpoints = g.dset[201].copy()
    for id in range(202, 213):
        DISpoints += g.dset[id]
    f = g.fitter.MinuitFitter(DISpoints, th_dis_nlo)
    f.release_parameters('ns', 'al0s', 'al0g')
    f.minuit.migrad()
    assert f.minuit.fval == approx(71.6184113449641)
    assert th_dis_nlo.parameters['al0s'] == approx(1.1283626215013125)


def test_gepardfitDVCSnlso3_short(th):
    """Test fitting nl-so3 model to HERA DVCS data.

    This is reduced faster version of old test below.
    """
    # To pre-calculate wce for parallel execution:
    th.chisq_single(g.dset[39])
    f = g.fitter.MinuitFitter(g.dset[39], th)
    f.release_parameters('ms2')
    f.minuit.migrad()
    assert f.minuit.fval == approx(19.69615, rel=1.e-2)


@mark.slow
def test_errorprop(th):
    """Test uncertainty propagation to CFF."""
    f = g.fitter.MinuitFitter(g.dset[36], th)
    f.release_parameters('ns', 'ms2', 'secs')
    f.fit()
    res = th.predict(g.dset[36][0], observable='ImH', uncertainty=True)
    # number from gepard-fortran
    assert res[1] == approx(187.684, rel=1.e-3)


@mark.slow
def test_gepardfitDVCSnlso3_lo_long(th):
    """Test fitting nl-SO3 LO model to HERA DVCS data.

    This should give same results as in gepard's
    'fit dvcs dvcs dvcs' or
    smallx-final.nb, section 1-[nlo]-LO.
    """
    # To pre-calculate wce for parallel execution:
    DVCSpoints = g.dset[36].copy()
    for id in range(37, 46):
        DVCSpoints += g.dset[id]
    th.chisq_single(DVCSpoints)
    f = g.fitter.MinuitFitter(DVCSpoints, th)
    f.release_parameters('ms2', 'secs', 'secg')
    f.minuit.migrad()
    # These are old value from NPB paper days:
    # assert th.chisq(f.fitpoints) == approx(95.92, rel=1.e-2)
    # assert th.parameters['ms2'] == approx(0.47839, rel=1e-2)
    # assert th.parameters['secs'] == approx(-0.15152, rel=1e-2)
    # assert th.parameters['secg'] == approx(-0.81216, rel=1e-2)
    # Today's values, from pyfortran 2022 branch:
    assert th.chisq(f.fitpoints) == approx(96.52722, rel=1.e-4)
    assert th.parameters['ms2'] == approx(0.477888357, rel=1e-4)
    assert th.parameters['secs'] == approx(-0.1513351, rel=1e-4)
    assert th.parameters['secg'] == approx(-0.8146066876, rel=1e-4)


@mark.slow
def test_gepardfitDVCSnlso3_nlo_long(th_nlo):
    """Test fitting nl-SO3 NLO model to HERA DVCS data.

    This should give same results as in
    smallx-final.nb, section 1-[nlo]-NLO MSBAR.
    """
    # To pre-calculate wce for parallel execution:
    DVCSpoints = g.dset[36].copy()
    for id in range(37, 46):
        DVCSpoints += g.dset[id]
    th_nlo.chisq_single(DVCSpoints)
    f = g.fitter.MinuitFitter(DVCSpoints, th_nlo)
    f.release_parameters('ms2', 'secs', 'secg')
    f.minuit.migrad()
    # These are old value from NPB paper days:
    # assert th_nlo.chisq(f.fitpoints) == approx(101.581, rel=1e-2)
    # assert th_nlo.parameters['ms2'] == approx(0.590339, rel=1e-2)
    # Today's values, from pyfortran 2022 branch:
    assert th_nlo.chisq(f.fitpoints) == approx(104.8157642, rel=1e-4)
    assert th_nlo.parameters['ms2'] == approx(0.586664, rel=1e-4)
