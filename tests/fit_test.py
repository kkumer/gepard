"""Testing code for MINUIT fitting."""

import gepard as g
from pytest import approx, mark

# The following is not resulting in paralelization
# os.environ["OPENBLAS_MAIN_FREE"] = "1"

@mark.slow
def test_fit_DIS_LO():
    """Test LO fitting to HERA DIS F2 data."""
    fit_gpd = g.model.PWNormGPD(p=0)
    m = g.model.MellinBarnesModel(gpds=fit_gpd)
    th = g.theory.BMK(model=m)
    th.m.parameters.update({'ns': 0.15, 'al0s': 1., 'alps': 0.15, 'ms2': 1.,
                            'secs': 0., 'al0g': 1.1, 'alpg': 0.15, 'mg2': 0.7})
    DISpoints = g.dset[201].copy()
    for id in range(202, 213):
        DISpoints += g.dset[id]
    f = g.fitter.FitterMinuit(DISpoints, th)
    f.fix_parameters('ALL')
    f.release_parameters('ns', 'al0s', 'al0g')
    f.minuit.migrad()
    assert f.minuit.fval == approx(49.731197679982)
    assert th.m.parameters['al0s'] == approx(1.1575179114289278)


@mark.slow
def test_fit_DIS_NLO():
    """Test NLO fitting to HERA DIS F2 data."""
    fit_gpd = g.model.PWNormGPD(p=1)
    m = g.model.MellinBarnesModel(gpds=fit_gpd)
    th = g.theory.BMK(model=m)
    th.m.parameters.update({'ns': 0.15, 'al0s': 1., 'alps': 0.15, 'ms2': 1.,
                            'secs': 0., 'al0g': 1.1, 'alpg': 0.15, 'mg2': 0.7})
    DISpoints = g.dset[201].copy()
    for id in range(202, 213):
        DISpoints += g.dset[id]
    f = g.fitter.FitterMinuit(DISpoints, th)
    f.fix_parameters('ALL')
    f.release_parameters('ns', 'al0s', 'al0g')
    f.minuit.migrad()
    assert f.minuit.fval == approx(71.6184113449641)
    assert th.m.parameters['al0s'] == approx(1.1283626215013125)


def test_gepardfitDVCSnlso3_short():
    """Test fitting nl-so3 model to HERA DVCS data.

    This is reduced faster version of old test below.
    """
    fit_gpd = g.model.PWNormGPD()
    m = g.model.MellinBarnesModel(gpds=fit_gpd)
    th = g.theory.BMK(model=m)
    th.m.parameters.update({'ns': 0.15203911208796006,
                            'al0s': 1.1575060246398083,
                            'alps': 0.15,
                            'ms2': 1.,
                            'secs': 0.,
                            'al0g': 1.247316701070471,
                            'alpg': 0.15,
                            'mg2': 0.7})
    # To pre-calculate wce for parallel execution:
    th.chisq_single(g.dset[39])
    f = g.fitter.FitterMinuit(g.dset[39], th)
    f.fix_parameters('ALL')
    f.release_parameters('ms2')
    f.minuit.migrad()
    assert f.minuit.fval == approx(19.69615, rel=1.e-2)


@mark.slow
def test_gepardfitDVCSnlso3_long():
    """Test fitting nl-so3 model to HERA DVCS data.

    This should give same results as in gepard's
    'fit dvcs dvcs dvcs' or
    smallx-final.nb, section 1-[nlo]-LO.
    """
    fit_gpd = g.model.PWNormGPD()
    m = g.model.MellinBarnesModel(gpds=fit_gpd)
    th = g.theory.BMK(model=m)
    th.m.parameters.update({'ns': 0.15203911208796006,
                            'al0s': 1.1575060246398083,
                            'alps': 0.15,
                            'ms2': 1.,
                            'secs': 0.,
                            'al0g': 1.247316701070471,
                            'alpg': 0.15,
                            'mg2': 0.7})
    # To pre-calculate wce for parallel execution:
    DVCSpoints = g.dset[36].copy()
    for id in range(37, 46):
        DVCSpoints += g.dset[id]
    th.chisq_single(DVCSpoints)
    f = g.fitter.FitterMinuit(DVCSpoints, th)
    f.fix_parameters('ALL')
    f.release_parameters('ms2', 'secs', 'secg')
    print(th.m.parameters_limits)
    print(f.minuit.limits)
    f.minuit.migrad()
    assert th.chisq(f.fitpoints) == approx(95.92, rel=1.e-2)
    assert th.m.parameters['ms2'] == approx(0.47839, rel=1e-2)
    assert th.m.parameters['secs'] == approx(-0.15152, rel=1e-2)
    assert th.m.parameters['secg'] == approx(-0.81216, rel=1e-2)
