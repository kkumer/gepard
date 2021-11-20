"""Testing code for CFFs."""

import gepard as g
import numpy as np
from pytest import approx, mark

par_test = {'ns': 2./3. - 0.4, 'al0s': 1.1, 'alps': 0.25, 'ms2': 1.1**2,
            'ng': 0.4, 'al0g': 1.2, 'alpg': 0.25, 'mg2': 1.2**2}

par_fit = {'ns':  0.152039, 'al0s': 1.15751, 'alps': 0.15, 'ms2': 0.478391,
           'secs': -0.15152, 'this': 0.,  # 'ng': 0.4,  # provided by ns
           'al0g': 1.24732, 'alpg': 0.15, 'mg2': 0.7, 'secg': -0.81217, 'thig': 0.,
           'kaps': 0.7, 'kapg': -0.2}

par_KM15 = {'tMv': 3.992860161655587, 'rS': 1.0, 'alv': 0.43, 'tal': 0.43,
            'Mpi': 3.999999852084612, 'Nv': 1.35, 'rv': 0.918393047884448,
            'Nsea': 0.0, 'alS': 1.13, 'rpi': 2.6463144464701536, 'alpS': 0.15,
            'C': 2.7678681812890016, 'tNv': 0.6, 'bS': 2.0, 'tbv': 0.4000000003259146,
            'bv': 0.40000206775282354, 'Mv': 0.7892078914770488, 'alpv': 0.85,
            'talp': 0.85, 'MC': 1.2040519863464496, 'trv': 0.881085721967267,
            'EPS': 2.0, 'EPG': 2.0, 'EM02G': 0.7, 'EDELM2S': 0.0,
            'SKEWG': 0.0, 'SKEWS': 0.0, 'PS': 2.0, 'EALPG': 0.15, 'M02D': 1.0,
            'EALPS': 0.15, 'M02S': 0.4818827240886959, 'PG': 2.0, 'EDELM2G': 0.0,
            'AL0D': 0.5, 'DELM2S': 0.0, 'DELM2G': 0.0, 'EAL0G': 1.1, 'EAL0S': 1.0,
            'KAPG': 0.0, 'ND': 1.0, 'NG': 0.5, 'SECS': 1.0707825621025808,
            'ESKEWG': 0.0, 'ETHIG': 0.0, 'ESECS': 0.0, 'ETHIS': 0.0, 'ESECG': 0.0,
            'ESKEWS': 0.0, 'DELB': 0.0, 'EM02S': 1.0, 'ALPD': 1.0,
            # This is duplicated, consolidate please
            'MS': 0.4818827240886959,
            # These are renamed already:
            'ns':  0.15203911208796006, 'al0s': 1.1575060246398083, 'alps': 0.15,
            'ms2': 0.4818827240886959, 'secs': 1.0707825621025808,
            'this': -0.36618269477432946,
            'al0g': 1.247316701070471, 'alpg': 0.15, 'mg2': 0.7,
            'secg': -2.990809378821039, 'thig': 0.9052207712570559,
            'kaps': 0.0, 'kapg': 0.0}


pt_test = g.data.DataPoint({'W': 82., 'Q2': 1., 't': 0.})
pt_test.xi = pt_test.Q2 / (2.0 * pt_test.W * pt_test.W + pt_test.Q2)

pt0_fit = g.data.DataPoint({'xi': 0.01, 'Q2': 4., 't': -0.2})  # noevol
pt_fit = g.data.DataPoint({'xi': 0.01, 'Q2': 8., 't': -0.2})  # evol

pt_bp = g.data.DataPoint({'xi': 1.e-5, 'Q2': 2.5, 't': -0.25})  # noevol
pt_evol = g.data.DataPoint({'xi': 1.e-5, 'Q2': 25., 't': -0.25})  # evol
pt_evolNS = g.data.DataPoint({'xi': 0.01, 'Q2': 10., 't': -0.25})  # evol


class gpd_model(g.model.ConformalSpaceGPD):
    """GPD model from the paper - singlet and nonsinglet part."""

    def __init__(self, type, **kwargs) -> None:
        """Defaults from old Fortran's GEPARD.INI."""
        kwargs.setdefault('p', 1)
        kwargs.setdefault('scheme', 'csbar')
        kwargs.setdefault('nf', 4)
        kwargs.setdefault('q02', 2.5)
        kwargs.setdefault('asp', np.array([0.05, 0.05, 0.05]))
        kwargs.setdefault('phi', 1.9)
        self.type = type
        super().__init__(**kwargs)
        R = 0.5  # ratio sbar/ubar
        self.frot = np.array([[1, 0, 1, 1],
                              [0, 1, 0, 0],
                              [-R/(2+R), 0, 1, -1]])

    def gpd_H(self, eta: float, t: float) -> np.ndarray:
        """GPD H from hep-ph/0703179."""
        return g.gpdj.ansatz07_fixed(self.jpoints, t, self.type).transpose()


def test_wc_LO():
    """Test LO DVCS Wilson coef."""
    test_gpd = g.model.Test(p=0)
    m_test = g.model.MellinBarnesModel(gpds=test_gpd)
    m_test.parameters.update(par_test)
    assert g.evolc.calc_wc(m_test, m_test.jpoints, 'DVCS')[0, 0, :2] == approx(
            np.array([1.7798226558761627+0.00017759121554287j, 0+0j]))


def test_wc_NLO():
    """Test NLO DVCS Wilson coef."""
    test_gpd = g.model.Test(p=1)
    m_test = g.model.MellinBarnesModel(gpds=test_gpd)
    m_test.parameters.update(par_test)
    assert g.evolc.calc_wc(m_test, m_test.jpoints, 'DVCS')[0, 1, :2] == approx(
            np.array([-0.88174829594212023+0.00093822077679447j,
                      -5.9050162592671382-0.00044618938685837j]))


def test_wce_LO():
    """Test LO DVCS evolved Wilson coef."""
    test_gpd = g.model.Test(p=0)
    m_test = g.model.MellinBarnesModel(gpds=test_gpd)
    m_test.parameters.update(par_test)
    assert g.evolc.calc_wce(m_test, 3.0, 'DVCS')[0, 0, :2] == approx(
            np.array([1.7328455630029231+0.00009701899018317j,
                      0.21666921098074668-0.0000851087087619j]), rel=1.e-12)


def test_wce_NLO():
    """Test NLO DVCS evolved Wilson coef."""
    test_gpd = g.model.Test(p=1)
    m_test = g.model.MellinBarnesModel(gpds=test_gpd)
    m_test.parameters.update(par_test)
    assert g.evolc.calc_wce(m_test, 3.0, 'DVCS')[0, 0, :2] == approx(
            np.array([1.6127545996599677+0.00014769567470216j,
                      -0.09044960485326564-0.00003265190306802j]), rel=1.e-10)


def test_cff_H_noevol():
    """Test testing (ReH, ImH) (no evolution)."""
    test_gpd = g.model.Test()
    m_test = g.model.MellinBarnesModel(gpds=test_gpd)
    m_test.parameters.update(par_test)
    assert m_test.cff(pt_test)[:2] == approx(
            [9839.566, 61614.9])


def test_cff_radLO():
    """Singlet LO CFF H (no evol)."""
    gpds = gpd_model('hard', p=0, scheme='csbar')
    m = g.model.MellinBarnesModel(gpds=gpds)
    qs = 5/18  # for nf=4
    m.dvcs_charges = (qs, qs, 0)  # select only singlet part of CFF
    assert m.cff(pt_bp)[:2] == approx(
            [39544.823112887607, 402367.23596533033])


def test_cff_radNLO():
    """Singlet NLO CFF H (no evol)."""
    gpds = gpd_model('hard', p=1, scheme='csbar')
    m = g.model.MellinBarnesModel(gpds=gpds)
    qs = 5/18  # for nf=4
    m.dvcs_charges = (qs, qs, 0)  # select only singlet part of CFF
    assert m.cff(pt_bp)[:2] == approx(
            [5747.0424614455933, 201256.45352582674])


def test_cff_radLO_evol_NS():
    """Non-singlet LO CFF H (evol)."""
    gpds = gpd_model('hardNS', p=0, scheme='csbar')
    m = g.model.MellinBarnesModel(gpds=gpds)
    qns = 1/6  # for nf=4
    m.dvcs_charges = (0, 0, qns)  # select only NS part of CFF
    assert m.cff(pt_evolNS)[:2] == approx(
            [-9.554282705911017, -25.83697714227015])


def test_cff_radNLO_CSBARevol_NS():
    """Non-singlet NLO CFF H (CSBAR evol)."""
    gpds = gpd_model('hardNS', p=1, scheme='csbar')
    m = g.model.MellinBarnesModel(gpds=gpds)
    qns = 1/6  # for nf=4
    m.dvcs_charges = (0, 0, qns)  # select only NS part of CFF
    assert m.cff(pt_evolNS)[:2] == approx(
            [-8.666225662660972, -24.309344967448652])


@mark.slow
def test_cff_radNLO_MSBARevol_NS():
    """Non-singlet NLO CFF H (CSBAR evol)."""
    gpds = gpd_model('hardNS', p=1, scheme='msbar')
    m = g.model.MellinBarnesModel(gpds=gpds)
    qns = 1/6  # for nf=4
    m.dvcs_charges = (0, 0, qns)  # select only NS part of CFF
    assert m.cff(pt_evolNS)[:2] == approx(
            [-8.26592013777189, -22.86536625486948])


def test_cff_rad_LOevol():
    """Singlet LO CFF H evolved."""
    gpds = gpd_model('hard', p=0, scheme='csbar')
    m = g.model.MellinBarnesModel(gpds=gpds)
    qs = 5/18  # for nf=4
    m.dvcs_charges = (qs, qs, 0)  # select only singlet part of CFF
    assert m.cff(pt_evol)[:2] == approx(
            [251460.03959908773, 1015357.1865059549])


def test_cff_radNLO_CSBAR_evol():
    """Singlet NLO CSBAR CFF H evolved."""
    gpds = gpd_model('hard', p=1, scheme='csbar')
    m = g.model.MellinBarnesModel(gpds=gpds)
    qs = 5/18  # for nf=4
    m.dvcs_charges = (qs, qs, 0)  # select only singlet part of CFF
    assert m.cff(pt_evol)[:2] == approx(
            [162468.6375694127, 729219.8768488609])


@mark.slow
def test_cff_radNLO_MSBAR_evol():
    """Singlet NLO MSBAR CFF H evolved."""
    gpds = gpd_model('hard', p=1, scheme='msbar')
    m = g.model.MellinBarnesModel(gpds=gpds)
    qs = 5/18  # for nf=4
    m.dvcs_charges = (qs, qs, 0)  # select only singlet part of CFF
    # Result of wrong ND-evolution fortran-gepard code
    # assert m.cff(pt_evol)[:2] == approx(
    #         [142867.21556625995, 653095.26655367797/1e5])
    # Result after correcting ND-evolution fortran-gepard code
    assert m.cff(pt_evol)[:2] == approx(
            [156576.80414343436, 686720.2142542489], rel=1.e-5)


def test_cff_H_nlso3():
    """Test nl-so3 (ReH, ImH) (LO evolved to multiple Q2)."""
    fit_gpd = g.model.PWNormGPD()
    m_fit = g.model.MellinBarnesModel(gpds=fit_gpd)
    m_fit.parameters.update(par_fit)
    # Q2 can be changed during calls:
    assert m_fit.cff(pt0_fit)[:2] == approx(
            [18.9540, 60.1622], abs=0.001)
    assert m_fit.cff(pt_fit)[:2] == approx(
            (26.9984, 66.5255), abs=0.001)


def test_cff_H_nlso3_separate():
    """Test nl-so3 ReH, ImH (LO evolved to multiple Q2)."""
    fit_gpd = g.model.PWNormGPD()
    m_fit = g.model.MellinBarnesModel(gpds=fit_gpd)
    m_fit.parameters.update(par_fit)
    # Q2 can be changed during calls:
    assert m_fit.ReH(pt0_fit) == approx(
            18.9540, abs=0.001)
    assert m_fit.ImH(pt0_fit) == approx(
            60.1622, abs=0.001)
    assert m_fit.ReH(pt_fit) == approx(
            26.9984, abs=0.001)
    assert m_fit.ImH(pt_fit) == approx(
            66.5255, abs=0.001)


def test_cff_E_nlso3():
    """Testing nl-so3 (ReE, ImE) (LO evolved to multiple Q2)."""
    fit_gpd = g.model.PWNormGPD()
    m_fit = g.model.MellinBarnesModel(gpds=fit_gpd)
    m_fit.parameters.update(par_fit)
    assert m_fit.cff(pt0_fit)[2:4] == approx(
            [13.2678, 42.11355])
    assert m_fit.cff(pt_fit)[2:4] == approx(
            [13.5297, 43.6024], rel=1e-4)


def test_cff_E_nlso3_separate():
    """Testing nl-so3 ReE, ImE (LO evolved to multiple Q2)."""
    fit_gpd = g.model.PWNormGPD()
    m_fit = g.model.MellinBarnesModel(gpds=fit_gpd)
    m_fit.parameters.update(par_fit)
    assert m_fit.ReE(pt0_fit) == approx(
            13.2678)
    assert m_fit.ImE(pt0_fit) == approx(
            42.11355)
    assert m_fit.ReE(pt_fit) == approx(
            13.5297, rel=1e-4)
    assert m_fit.ImE(pt_fit) == approx(
            43.6024, rel=1e-4)


# def test_cff_H_KM15():
    # """Test of MB part of KM15 - temp.."""
    # fit_gpd = g.model.PWNormGPD()
    # m_fit = g.model.MellinBarnesModel(gpds=fit_gpd)
    # m_fit.parameters.update(par_KM15)
    # # Q2 can be changed during calls:
    # assert m_fit.ImH(pt0_fit) == approx(66.90122688611703)
    # assert m_fit.cff(pt_fit)[1] == approx(
            # 66.90122688611703)


def test_KM15_cffs():
    """Test CFFs of KM15 model fit."""
    fit_gpd = g.model.PWNormGPD()
    mMB = g.model.MellinBarnesModel(gpds=fit_gpd)
    mDR = g.model.ComptonModelDRPP()
    m = g.model.ComptonHybrid(mMB, mDR)
    m.parameters.update(par_KM15)
    assert m.ImH(pt0_fit) == approx(
            70.3619128283078)
