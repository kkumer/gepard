"""Testing code for CFFs."""

import gepard as g
import numpy as np
from gepard.fits import par_KM15
from pytest import approx, fixture, mark, raises

par_test = {'ns': 2./3. - 0.4, 'al0s': 1.1, 'alps': 0.25, 'ms2': 1.1**2,
            'ng': 0.4, 'al0g': 1.2, 'alpg': 0.25, 'mg2': 1.2**2}

par_fit = {'ns':  0.152039, 'al0s': 1.15751, 'alps': 0.15, 'ms2': 0.478391,
           'secs': -0.15152, 'this': 0.,  # 'ng': 0.4,  # provided by ns
           'al0g': 1.24732, 'alpg': 0.15, 'mg2': 0.7, 'secg': -0.81217, 'thig': 0.,
           'kaps': 0.7, 'kapg': -0.2}

pt_test = g.data.DataPoint({'W': 82., 'Q2': 1., 't': 0.})
pt_test.xi = pt_test.Q2 / (2.0 * pt_test.W * pt_test.W + pt_test.Q2)

pt0_fit = g.data.DataPoint({'xi': 0.01, 'Q2': 4., 't': -0.2})  # noevol
pt_fit = g.data.DataPoint({'xi': 0.01, 'Q2': 8., 't': -0.2})  # evol

# DataPoint should be in gepard namespace as well:
pt_bp = g.DataPoint({'xi': 1.e-5, 'Q2': 2.5, 't': -0.25})  # noevol
pt_evol = g.DataPoint({'xi': 1.e-5, 'Q2': 25., 't': -0.25})  # evol
pt_evolNS = g.DataPoint({'xi': 0.01, 'Q2': 10., 't': -0.25})  # evol


class gpd_model(g.gpd.ConformalSpaceGPD):
    """GPD model from the paper - singlet and nonsinglet part."""

    def __init__(self, type, **kwargs) -> None:
        """Defaults from old Fortran's GEPARD.INI."""
        kwargs.setdefault('p', 1)
        kwargs.setdefault('scheme', 'csbar')
        kwargs.setdefault('nf', 4)
        kwargs.setdefault('Q02', 2.5)
        kwargs.setdefault('asp', np.array([0.05, 0.05, 0.05]))
        kwargs.setdefault('phi', 1.9)
        self.type = type
        super().__init__(**kwargs)
        R = 0.5  # ratio sbar/ubar
        self.frot = np.array([[1, 0, 1, 1],
                              [0, 1, 0, 0],
                              [-R/(2+R), 0, 1, -1]])

    def H(self, eta: float, t: float) -> np.ndarray:
        """GPD H from hep-ph/0703179."""
        return g.gpd.ansatz07_fixed(self.jpoints, t, self.type).transpose()


class CFFTest1(g.gpd.TestGPD, g.cff.MellinBarnesCFF):
	pass

@fixture
def th1_lo():
    th = CFFTest1(p=0)
    th.parameters.update(par_test)
    return th

@fixture
def th1_nlo():
    th = CFFTest1(p=1)
    th.parameters.update(par_test)
    return th


class CFFTest2(gpd_model, g.cff.MellinBarnesCFF):
	pass


class CFFTest3(g.gpd.PWNormGPD, g.cff.MellinBarnesCFF):

    def pw_strengths_E(self):
        """For this model, PW strenghts are same as for H."""
        return self.pw_strengths()

    def E(self, eta: float, t: float) -> np.ndarray:
        """Return (npts, 4) array E_j^a for all j-points and 4 flavors."""
        kappa = np.array([self.parameters['kaps'], self.parameters['kapg'], 0, 0])
        return kappa * g.gpd.singlet_ng_constrained(self.jpoints, t, self.parameters,
                self.residualt).transpose()


class CFFTest4(g.gpd.PWNormGPD, g.cff.HybridFixedPoleCFF):
	pass


def test_wc_LO(th1_lo):
    """Test LO DVCS Wilson coef."""
    assert g.wilson.calc_wc(th1_lo, th1_lo.jpoints, 'DVCS')[0, 0, :2] == approx(
            np.array([1.7798226558761627+0.00017759121554287j, 0+0j]))


def test_wc_NLO(th1_nlo):
    """Test NLO DVCS Wilson coef."""
    assert g.wilson.calc_wc(th1_nlo, th1_nlo.jpoints, 'DVCS')[0, 1, :2] == approx(
            np.array([-0.88174829594212023+0.00093822077679447j,
                      -5.9050162592671382-0.00044618938685837j]))


def test_wce_LO(th1_lo):
    """Test LO DVCS evolved Wilson coef."""
    assert g.wilson.calc_wce(th1_lo, 3.0, 'DVCS')[0, 0, :2] == approx(
            np.array([1.7328455630029231+0.00009701899018317j,
                      0.21666921098074668-0.0000851087087619j]), rel=1.e-12)


def test_wce_NLO(th1_nlo):
    """Test NLO DVCS evolved Wilson coef."""
    assert g.wilson.calc_wce(th1_nlo, 3.0, 'DVCS')[0, 0, :2] == approx(
            np.array([1.6127545996599677+0.00014769567470216j,
                      -0.09044960485326564-0.00003265190306802j]), rel=1.e-10)


def test_cff_H_noevol(th1_lo):
    """Test testing (ReH, ImH) (no evolution)."""
    assert th1_lo.cff(pt_test)[:2] == approx(
            [9839.566, 61614.9])


def test_cff_radLO():
    """Singlet LO CFF H (no evol)."""
    m = CFFTest2(type='hard', p=0, scheme='csbar')
    qs = 5/18  # for nf=4
    m.dvcs_charges = (qs, qs, 0)  # select only singlet part of CFF
    assert m.cff(pt_bp)[:2] == approx(
            [39544.823112887607, 402367.23596533033])


def test_cff_radNLO():
    """Singlet NLO CFF H (no evol)."""
    m = CFFTest2(type='hard', p=1, scheme='csbar')
    qs = 5/18  # for nf=4
    m.dvcs_charges = (qs, qs, 0)  # select only singlet part of CFF
    assert m.cff(pt_bp)[:2] == approx(
            [5747.0424614455933, 201256.45352582674])


def test_cff_proc_class_except():
    """Trigger process_class exception."""
    m = CFFTest2(type='hard', p=1, scheme='csbar')
    qs = 5/18  # for nf=4
    m.dvcs_charges = (qs, qs, 0)  # select only singlet part of CFF
    with raises(Exception, match='process_class SIDIS is not DIS, DVCS or DVMP!'):
        g.wilson.calc_wc(m, m.jpoints, 'SIDIS')
    with raises(Exception, match='process_class SIDIS is neither DVCS nor DIS!'):
        g.c1dvcs.shift1(m, m.jpoints, 'SIDIS')


def test_cff_scheme_except():
    """Trigger scheme exception."""
    m = CFFTest2(type='hard', p=1, scheme='ms')
    qs = 5/18  # for nf=4
    m.dvcs_charges = (qs, qs, 0)  # select only singlet part of CFF
    with raises(Exception, match='Scheme ms is neither msbar nor csbar!'):
        m.cff(pt_bp)


def test_cff_radLO_evol_NS():
    """Non-singlet LO CFF H (evol)."""
    m = CFFTest2(type='hardNS', p=0, scheme='csbar')
    qns = 1/6  # for nf=4
    m.dvcs_charges = (0, 0, qns)  # select only NS part of CFF
    assert m.cff(pt_evolNS)[:2] == approx(
            [-9.554282705911017, -25.83697714227015])


def test_cff_radNLO_CSBARevol_NS():
    """Non-singlet NLO CFF H (CSBAR evol)."""
    m = CFFTest2(type='hardNS', p=1, scheme='csbar')
    qns = 1/6  # for nf=4
    m.dvcs_charges = (0, 0, qns)  # select only NS part of CFF
    assert m.cff(pt_evolNS)[:2] == approx(
            [-8.666225662660972, -24.309344967448652])


@mark.slow
def test_cff_radNLO_MSBARevol_NS():
    """Non-singlet NLO CFF H (CSBAR evol)."""
    m = CFFTest2(type='hardNS', p=1, scheme='msbar')
    qns = 1/6  # for nf=4
    m.dvcs_charges = (0, 0, qns)  # select only NS part of CFF
    assert m.cff(pt_evolNS)[:2] == approx(
            [-8.26592013777189, -22.86536625486948])


def test_cff_rad_LOevol():
    """Singlet LO CFF H evolved."""
    m = CFFTest2(type='hard', p=0, scheme='csbar')
    qs = 5/18  # for nf=4
    m.dvcs_charges = (qs, qs, 0)  # select only singlet part of CFF
    assert m.cff(pt_evol)[:2] == approx(
            [251460.03959908773, 1015357.1865059549])


def test_cff_radNLO_CSBAR_evol():
    """Singlet NLO CSBAR CFF H evolved."""
    m = CFFTest2(type='hard', p=1, scheme='csbar')
    qs = 5/18  # for nf=4
    m.dvcs_charges = (qs, qs, 0)  # select only singlet part of CFF
    assert m.cff(pt_evol)[:2] == approx(
            [162468.6375694127, 729219.8768488609])


@mark.slow
def test_cff_radNLO_MSBAR_evol():
    """Singlet NLO MSBAR CFF H evolved."""
    m = CFFTest2(type='hard', p=1, scheme='msbar')
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
    m_fit = CFFTest3()
    m_fit.parameters.update(par_fit)
    # Q2 can be changed during calls:
    assert m_fit.cff(pt0_fit)[:2] == approx(
            [18.9540, 60.1622], abs=0.001)
    assert m_fit.cff(pt_fit)[:2] == approx(
            (26.9984, 66.5255), abs=0.001)


def test_cff_H_nlso3_separate():
    """Test nl-so3 ReH, ImH (LO evolved to multiple Q2)."""
    m_fit = CFFTest3()
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
    m_fit = CFFTest3()
    m_fit.parameters.update(par_fit)
    assert m_fit.cff(pt0_fit)[2:4] == approx(
            [13.2678, 42.11355])
    assert m_fit.cff(pt_fit)[2:4] == approx(
            [13.5297, 43.6024], rel=1e-4)


def test_cff_E_nlso3_separate():
    """Testing nl-so3 ReE, ImE (LO evolved to multiple Q2)."""
    m_fit = CFFTest3()
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
    # m_fit = CFFTest3()
    # m_fit.parameters.update(par_KM15)
    # # Q2 can be changed during calls:
    # assert m_fit.ImH(pt0_fit) == approx(66.90122688611703)
    # assert m_fit.cff(pt_fit)[1] == approx(
    #       # 66.90122688611703)


def test_KM15_cffs():
    """Test CFFs of KM15 model fit."""
    m = CFFTest4()
    m.parameters.update(par_KM15)
    assert m.ImH(pt0_fit) == approx(
            70.3619128283078)
    assert m.ReH(pt0_fit) == approx(
            7.764977856195072)


def test_print_cffs(capfd):
    """Test printing CFFs at a point."""
    m = CFFTest4()
    m.parameters.update(par_KM15)
    m.print_CFFs(pt0_fit)
    out, err = capfd.readouterr()
    assert out[:25] == " ImH = 70.36\n ReH =  7.76"


def test_print_cffs_mma(capfd):
    """Test printing CFFs at a point (mma style)."""
    m = CFFTest4()
    m.parameters.update(par_KM15)
    m.print_CFFs(pt0_fit, format='mma')
    out, err = capfd.readouterr()
    assert out[:51] == "{ImH -> 70.361913, ReH -> 7.764978, ImE -> 0.000000"
