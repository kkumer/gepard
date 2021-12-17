"""Testing code for observables."""

import gepard as g
import numpy as np
from pytest import approx, fixture, mark

par_DMGLO1 = {'Nsea': 1.5, 'alS': 1.13,  'alpS': 0.15,  'MS': 0.707107, 'rS': 1.0,
              'bS': 2.00203, 'Nv': 1.35, 'alv': 0.43, 'alpv': 0.85, 'Mv': 1.01097,
              'rv': 0.496383, 'bv': 2.15682, 'C': 6.90484, 'MC': 1.33924, 'tNv': 0.6,
              'tMv': 2.69667, 'trv': 5.97923, 'tbv': 3.25607}

par_DMepsGLO1 = {'C': 6.042, 'MC': 1.5238, 'MS': 0.7071067811865476, 'Mv': 0.8,
                 'Nsea': 1.5, 'Nv': 1.3500000000000001, 'alS': 1.1299999999999999,
                 'alpS': 0.14999999999999999, 'alpv': 0.84999999999999998,
                 'alv': 0.42999999999999999, 'bS': 4.5845, 'bv': 2.398, 'rS': 1.0,
                 'rv': 1.1064, 'tMv': 0.8, 'tNv': 0.6, 'tbv': 1.5, 'trv': 4.8966}


par_DM12 = {"NS": 0.152, "AL0S": 1.1575, "ALPS": 0.1028, "M02S": 0.1814, "DELM2S": 0.,
            "PS": 2., "SECS": 0.5142, "THIS": -0.2106, "SKEWS": 0., "AL0G": 1.2473,
            "ALPG": -0.023, "M02G": 0.2272, "DELM2G": 0., "PG": 2., "SECG": -4.587,
            "THIG": 1.7607, "SKEWG": 0., "KAPS": 1.5265, "EAL0S": 1.1754,
            "EALPS": 0.0239, "EM02S": 0.1808, "EDELM2S": 0., "EPS": 2., "ESECS": 0.5235,
            "ETHIS": -0.2194, "ESKEWS": 0., "EAL0G": 1.8486, "EALPG": 0.05,
            "EM02G": 0.1487, "EDELM2G": 0., "EPG": 2., "ESECG": -4.6366,
            "ETHIG": 1.8638, "ESKEWG": 0.,
            'fix_ALPS': False, 'fix_M02S': False, 'fix_SECS': False, 'fix_THIS': False,
            'fix_KAPS': False, 'fix_ALPG': False, 'fix_M02G': False, 'fix_SECG': False,
            'fix_THIG': False, 'fix_EAL0S': False, 'fix_EALPS': False,
            'fix_EM02S': False, 'fix_ESECS': False, 'fix_ETHIS': False,
            'fix_EAL0G': False, 'fix_EM02G': False, 'fix_ESECG': False}


# testing data point for hotfixedBMK
pt0 = g.data.dset[31][12].copy()
pt0.in1polarization = 1
pt0.in1charge = -1
pt0.FTn = -1
# pt0.prepare(Approach.hotfixedBMK)

# testing data point for BM10
pt1 = g.data.dset[33][-1].copy()
pt1.in2polarization = 1
pt1.phi = 2.*np.pi/5.
# pt1.prepare(Approach.BM10)

# testing data points for BM10 and long. TSA and BTSA
ptt = g.data.dset[52][6].copy()
# ptt.prepare(Approach.BM10)
ptb = g.data.dset[53][3].copy()
# ptb.prepare(Approach.BM10)

# testing data point for transversal TSA in BMK
pttrans = ptt.copy()
pttrans.phi = 0.5
pttrans.in2polarizationvector = 'T'
pttrans.varFTn = 1

# testing data point for wBSS
ptw = g.data.dset[56][0]
# testing data point for wBSD
ptwd = g.data.dset[55][0]


class BMK(g.eff.DipoleEFF, g.cff.DispersionFixedPoleCFF, g.dvcs.hotfixedBMK):
	pass


class BM10(g.eff.DipoleEFF, g.cff.DispersionFixedPoleCFF, g.dvcs.BM10):
	pass


@fixture
def th_BMK():
    th = BMK()
    th.parameters.update(par_DMGLO1)
    return th

@fixture
def th_BM10():
    th = BM10()
    th.parameters.update(par_DMepsGLO1)
    return th

def test_CFF(th_BMK):
    """Calculate CFF H."""
    assert th_BMK.ImH(pt0) == approx(17.67971592396648)
    assert th_BMK.ReH(pt0) == approx(-2.4699741916859592)


def test_CFF2(th_BM10):
    """Calculate CFF H."""
    assert th_BM10.ImH(pt1) == approx(1.3213158482535692)
    assert th_BM10.ReH(pt1) == approx(-3.8889361918326872)


def test_Xunp(th_BMK):
    """Calculate basic cross section Xunp."""
    assert th_BMK.Xunp(pt0, vars={'phi': 1.}) == approx(1.8872666478756337)
    ar = th_BMK.Xunp(pt0, vars={'phi': np.array([0.3, 1.1])})
    assert isinstance(ar, np.ndarray)
    assert ar[0] == approx(1.4413231946120821)
    assert ar[1] == approx(1.9763350136864286)


def test_Xunp2(th_BMK):
    """Any kinematic variable can be in vars."""
    assert th_BMK.Xunp(pt0,
                  vars={'phi': 1., 'xB': 0.07}) == approx(3.0168274215074025)


@mark.skip(reason='Feature not implemented.')
def test_Xunp3(th_BMK):
    """A ndarray of Q2 could be in vars."""
    ar = th_BMK.Xunp(pt0, vars={'phi': 1, 'Q2': np.array([2.3, 2.5])})
    assert isinstance(ar, np.ndarray)
    assert ar[0] == approx(1.9769930014185824)
    assert ar[1] == approx(2.0929323473733308)


@mark.skip(reason='Feature not implemented.')
def test_Xunp4(th_BMK):
    """New feature: ndarray of any kinematical variable could be in vars."""
    ar = th_BMK.Xunp(pt0, vars={'phi': 1, 'xB': np.array([0.07, 0.11])})
    assert isinstance(ar, np.ndarray)
    assert ar[1] == approx(0.68881435298587579)


def test_XunpBM10(th_BM10):
    """Calculate unpolarized cross section Xunp in BM10 Approach."""
    pt1.prepare()
    assert th_BM10.PreFacSigma(pt1)*th_BM10.TBH2unp(pt1) == approx(0.01502937358336803)
    assert th_BM10.PreFacSigma(pt1)*th_BM10.TDVCS2unp(pt1) == approx(0.012565093106990456)
    assert th_BM10.PreFacSigma(pt1)*th_BM10.TINTunp(pt1) == approx(0.0011255158978939425)
    assert th_BM10.Xunp(pt1) == approx(0.028719982588252427)


def test_XLP(th_BM10):
    """Calculate long. polarized cross section XLP in BM10 Approach."""
    pt1.prepare()
    assert th_BM10.PreFacSigma(pt1)*th_BM10.TBH2LP(pt1) == approx(0.009495908777414035)
    assert th_BM10.PreFacSigma(pt1)*th_BM10.TDVCS2LP(pt1) == approx(-0.0032470111398419628)
    assert th_BM10.PreFacSigma(pt1)*th_BM10.TINTLP(pt1) == approx(0.0085102074298275109)
    assert th_BM10.XLP(pt1) == approx(0.014759105067399584)


def test_XTP(th_BMK):
    """Calculate transv. polarized cross section XTP in BMK Approach."""
    pt1.varphi = 1.
    pt1.prepare()
    assert th_BMK.PreFacSigma(pt1)*th_BMK.TBH2TP(pt1) == approx(0.0021662047872426475)
    assert th_BMK.PreFacSigma(pt1)*th_BMK.TDVCS2TP(pt1) == approx(-0.0021305191589529025)
    assert th_BMK.PreFacSigma(pt1)*th_BMK.TINTTP(pt1) == approx(-0.0098483713204748375)
    assert th_BMK.XTP(pt1) == approx(-0.009812685692185092)


def test_BSA(th_BMK):
    """Calculate BSA in BMK Approach."""
    assert th_BMK.BSA(pt0) == approx(0.1845304070958366)
    assert th_BMK.ALUI(pt0) == approx(0.1819945876282851)


def test_TSA(th_BM10):
    """Calculate longitudinal TSA in BM10 Approach."""
    assert th_BM10.TSA(ptt) == approx(-0.47969623208934847)


def test_BTSA(th_BM10):
    """Calculate longitudinal BTSA in BM10 Approach."""
    assert th_BM10.BTSA(ptb) == approx(0.25592806446362842)


def test_TTSA(th_BMK):
    """Calculate transversal TSA in BMK Approach and frame."""
    assert th_BMK.TSA(pttrans) == approx(0.080468284490077271)


def test_BSSw(th_BMK):
    """Calculate weighted BSS."""
    assert th_BMK.BSSw(ptw) == approx(0.056334042569159554)


def test_BSDw(th_BMK):
    """Calculate weighted BSD."""
    assert th_BMK.BSDw(ptwd)*1e3 == approx(8.91853141911449)


def test_AUTI(th_BMK):
    """Calculate transversal TSA - INT part - in BMK Approach and frame."""
    assert th_BMK.AUTI(pttrans) == approx(0.075300023640416394)


@mark.skip(reason='ansatz EFLEXP not yet transferred.')
def test_AUTDVCS():
    """Calculate transversal TSA - DVCS part - in BMK Approach and frame."""
    # Gepard model
    mGepard = g.cff.Gepard(ansatz='EFLEXP')
    mGepard.parameters.update(par_DM12)
    th = g.dvcs.BMK(mGepard)
    pttrans.varFTn = -1
    assert th.AUTDVCS(pttrans)*1e3 == approx(-1.5171462298928092)
