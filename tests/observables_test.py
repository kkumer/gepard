"""Testing code for observables."""

import gepard as g
import numpy as np
from pytest import approx, fixture, mark

par_DMGLO1 = {'Nsea': 1.5, 'alS': 1.13,  'alpS': 0.15,  'mS2': 0.5, 'rS': 1.0,
              'bS': 2.00203, 'Nv': 1.35, 'alv': 0.43, 'alpv': 0.85, 'mv2': 1.0220603409,
              'rv': 0.496383, 'bv': 2.15682, 'C': 6.90484, 'mC2': 1.7935637776, 'tNv': 0.6,
              'tmv2': 7.2720290889, 'trv': 5.97923, 'tbv': 3.25607}

par_DMepsGLO1 = {'C': 6.042, 'mC2': 2.32196644, 'mS2': 0.5, 'mv2': 0.64,
                 'Nsea': 1.5, 'Nv': 1.35, 'alS': 1.13,
                 'alpS': 0.15, 'alpv': 0.85,
                 'alv': 0.43, 'bS': 4.5845, 'bv': 2.398, 'rS': 1.0,
                 'rv': 1.1064, 'tmv2': 0.64, 'tNv': 0.6, 'tbv': 1.5, 'trv': 4.8966}

par_DM12 = {"ns": 0.152, "al0s": 1.1575, "alps": 0.1028, "ms2": 0.1814, "delms2": 0.,
            "pows": 2., "secs": 0.5142, "this": -0.2106, # "SKEWS": 0.,
            "al0g": 1.2473, "alpg": -0.023, "mg2": 0.2272, "delmg2": 0.,
            "powg": 2., "secg": -4.587, "thig": 1.7607, # "SKEWG": 0.,
            "kaps": 1.5265, "kapg": 0,  # kapg not present in orig
            "Ens": 0.152,  # for fitting should be pegged to ns
            "Eal0s": 1.1754, "Ealps": 0.0239, "Ems2": 0.1808, "Edelms2": 0.,
            "Epows": 2., "Esecs": 0.5235, "Ethis": -0.2194, # "ESKEWS": 0.,
            "Eal0g": 1.8486, "Ealpg": 0.05, "Emg2": 0.1487, "Edelmg2": 0.,
            "Epowg": 2., "Esecg": -4.6366, "Ethig": 1.8638} #, "ESKEWG": 0.}


# testing data point for hotfixedBMK
pt0 = g.data.dset[31][12].copy()
pt0.in1polarization = 1
pt0.in1charge = -1
pt0.FTn = -1
# pt0.prepare(Approach.hotfixedBMK)

# testing data point for BM10
pt1 = g.data.dset[33][-1].copy()
pt1.phi = 2.*np.pi/5.

# unpolarized testing data point for BM10
pt2 = g.data.dset[33][-1].copy()
pt2.phi = 2.*np.pi/5.

# testing data points for BM10 and long. TSA and BTSA
ptt = g.data.dset[52][6].copy()
# ptt.prepare(Approach.BM10)
ptb = g.data.dset[53][3].copy()
# ptb.prepare(Approach.BM10)

# testing data point for transversal TSA in BMK
pttrans = ptt.copy()
pttrans.phi = 0.5
pttrans.in2polarizationvector = 'T'
pttrans.in2polarization = 1
pttrans.varFTn = 1

# testing data point for wXUU
ptw = g.data.dset[56][0]
# testing data point for wXLU
ptwd = g.data.dset[55][0]


class BMK(g.eff.DipoleEFF, g.cff.DispersionFixedPoleCFF, g.dvcs.hotfixedBMK):
	pass


class BM10(g.eff.DipoleEFF, g.cff.DispersionFixedPoleCFF, g.dvcs.BM10):
	pass

class DM12(g.eff.KellyEFF, g.gpd.PWNormGPD, g.cff.MellinBarnesCFF, g.dvcs.BMK):

    def E(self, eta: float, t: float) -> np.ndarray:
        """Return (npts, 4) array E_j^a for all j-points and 4 flavors."""
        # Implement BS+BG=0 sum rule that fixes 'kapg'
        self.parameters['kapg'] = - self.parameters['kaps'] * self.parameters['ns'] / (
                0.6 - self.parameters['ns'])
        kappa = np.array([self.parameters['kaps'], self.parameters['kapg'], 0, 0])
        return kappa * g.gpd.singlet_ng_constrained_E(self.jpoints, t,
                            self.parameters, self.residualt).transpose()

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

@fixture
def th_DM12():
    th = DM12(residualt='exp')
    th.parameters.update(par_DM12)
    return th

def test_CFF(th_BMK):
    """Calculate CFF H."""
    assert th_BMK.ImH(pt0) == approx(17.67971592396648)
    assert th_BMK.ReH(pt0) == approx(-2.4699741916859592)


def test_CFF2(th_BM10):
    """Calculate CFF H."""
    assert th_BM10.ImH(pt1) == approx(1.3213158482535692)
    assert th_BM10.ReH(pt1) == approx(-3.8889361918326872)


def test_XS(th_BMK):
    """Calculate basic cross section XS."""
    assert th_BMK.XS(pt0, vars={'phi': 1.}) == approx(1.8872666478756337)
    ar = th_BMK.XS(pt0, vars={'phi': np.array([0.3, 1.1])})
    assert isinstance(ar, np.ndarray)
    assert ar[0] == approx(1.4413231946120821)
    assert ar[1] == approx(1.9763350136864286)


def test_XS_2(th_BMK):
    """Any kinematic variable can be in vars."""
    assert th_BMK.XS(pt0,
                  vars={'phi': 1., 'xB': 0.07}) == approx(3.0168274215074025)


@mark.skip(reason='Feature not implemented.')
def test_XS_3(th_BMK):
    """A ndarray of Q2 could be in vars."""
    ar = th_BMK.XS(pt0, vars={'phi': 1, 'Q2': np.array([2.3, 2.5])})
    assert isinstance(ar, np.ndarray)
    assert ar[0] == approx(1.9769930014185824)
    assert ar[1] == approx(2.0929323473733308)


@mark.skip(reason='Feature not implemented.')
def test_XS_4(th_BMK):
    """New feature: ndarray of any kinematical variable could be in vars."""
    ar = th_BMK.XS(pt0, vars={'phi': 1, 'xB': np.array([0.07, 0.11])})
    assert isinstance(ar, np.ndarray)
    assert ar[1] == approx(0.68881435298587579)


def test_XS_unp_BM10(th_BM10):
    """Calculate unpolarized cross section in BM10 Approach."""
    pt2.prepare()
    assert th_BM10.PreFacSigma(pt2)*th_BM10.TBH2unp(pt2) == approx(0.01502937358336803)
    assert th_BM10.PreFacSigma(pt2)*th_BM10.TDVCS2unp(pt2) == approx(0.012565093106990456)
    assert th_BM10.PreFacSigma(pt2)*th_BM10.TINTunp(pt2) == approx(0.0011255158978939425)
    assert th_BM10.XS(pt2) == approx(0.028719982588252427)


def test_XUL(th_BM10):
    """Calculate long. polarized target cross section XUL in BM10 Approach."""
    pt = pt1.copy()
    pt.in2polarizationvector = 'L'
    pt.in2polarization = 1
    pt.prepare()
    assert th_BM10.PreFacSigma(pt)*th_BM10.TBH2LP(pt) == approx(0.009495908777414035)
    assert th_BM10.PreFacSigma(pt)*th_BM10.TDVCS2LP(pt) == approx(-0.0032470111398419628)
    assert th_BM10.PreFacSigma(pt)*th_BM10.TINTLP(pt) == approx(0.0085102074298275109)
    assert th_BM10.XUL(pt) == approx(0.014759105067399584)


def test_XUT(th_BMK):
    """Calculate transv. polarized target cross section XUT in BMK Approach."""
    pt = pt1.copy()
    pt.varphi = 1.
    pt.in2polarizationvector = 'T'
    pt.in2polarization = 1
    pt.prepare()
    assert th_BMK.PreFacSigma(pt)*th_BMK.TBH2TP(pt) == approx(0.0021662047872426475)
    assert th_BMK.PreFacSigma(pt)*th_BMK.TDVCS2TP(pt) == approx(-0.0021305191589529025)
    assert th_BMK.PreFacSigma(pt)*th_BMK.TINTTP(pt) == approx(-0.0098483713204748375)
    assert th_BMK.XUT(pt) == approx(-0.009812685692185092)

def test_ALU_approx(th_BMK):
    """Test approximate formula for ALU in BMK Approach."""
    assert th_BMK._ALUapprox(pt0) == approx(0.1845304070958366)

def test_ALU_exact(th_BMK):
    """Test ALU in BMK Approach."""
    assert th_BMK._ALUexact(pt0) == approx(0.1847045934136741)

def test_ALUI(th_BMK):
    """Test ALUI in BMK Approach."""
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


def test_XUUw(th_BMK):
    """Calculate weighted XUU."""
    assert th_BMK.XUUw(ptw) == approx(0.056334042569159554)


def test_XLUw(th_BMK):
    """Calculate weighted XLU."""
    assert th_BMK.XLUw(ptwd)*1e3 == approx(8.91853141911449)


def test_AUTI(th_BMK):
    """Calculate transversal TSA - INT part - in BMK Approach and frame."""
    assert th_BMK.AUTI(pttrans) == approx(0.075300023640416394)


def test_AUTDVCS(th_DM12):
    """Calculate transversal TSA - DVCS part - in BMK Approach and frame."""
    # Gepard model
    pttrans.varFTn = -1
    assert th_DM12.AUTDVCS(pttrans) == approx(-0.0015171462298928092)


def test_repr_DataPoint():
    """Representation of DataPoint for printing."""
    assert repr(g.dset[32][0]) == 'DataPoint: AC = -0.027'


def test_DataSet_slice():
    """Slicing the DataSet."""
    assert repr(g.dset[32][1:4]) == 'DataSet with 3 points'
