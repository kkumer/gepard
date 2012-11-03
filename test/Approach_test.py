
import copy, sys
from nose.tools import *
import numpy as np

import utils, Model, Approach
from results import DMGLO1, DMepsGLO1, DM12  #use some testpars here?

data = utils.loaddata('data/ep2epgamma', approach=Approach.hotfixedBMK)  

m = Model.ModelDR()
m.parameters.update(DMGLO1)
t = Approach.hotfixedBMK(m, optimization = False)

mBM10 = Model.ModelDR()
mBM10.parameters.update(DMepsGLO1)
tBM10 = Approach.BM10(mBM10, optimization = False)

# Optimized model
mopt = Model.ModelDR(optimization = True)
mopt.parameters.update(DMepsGLO1)
mopt.ndparameters = np.array([mopt.parameters[name] for name in mopt.parameter_names]+[0.,0.])
topt = Approach.BM10(mopt, optimization = False)

# Gepard model
mGepard = Model.Gepard(ansatz='EFLEXP')
th = Approach.BMK(mGepard)
mGepard.parameters.update(DM12)


# testing data point for hotfixedBMK
pt0 = copy.deepcopy(data[31][12])  # was data[1][0]
pt0.in1polarization = 1
pt0.in1charge = -1
pt0.FTn = -1
pt0.prepare(Approach.hotfixedBMK)

# testing data point for BM10
pt1 = copy.deepcopy(data[33][-1])
pt1.in2polarization = 1
pt1.phi = 2.*np.pi/5.
pt1.prepare(Approach.BM10)

# testing data points for BM10 and long. TSA and BTSA
ptt = copy.deepcopy(data[52][6])
ptt.prepare(Approach.BM10)
ptb = copy.deepcopy(data[53][3])
ptb.prepare(Approach.BM10)

# testing data point for transversal TSA in BMK
pttrans = copy.deepcopy(ptt)
pttrans.phi = 0.5
pttrans.in2polarizationvector = 'T'
pttrans.varFTn = 1

# testing data point for wBSS
ptw = data[56][0]
# testing data point for wBSD
ptwd = data[55][0]

def test_CFF():
    """Calculate CFF H."""
    assert_almost_equal(m.ImH(pt0), 17.67971592396648)
    assert_almost_equal(m.ReH(pt0), -2.4699741916859592)

test_CFF.one = 1

def test_CFF2():
    """Calculate CFF H."""
    assert_almost_equal(mBM10.ImH(pt1), 1.3213158482535692)
    assert_almost_equal(mBM10.ReH(pt1), -3.8889361918326872)


def test_CFFopt():
    """Calculate CFF H."""
    assert_almost_equal(mopt.ImH(pt1), 1.3213158482535692)
    assert_almost_equal(mopt.ReH(pt1), -3.8889361918326872)

test_CFFopt.optimization = 1

def test_Xunp():
    """Calculate basic cross section Xunp."""
    assert_almost_equal(t.Xunp(pt0, vars={'phi':1.}), 1.8872666478756337)
    ar = t.Xunp(pt0, vars={'phi':np.array([0.3, 1.1])})
    assert isinstance(ar, np.ndarray)
    assert_almost_equal(ar[0], 1.4413231946120821)
    assert_almost_equal(ar[1], 1.9763350136864286)

def test_Xunp2():
    """Any kinematic variable can be in vars."""
    assert_almost_equal(t.Xunp(pt0, 
        vars={'phi':1., 'xB':0.07}),  3.0168274215074025)

def test_Xunp3():
    """ndarray of Q2 could be in vars."""
    ar = t.Xunp(pt0, vars={'phi':1, 'Q2':np.array([2.3, 2.5])})
    assert isinstance(ar, np.ndarray)
    assert_almost_equal(ar[0], 1.9769930014185824)
    assert_almost_equal(ar[1], 2.0929323473733308)

def test_Xunp4():
    """New feature: ndarray of any kinematical variable could be in vars."""
    ar = t.Xunp(pt0, vars={'phi':1, 'xB':np.array([0.07, 0.11])})
    assert isinstance(ar, np.ndarray)
    assert_almost_equal(ar[1], 0.68881435298587579)

test_Xunp4.newfeature = 1

def test_XunpBM10():
    """Calculate unpolarized cross section Xunp in BM10 Approach."""
    assert_almost_equal(tBM10.PreFacSigma(pt1)*tBM10.TBH2unp(pt1), 0.01502937358336803)
    assert_almost_equal(tBM10.PreFacSigma(pt1)*tBM10.TDVCS2unp(pt1), 0.012565093106990456)
    assert_almost_equal(tBM10.PreFacSigma(pt1)*tBM10.TINTunp(pt1), 0.0011255158978939425)
    assert_almost_equal(tBM10.Xunp(pt1), 0.028719982588252427)

def test_XLP():
    """Calculate long. polarized cross section XLP in BM10 Approach."""
    assert_almost_equal(tBM10.PreFacSigma(pt1)*tBM10.TBH2LP(pt1), 0.009495908777414035)
    assert_almost_equal(tBM10.PreFacSigma(pt1)*tBM10.TDVCS2LP(pt1), -0.0032470111398419628)
    assert_almost_equal(tBM10.PreFacSigma(pt1)*tBM10.TINTLP(pt1), 0.0085102074298275109)
    assert_almost_equal(tBM10.XLP(pt1), 0.014759105067399584)

def test_XunpBM10opt():
    """Calculate optimally unpolarized cross section Xunp in BM10 Approach."""

    assert_almost_equal(topt.PreFacSigma(pt1)*topt.TBH2unp(pt1), 0.01502937358336803)
    assert_almost_equal(topt.PreFacSigma(pt1)*topt.TDVCS2unp(pt1), 0.012565093106990456)
    assert_almost_equal(topt.PreFacSigma(pt1)*topt.TINTunp(pt1), 0.0011255158978939425)
    assert_almost_equal(topt.Xunp(pt1), 0.028719982588252427)

test_XunpBM10opt.optimization = 1

def test_XTP():
    """Calculate transv. polarized cross section XTP in BMK Approach."""
    pt1.varphi = 1.
    assert_almost_equal(t.PreFacSigma(pt1)*t.TBH2TP(pt1), 0.0021662047872426475)
    assert_almost_equal(t.PreFacSigma(pt1)*t.TDVCS2TP(pt1), -0.0021305191589529025)
    assert_almost_equal(t.PreFacSigma(pt1)*t.TINTTP(pt1), -0.0098483713204748375)
    assert_almost_equal(t.XTP(pt1), -0.009812685692185092)

def test_BSA():
    """Calculate BSA in BMK Approach."""
    assert_almost_equal(t.BSA(pt0), 0.1845304070958366)
    assert_almost_equal(t.ALUI(pt0), 0.1819945876282851)

def test_TSA():
    """Calculate longitudinal TSA in BM10 Approach."""
    assert_almost_equal(tBM10.TSA(ptt), -0.47969623208934847)

def test_BTSA():
    """Calculate longitudinal BTSA in BM10 Approach."""
    assert_almost_equal(tBM10.BTSA(ptb), 0.25592806446362842)

def test_TTSA():
    """Calculate transversal TSA in BMK Approach and frame."""
    assert_almost_equal(t.TSA(pttrans), 0.080468284490077271)

def test_BSSw():
    """Calculate weighted BSS"""
    assert_almost_equal(t.BSSw(ptw), 0.056334042569159554)

def test_BSDw():
    """Calculate weighted BSD"""
    assert_almost_equal(t.BSDw(ptwd)*1e3, 8.91853141911449)

def test_AUTI():
    """Calculate transversal TSA - INT part - in BMK Approach and frame."""
    assert_almost_equal(t.AUTI(pttrans), 0.075300023640416394)

def test_AUTDVCS():
    """Calculate transversal TSA - DVCS part - in BMK Approach and frame."""
    pttrans.varFTn = -1
    assert_almost_equal(th.AUTDVCS(pttrans)*1e3, -1.5171462298928092)


