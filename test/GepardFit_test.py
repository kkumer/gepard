
import shutil, copy, math
from nose.tools import *

import numpy as np

import utils, Model, Approach, Fitter

data = utils.loaddata('data/ep2epgamma', approach=Approach.hotfixedBMK)  
data.update(utils.loaddata('data/gammastarp2gammap', approach=Approach.hotfixedBMK))

# Gepard only
mGepard = Model.ComptonGepard()
tGepard = Approach.hotfixedBMK(mGepard)
# removing some limits for compatibility with old Gepard
del mGepard.parameters['limit_M02S']
del mGepard.parameters['limit_M02G']

# Hybrid: Gepard+DR (can reuse above Gepard)
mDRsea = Model.ComptonModelDRsea()
mDRseaopt = Model.ComptonModelDRsea(optimization=True)
m = Model.Hybrid(mGepard, mDRsea)
mopt = Model.Hybrid(mGepard, mDRseaopt)
t = Approach.hotfixedBMK(m)
topt = Approach.hotfixedBMK(mopt)

# Shortcut for seting gepard model parameters
def setpar(i, val):
    mGepard.parameters[mGepard.parameters_index[i]] = val

# DVCSpoints = all data from gepard's dvcs.dat
# model parameters from DIS fit should be fixed
DVCSpoints = data[36] + data[37] + data[38] + data[39] + \
  data[40] + data[41] + data[42] + data[43] + data[44] + \
  data[45]

# testing data set for fits
fitpoints = data[31][12:14] + data[8][1:3] + data[30][2:4]

def test_gepardfitsimple():
    """Test simple fitting to H1 DVCS (one dataset) via gepard"""
    fitpoints = data[36][:4]
    setpar(21, 0.5)   
    setpar(11, 0.15)  
    setpar(12, 1.16)  
    setpar(13, 0.15)  
    setpar(14, 0.4)   
    setpar(16, 2.0)   
    setpar(19, -10.0) 
    setpar(22, 1.25)  
    setpar(23, 0.15)  
    setpar(24, 0.5)   
    setpar(26, 2.0)   
    setpar(29, -32.0) 
    setpar(15, 0.0)   
    setpar(25, 0.0)   
    setpar(17, 0.0)   
    setpar(27, 0.0)   
    setpar(18, 0.0)   
    setpar(28, 0.0)   
    setpar(32, 0.0)
    setpar(42, 0.0)
    mGepard.g.parint.p = 0
    mGepard.g.init()
    tGepard.model.release_parameters('M02S', 'SKEWS')
    f = Fitter.FitterMinuit(fitpoints, tGepard)
    f.fit()
    chisq = tGepard.chisq(fitpoints)[0]
    assert_almost_equal(chisq, 5.1423271052023196, 3)
    tGepard.model.fix_parameters('ALL')

def test_gepardfitDVCSsumso3():
    """Test fitting to HERA DVCS via gepard in sum-so3 model
    
    This should give same results as in smallx-final.nb,
    section 1-[sum]."""
    setpar(11,  0.15203911208796006)
    setpar(12,  1.1575060246398083)
    setpar(13,  0.15)
    setpar(14,  1.)
    setpar(15,  0.)
    setpar(16,  2.)
    setpar(17,  0.)
    setpar(18,  0.)
    setpar(19,  0.)
    setpar(22,  1.247316701070471)
    setpar(23,  0.15)
    setpar(24,  0.7)
    setpar(25,  0.)
    setpar(26,  2.)
    setpar(27,  0.)
    setpar(28,  0.)
    setpar(29,  0.)
    setpar(32,  0.)
    setpar(42,  0.)
    mGepard.g.parint.p = 0
    mGepard.g.init()
    tGepard.model.release_parameters('M02S','SKEWS','SKEWG')
    f = Fitter.FitterMinuit(DVCSpoints, tGepard)
    f.fit()
    chisq = tGepard.chisq(DVCSpoints)[0]
    assert_almost_equal(chisq, 101.094, 1)
    tGepard.model.fix_parameters('ALL')

test_gepardfitDVCSsumso3.long = 1
# It's not a 'new feature' - test passes but the test
# from GepardFitNLO_test.py is more comprehensive
test_gepardfitDVCSsumso3.newfeature = 1

def test_gepardfitDVCSnlso3():
    """Test fitting to HERA DVCS via gepard in nlso3 model
    
    This should give same results as in gepard's
	'fit dvcs dvcs dvcs' or
	smallx-final.nb, section 1-[nlo]-LO."""
    setpar(11,  0.15203911208796006)
    setpar(12,  1.1575060246398083)
    setpar(13,  0.15)
    setpar(14,  1.)
    setpar(15,  0.)
    setpar(16,  2.)
    setpar(17,  0.)
    setpar(18,  0.)
    setpar(19,  0.)
    setpar(22,  1.247316701070471)
    setpar(23,  0.15)
    setpar(24,  0.7)
    setpar(25,  0.)
    setpar(26,  2.)
    setpar(27,  0.)
    setpar(28,  0.)
    setpar(29,  0.)
    setpar(32,  0.)
    setpar(42,  0.)
    mGepard.g.parint.p = 0
    mGepard.g.init()
    tGepard.model.release_parameters('M02S','SECS','SECG')
    f = Fitter.FitterMinuit(DVCSpoints, tGepard)
    f.fit()
    chisq = tGepard.chisq(DVCSpoints)[0]
    assert_almost_equal(chisq/1000, 95.9/1000, 3)
    tGepard.model.fix_parameters('ALL')

test_gepardfitDVCSnlso3.long = 1

def test_gepardfitDVCSthi():
    """Test fitting to HERA DVCS via gepard in 3-PW nlso3 model."""
    setpar(11,  0.15203911208796006)
    setpar(12,  1.1575060246398083)
    setpar(13,  0.15)
    setpar(14,  0.5)
    setpar(15,  0.)
    setpar(16,  2.)
    setpar(17, -0.15)
    setpar(18,  0.)
    setpar(19,  0.)
    setpar(22,  1.247316701070471)
    setpar(23,  0.15)
    setpar(24,  0.7)
    setpar(25,  0.)
    setpar(26,  2.)
    setpar(27, -0.81)
    setpar(28,  0.)
    setpar(29,  0.)
    setpar(32,  0.)
    #setpar(42,  0.1)
    setpar(42,  0.0)
    mGepard.g.parint.p = 0
    mGepard.g.init()
    tGepard.model.release_parameters('M02S','SECG','THIG')
    f = Fitter.FitterMinuit(DVCSpoints, tGepard)
    f.fit()
    chisq = tGepard.chisq(DVCSpoints)[0]
    assert_almost_equal(chisq/1000, 68.888/1000, 3)
    #assert_almost_equal(chisq/1000, 234.597/1000, 3)
    tGepard.model.fix_parameters('ALL')

test_gepardfitDVCSthi.long = 1
    
def test_hybridfitGepard():
    """Test simple hybrid fitting via gepard + DR (DR=0)"""
    fitpoints = data[36][:4]
    setpar(11, 0.15)  
    setpar(12, 1.16)  
    setpar(13, 0.15)  
    setpar(14, 0.4)   
    setpar(16, 2.0)   
    setpar(19, -10.0) 
    setpar(22, 1.25)  
    setpar(23, 0.15)  
    setpar(24, 0.5)   
    setpar(26, 2.0)   
    setpar(29, -32.0) 
    setpar(15, 0.0)   
    setpar(25, 0.0)   
    setpar(17, 0.0)   
    setpar(27, 0.0)   
    setpar(18, 0.0)   
    setpar(28, 0.0)   
    setpar(32, 0.0)
    setpar(42, 0.0)
    t.m.parameters['Nv'] = 0
    t.m.parameters['C'] = 0
    t.m.g.parint.p = 0
    t.m.g.parchr.scheme = np.array([c for c in 'CSBAR'])
    t.m.g.init()
    t.m.release_parameters('M02S', 'SKEWS')
    f = Fitter.FitterMinuit(fitpoints, t)
    f.fit()
    chisq = t.chisq(fitpoints)[0]
    assert_almost_equal(chisq, 5.1423271052023196, 2)
    tGepard.model.fix_parameters('ALL')

    
def test_hybridfitDVCS():
    """Test fitting to HERA DVCS with DR part off.
    
    This should give same results as in smallx-final.nb,
    section 1.-[sum]."""
    setpar(11,  0.15203911208796006)
    setpar(12,  1.1575060246398083)
    setpar(13,  0.15)
    setpar(14,  1.)
    setpar(15,  0.)
    setpar(16,  2.)
    setpar(17,  0.)
    setpar(18,  0.)
    setpar(19,  0.)
    setpar(22,  1.247316701070471)
    setpar(23,  0.15)
    setpar(24,  0.7)
    setpar(25,  0.)
    setpar(26,  2.)
    setpar(27,  0.)
    setpar(28,  0.)
    setpar(29,  0.)
    setpar(32, 0.0)
    setpar(42, 0.0)
    t.m.parameters['Nv'] = 0
    t.m.parameters['C'] = 0
    t.m.g.parint.p = 0
    t.m.g.parchr.scheme = np.array([c for c in 'CSBAR'])
    t.m.g.init()
    t.m.release_parameters('M02S','SKEWS','SKEWG')
    f = Fitter.FitterMinuit(DVCSpoints, t)
    f.fit()
    chisq = t.chisq(DVCSpoints)[0]
    assert_almost_equal(chisq, 101.094, 1)
    tGepard.model.fix_parameters('ALL')

test_hybridfitDVCS.long = 1

def test_gepardfitDIS():
    """Test fitting to H1 DIS via gepard"""
    # DISpoints = all data from gepard's dis.dat
    t.model.release_parameters('NS', 'AL0S', 'AL0G')
    f = Fitter.FitterMinuit(DISpoints, t)
    f.fit()
    chisq = t.chisq(fitpoints)[0]
    assert_almost_equal(chisq, 49.7312, 5)

test_gepardfitDIS.newfeature = 1

def test_hybridfit():
    """Test fitting to large- and small-x data

    with both parts of hybrid model.
    
    """
    setpar(11,  0.15203911208796006)
    setpar(12,  1.1575060246398083)
    setpar(13,  0.15)
    setpar(14,  1.)
    setpar(15,  0.)
    setpar(16,  2.)
    setpar(17,  0.)
    setpar(18,  0.)
    setpar(19,  0.)
    setpar(22,  1.247316701070471)
    setpar(23,  0.15)
    setpar(24,  0.7)
    setpar(25,  0.)
    setpar(26,  2.)
    setpar(27,  0.)
    setpar(28,  0.)
    setpar(29,  0.)
    setpar(32, 0.0)
    setpar(42, 0.0)
    t.m.parameters['Nv'] = 0
    t.m.parameters['C'] = 0
    t.m.g.parint.p = 0
    t.m.g.parchr.scheme = np.array([c for c in 'CSBAR'])
    t.m.g.init()
    t.m.release_parameters('M02S','SKEWG', 'Nv', 'C')
    f = Fitter.FitterMinuit(fitpoints+DVCSpoints[:6], t)
    f.fit()
    chisq = t.chisq(DVCSpoints)[0]
    assert_almost_equal(chisq, 358.52558417612039, 2)
    tGepard.model.fix_parameters('ALL')

test_hybridfit.long = 1

def test_hybridfitopt():
    """Test optimized fitting to large- and small-x data

    with both parts of hybrid model.
    
    """
    setpar(11,  0.15203911208796006)
    setpar(12,  1.1575060246398083)
    setpar(13,  0.15)
    setpar(14,  1.)
    setpar(15,  0.)
    setpar(16,  2.)
    setpar(17,  0.)
    setpar(18,  0.)
    setpar(19,  0.)
    setpar(22,  1.247316701070471)
    setpar(23,  0.15)
    setpar(24,  0.7)
    setpar(25,  0.)
    setpar(26,  2.)
    setpar(27,  0.)
    setpar(28,  0.)
    setpar(29,  0.)
    setpar(32, 0.0)
    setpar(42, 0.0)
    topt.m.parameters['Nv'] = 0
    topt.m.DR.ndparameters[0] = 0
    topt.m.parameters['C'] = 0
    topt.m.DR.ndparameters[12] = 0
    topt.m.g.parint.p = 0
    topt.m.g.parchr.scheme = np.array([c for c in 'CSBAR'])
    topt.m.g.init()
    topt.m.release_parameters('M02S','SKEWG', 'Nv', 'C')
    f = Fitter.FitterMinuit(fitpoints+DVCSpoints[:6], topt)
    f.fit()
    chisq = topt.chisq(DVCSpoints)[0]
    assert_almost_equal(chisq, 358.52558417612039, 2)
    tGepard.model.fix_parameters('ALL')

test_hybridfitopt.long = 1
test_hybridfitopt.newfeature = 1

pt0 = fitpoints[0]
print 'nonoptimized = %f' % m.DR.ReH(pt0)
print 'optimized = %f' % mopt.DR.ReH(pt0)
