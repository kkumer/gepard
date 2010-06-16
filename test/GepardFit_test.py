
import shutil, copy, math
from nose.tools import *
import numpy as np

import gepard as g

import utils, Model, Approach, Fitter

data = utils.loaddata('data/ep2epgamma')  
data.update(utils.loaddata('data/gammastarp2gammap'))

shutil.copy2('test/GEPARD.INI.FIT', 'GEPARD.INI')

# Gepard only
mGepard = Model.ComptonGepard()
tGepard = Approach.hotfixedBMK(mGepard)

# Hybrid: Gepard+DR (can reuse above Gepard)
mDRsea = Model.ComptonModelDRsea()
m = Model.Hybrid(mGepard, mDRsea)
t = Approach.hotfixedBMK(m)

# Shortcut for seting gepard model parameters
def setpar(i, val):
    mGepard.parameters[mGepard.parameters_index[i]] = val



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
    mGepard.g.parint.p = 0
    mGepard.g.init()
    tGepard.model.release_parameters('M02S', 'SKEWS')
    f = Fitter.FitterMinuit(fitpoints, tGepard)
    f.fit()
    chisq = tGepard.chisq(fitpoints)[0]
    assert_almost_equal(chisq, 5.1423271052023196, 3)

def test_gepardfitDVCS():
    """Test fitting to HERA DVCS via gepard
    
    This should give same results as in smallx-final.nb,
    section 1-[sum]."""
    # DVCSpoints = all data from gepard's dvcs.dat
    # model parameters from DIS fit should be fixed
    DVCSpoints = data[36] + data[37] + data[38] + data[39] + \
      data[40] + data[41] + data[42] + data[43] + data[44] + \
      data[45]
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
    mGepard.g.parint.p = 0
    mGepard.g.init()
    tGepard.model.release_parameters('M02S','SKEWS','SKEWG')
    f = Fitter.FitterMinuit(DVCSpoints, tGepard)
    f.fit()
    chisq = tGepard.chisq(DVCSpoints)[0]
    assert_almost_equal(chisq, 101.094, 1)

test_gepardfitDVCS.long = 1
    
def test_hybridfitGepard():
    """Test simple hybrid fitting via gepard + DR (DR=0)"""
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
    t.m.parameters['Nv'] = 0
    t.m.parameters['C'] = 0
    t.m.g.parint.p = 0
    t.m.g.init()
    t.m.release_parameters('M02S', 'SKEWS')
    f = Fitter.FitterMinuit(fitpoints, t)
    f.fit()
    chisq = t.chisq(fitpoints)[0]
    assert_almost_equal(chisq, 5.1423271052023196, 2)

def test_hybridfitDR():
    """Test simple hybrid fitting via gepard + DR (gepard=0)"""
    fitpoints = data[31][12:14] + data[8][1:3] + data[30][2:4]
    f = Fitter.FitterMinuit(fitpoints, t)
    f.fit()
    chisq = t.chisq(fitpoints)[0]
    assert_almost_equal(chisq, 6.7638690027046771, 5)

test_hybridfitDR.newfeature = 1
    
def test_hybridfitDVCS():
    """Test fitting to HERA DVCS with DR part off.
    
    This should give same results as in smallx-final.nb,
    section 1.-[sum]."""
    # DVCSpoints = all data from gepard's dvcs.dat
    # model parameters from DIS fit should be fixed
    DVCSpoints = data[36] + data[37] + data[38] + data[39] + \
      data[40] + data[41] + data[42] + data[43] + data[44] + \
      data[45]
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
    t.m.parameters['Nv'] = 0
    t.m.parameters['C'] = 0
    t.m.g.parint.p = 0
    t.m.g.init()
    t.m.release_parameters('M02S','SKEWS','SKEWG')
    f = Fitter.FitterMinuit(DVCSpoints, t)
    f.fit()
    chisq = t.chisq(DVCSpoints)[0]
    assert_almost_equal(chisq, 101.094, 1)

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

