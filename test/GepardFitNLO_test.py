
import shutil, copy, math
from nose.tools import *

import numpy as np

import utils, Model, Approach, Fitter

data = utils.loaddata('data/ep2epgamma', approach=Approach.hotfixedBMK)  
data.update(utils.loaddata('data/gammastarp2gammap', approach=Approach.hotfixedBMK))
data.update(utils.loaddata('data/DIS', approach=Approach.BMK))

# Gepard only
mGepard = Model.ComptonGepard(scheme='MSBAR')
tGepard = Approach.hotfixedBMK(mGepard)
# removing some limits for compatibility with old Gepard
del mGepard.parameters['limit_M02S']
del mGepard.parameters['limit_M02G']

# Hybrid: Gepard+DR (can reuse above Gepard)
mDRsea = Model.ComptonModelDRsea()
m = Model.HybridDipole(mGepard, mDRsea)
t = Approach.hotfixedBMK(m)

# Shortcut for seting gepard model parameters
def setpar(i, val):
    mGepard.parameters[mGepard.parameters_index[i]] = val

# DISpoints = all data from gepard's dis.dat
DISpoints = data[201] + data[202] + data[203] + data[204] + \
            data[205] + data[206] + data[207] + data[208] + \
            data[209] + data[210] + data[211] + data[212]

# DVCSpoints = all data from gepard's dvcs.dat
# model parameters from DIS fit should be fixed
DVCSpoints = data[36] + data[37] + data[38] + data[39] + \
  data[40] + data[41] + data[42] + data[43] + data[44] + \
  data[45]


def test_gepardfitDISNLO():
    """Test fitting to H1 DIS F2 at NLO
    
    This should give same results as in smallx-final.nb,
    section 1- DIS fits - NLO."""
    setpar(11,  0.15)
    setpar(12,  1.0)
    setpar(13,  0.15)
    setpar(14,  1.)
    setpar(15,  0.)
    setpar(16,  2.)
    setpar(17,  0.)
    setpar(18,  0.)
    setpar(19,  0.)
    setpar(22,  1.1)
    setpar(23,  0.15)
    setpar(24,  0.7)
    setpar(25,  0.)
    setpar(26,  2.)
    setpar(27,  0.)
    setpar(28,  0.)
    setpar(29,  0.)
    setpar(37,  0.)
    setpar(47,  0.)
    mGepard.g.parint.p = 1
    mGepard.g.init()
    tGepard.model.fix_parameters('ALL')
    tGepard.model.release_parameters('NS','AL0S','AL0G')
    f = Fitter.FitterMinuit(DISpoints, tGepard)
    f.fit()
    chisq = tGepard.chisq(DISpoints)[0]
    assert_almost_equal(chisq/100, 71.61841134496/100, 5)
    # BTW: NNLODISChis = {80.97373402543556, 82}
    assert_almost_equal(mGepard.parameters['AL0G'], 1.0990016751643479, 3)

test_gepardfitDISNLO.long = 1
test_gepardfitDISNLO.gepardsuite = 1

def test_gepardfitDVCSNLO():
    """Test fitting to HERA DVCS via gepard nl-SO3 at NLO
    
    This should give same results as in smallx-final.nb,
    section 1-[nlo] - NLO - MSBAR."""
    setpar(11,  0.1678)
    setpar(12,  1.12835)
    setpar(13,  0.15)
    setpar(14,  1.)
    setpar(15,  0.)
    setpar(16,  2.)
    setpar(17,  0.)
    setpar(18,  0.)
    setpar(19,  0.)
    setpar(22,  1.099)
    setpar(23,  0.15)
    setpar(24,  0.7)
    setpar(25,  0.)
    setpar(26,  2.)
    setpar(27,  0.)
    setpar(28,  0.)
    setpar(29,  0.)
    setpar(37,  0.)
    setpar(47,  0.)
    mGepard.g.parint.p = 1
    mGepard.g.init()
    tGepard.model.fix_parameters('ALL')
    tGepard.model.release_parameters('M02S','SECS','SECG')
    f = Fitter.FitterMinuit(DVCSpoints, tGepard)
    f.fit()
    chisq = tGepard.chisq(DVCSpoints)[0]
    assert_almost_equal(chisq/1000, 101.581/1000, 2)
    assert_almost_equal(mGepard.parameters['M02S'], 0.590339, 2)

test_gepardfitDVCSNLO.long = 1
test_gepardfitDVCSNLO.gepardsuite = 1

