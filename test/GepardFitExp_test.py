
import shutil, copy, math
from nose.tools import *

import numpy as np

import utils, Model, Approach, Fitter, Data

data = utils.loaddata('data/ep2epgamma', approach=Approach.hotfixedBMK)  
data.update(utils.loaddata('data/gammastarp2gammap', approach=Approach.hotfixedBMK))

# Gepard only
mGepard = Model.ComptonGepard(ansatz='FITEXP')
tGepard = Approach.hotfixedBMK(mGepard)


# Shortcut for seting gepard model parameters
def setpar(i, val):
    mGepard.parameters[mGepard.parameters_index[i]] = val

# DVCSpoints = all data from gepard's dvcs.dat
# model parameters from DIS fit should be fixed
DVCSpoints = data[36] + data[37] + data[38] + data[39] + \
  data[40] + data[41] + data[42] + data[43] + data[44] + \
  data[45]

reducedpoints = DVCSpoints #[::6]


def test_gepardfitDVCSexp():
    """Test fitting to HERA DVCS via gepard in nlso3 exp-t model
    
    This should give same results as in 
	smallx-final.nb, section 2-[nlo]-LO."""
    setpar(11,  0.15203911208796006)
    setpar(12,  1.1575060246398083)
    setpar(13,  0.0)
    setpar(14,  1./4.63)
    setpar(15,  0.)
    setpar(16,  2.)
    setpar(17,  0.)
    setpar(18,  0.)
    setpar(19,  0.)
    setpar(22,  1.247316701070471)
    setpar(23,  0.0)
    setpar(24,  1./4.63)
    setpar(25,  0.)
    setpar(26,  2.)
    setpar(27,  0.)
    setpar(28,  0.)
    setpar(29,  0.)
    setpar(37,  0.)
    setpar(47,  0.)
    mGepard.g.parint.p = 0
    mGepard.g.init()
    tGepard.model.release_parameters('M02S','SECS','SECG')
    f = Fitter.FitterMinuit(reducedpoints, tGepard)
    f.fit()
    chisq = tGepard.chisq(reducedpoints)[0]
    assert_almost_equal(chisq/100, 97.93273/100, 1)
    assert_almost_equal(tGepard.model.parameters['M02S'], 0.16058, 2)
    assert_almost_equal(tGepard.model.parameters['SECS'], -0.18120, 2)
    assert_almost_equal(tGepard.model.parameters['SECG'], -0.86294, 2)
    tGepard.model.fix_parameters('ALL')

test_gepardfitDVCSexp.long = 1
test_gepardfitDVCSexp.extendedtesting = 1
test_gepardfitDVCSexp.gepardsuite = 1

