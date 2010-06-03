
import copy, math
from nose.tools import *
import numpy as np

import gepard as g

#import utils, Model, Approach, Fitter
import utils, Model, Approach, Data

from constants import Mp, Mp2

#data = utils.loaddata('data/ep2epgamma')  
# testing data set for fits
#fitpoints = data[36]

m = Model.ComptonGepard()
t = Approach.hotfixedBMK(m)

# Seting model parameters to be as in test.F
def setpar(i, val):
    t.m.parameters[t.m.parameters_index[i]] = val

setpar(21, 0.4)
setpar(11, 2./3. - 0.4)
setpar(12, 1.1)
setpar(13, 0.25)
setpar(14, 1.1)
setpar(22, 1.2)
setpar(23, 0.25)
setpar(24, 1.2)
setpar(17, 0.0)
setpar(27, 0.0)
setpar(18, 0.0)
setpar(28, 0.0)
setpar(19, 0.0)
setpar(29, 0.0)


pt = Data.DummyPoint()


def test_gepardXDVCSt():
    """Calculate LO DVCS partial cross section using gepard (no evolution)."""
    pt.W = 82.
    pt.Q2 = 1.
    pt.t = 0.0
    pt.xi = pt.Q2 / ( 2.0 * pt.W * pt.W + pt.Q2)
    t.m.g.parint.p = 0
    t.m.g.init()
    aux = t.PartialCrossSection(pt)
    assert_almost_equal(aux, 5607.5998541187819, 2)

def test_gepardXDVCSevol():
    """Calculate basic NLO DVCS cross section using gepard (+evolution)."""
    pt.W = 3.5
    pt.Q2 = 3.
    pt.t = 0.0
    pt.xi = pt.Q2 / ( 2.0 * pt.W * pt.W + pt.Q2)
    t.m.g.parint.p = 1
    t.m.g.init()
    aux = t.TotalCrossSection(pt)
    assert_almost_equal(aux, 6.8606314494041793)
    # To get complete agreement with Fortran take
    # (slower) tquadrature = quadSciPy10 and:
    # assert_almost_equal(aux, 6.8612469682766850, 5)

#test_gepardXDVCSevol.newfeature = 1

def test_gepardfitsimple():
    """Test simple fitting to H1 DVCS (one dataset) via gepard"""
    # fitpoints = data[36] dataset
    f = Fitter.FitterMinuit(fitpoints, t)
    f.fit()
    chisq = t.chisq(fitpoints)[0]
    assert_almost_equal(chisq, 0.2665403, 5)

test_gepardfitsimple.newfeature = 1

def test_gepardfitDIS():
    """Test fitting to H1 DIS via gepard"""
    # DISpoints = all data from gepard's dis.dat
    t.model.release_parameters('NS', 'AL0S', 'AL0G')
    f = Fitter.FitterMinuit(DISpoints, t)
    f.fit()
    chisq = t.chisq(fitpoints)[0]
    assert_almost_equal(chisq, 49.7312, 5)

test_gepardfitDIS.newfeature = 1

def test_gepardfitDVCS():
    """Test fitting to H1 DVCS via gepard"""
    # DVCSpoints = all data from gepard's dvcs.dat
    # model parameters from DIS fit should be fixed
    t.model.fix_parameters('NS', 'AL0S', 'AL0G')
    t.model.release_parameters('M02S','SKEWS','SKEWG')
    f = Fitter.FitterMinuit(DVCSpoints, t)
    f.fit()
    chisq = t.chisq(fitpoints)[0]
    assert_almost_equal(chisq, 101.094, 5)

test_gepardfitDVCS.newfeature = 1
