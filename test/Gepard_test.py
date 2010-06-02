
import copy, math
from nose.tools import *
import numpy as np

import gepard as g

#import utils, Model, Approach, Fitter
import utils, Model, Approach, Data

from constants import Mp, Mp2

#data = utils.loaddata('data/ep2epgamma')  

# TODO
m = Model.ComptonGepard()
t = Approach.hotfixedBMK(m)

# testing data set for fits
#fitpoints = data[36]

pt = Data.DummyPoint()
pt.exptype = 'collider'
pt.in1particle = 'e'
pt.in1charge = 1
pt.in1energy = 5.
pt.in1polarizationvector = 'L'
pt.in1polarization = 1.0
pt.in2particle = 'p'
pt.in2energy = 250.
pt.s = 2 * pt.in1energy * (pt.in2energy + math.sqrt(
    pt.in2energy**2 - Mp2)) + Mp2
pt.W = 82.
pt.Q2 = 1.
pt.t = 0.0

utils.fill_kinematics(pt)

## Setting up parameters and initializing gepard
g.readpar()
g.parchr.fftype = np.array([c for c in 'SINGLET   '])
g.init()

g.par.par[21-1] = 0.4
g.par.par[11-1] = 2./3. - 0.4
g.par.par[12-1] = 1.1
g.par.par[13-1] = 0.25
g.par.par[14-1] = 1.1
g.par.par[22-1] = 1.2
g.par.par[23-1] = 0.25
g.par.par[24-1] = 1.2

#   turn off second PW
g.par.par[17-1] = 0.0
g.par.par[27-1] = 0.0

#   turn off GPD E
g.par.par[18-1] = 0.0
g.par.par[28-1] = 0.0

#   turn off xi-dependence
g.par.par[19-1] = 0.0
g.par.par[29-1] = 0.0


def test_gepardF2():
    """Calculate LO DIS F2 using gepard."""
    g.parint.p = 0
    g.parchr.process = np.array([c for c in 'DIS   '])
    g.kinematics.del2 = 0.0

    g.init()
    g.getmbgpd()
    g.hgrid.hgrid[0] = g.mbgpd.mbgpd

    g.kinematics.xi = 0.002
    g.kinematics.q2 = 1.0
    g.nqs.nqsdis = 1
    g.qs.qsdis[0] = 1.0
    g.evolc(1,1)
    g.f2f()
    assert_almost_equal(g.f2.f2[g.parint.p], 0.65784282163635943)

def test_gepardXDVCSevol():
    """Calculate basic NLO DVCS cross section using gepard (+evolution)."""
    g.parint.p = 1
    g.parchr.process = np.array([c for c in 'DVCS  '])
    g.init()

    g.kinematics.w2 = 3.5**2
    g.kinematics.q2 = 3.0
    g.kinematics.xi = g.kinematics.q2 / ( 2.0 * g.kinematics.w2 + g.kinematics.q2)
    g.mt.mtind = 0
    g.nqs.nqs = 2
    g.qs.qs[1] = g.kinematics.q2
    g.evolc(1,2)
    aux = g.sigma()
    assert_almost_equal(aux, 6.8612469682766850, 5)


def test_gepardXDVCSt():
    aux = t.PartialCrossSection(pt)
    assert_almost_equal(aux, 5607.5998541187819, 2)

#def test_gepardXDVCSt():
#    g.parint.p = 0
#    g.parchr.process = np.array([c for c in 'DVCS  '])
#    g.init()
#    
#    g.kinematics.w2 = 82.0**2
#    g.kinematics.q2 = 1.0
#    g.kinematics.xi = g.kinematics.q2 / ( 2.0 * g.kinematics.w2 + g.kinematics.q2)

#    g.nqs.nqs = 1
#    g.qs.qs[0] = g.kinematics.q2
#    g.evolc(1,1)

#    g.mt.mtind = 0

#    aux = g.parsigma()
#    assert_almost_equal(aux, 5607.5998541187819, 2)


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
