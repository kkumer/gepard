
import copy
from nose.tools import *
import numpy as np

import gepard as g

#import utils, Model, Approach, Fitter

#data = utils.loaddata('data/ep2epgamma')  

# TODO
#m = Model.ModelSigma()
#t = Approach.Gepard(m)

# testing data set for fits
#fitpoints = data[36]


# Setting up parameters and initializing gepard
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
    assert_almost_equal(g.f2.f2[g.parint.p], 0.65783876286, 5)


def test_gepardXDVCSt():
    """Calculate basic NLO DVCS cross section using gepard."""
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
    assert_almost_equal(aux, 6.8612738, 4)


def test_gepardfit():
    """Test fitting to H1 DVCS via gepard"""
    f = Fitter.FitterMinuit(fitpoints, t)
    f.fit()
    chisq = t.chisq(fitpoints)[0]
    assert_almost_equal(chisq, 0.2665403, 5)


test_gepardfit.newfeature = 1
