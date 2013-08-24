
import shutil, copy, math
from nose.tools import *
import numpy as np

import utils, Model, Approach, Data, Fitter

from constants import Mp, Mp2

m = Model.ComptonGepard()
t = Approach.hotfixedBMK(m)

# Setting gepard to test values
t.m.g.parint.p = 1
t.m.g.parint.nf = 3
t.m.g.astrong.asp = np.array([0.05, 0.05, 0.05])
t.m.g.parflt.q02 = 1.0
t.m.g.parchr.ansatz = np.array([c for c in 'TEST  ']) # array(6)

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
    pt.xB = 2*pt.xi/(1.+pt.xi)
    t.m.g.parint.p = 0
    t.m.g.parint.pid = 1  # CSBAR
    t.m.g.init()
    aux = t.X(pt)
    assert_almost_equal(aux, 5608.42804256, 2)

test_gepardXDVCSt.gepardsuite = 1

def test_gepardXDVCSevol():
    """Calculate NLO DVCS cross section using gepard (+evolution)."""
    pt.W = 55          
    pt.Q2 = 3.          
    pt.xi = pt.Q2 / ( 2.0 * pt.W * pt.W + pt.Q2)
    pt.xB = 2*pt.xi/(1.+pt.xi)
    pt.t = 0  #in case this test is run alone
    del pt.t
    t.m.g.parint.p = 1
    t.m.g.init()
    aux = t.X(pt)
    assert_almost_equal(aux, 32.29728102)
    # To get complete agreement with Fortran take
    # (slower) tquadrature = quadSciPy10 and:
    # set pt.W=3.5  and use XDVCStApprox and
    # assert_almost_equal(aux, 6.8612469682766850, 5)

test_gepardXDVCSevol.gepardsuite = 1

def test_hybrid():
    """GepardDR with switched-off DR should be same as Gepard."""
    mDRsea = Model.ComptonModelDRsea()
    mBoth = Model.HybridDipole(m, mDRsea)
    tBoth = Approach.hotfixedBMK(mBoth)
    tBoth.m.parameters['Nv'] = 0
    tBoth.m.parameters['Nsea'] = 0
    pt.W = 82.
    pt.Q2 = 1.
    pt.t = 0.0
    pt.xi = pt.Q2 / ( 2.0 * pt.W * pt.W + pt.Q2)
    pt.xB = 2*pt.xi/(1.+pt.xi)
    tBoth.m.g.parint.p = 0
    tBoth.m.g.parint.pid = 1
    tBoth.m.g.init()
    aux = tBoth.X(pt)/1.e4
    assert_almost_equal(aux, 0.56084280, 3)

def test_hybridopt():
    """GepardDR optimized with switched-off DR should be same as Gepard."""
    mDRsea = Model.ComptonModelDRsea(optimization=True)
    mBoth = Model.HybridDipole(m, mDRsea)
    tBoth = Approach.hotfixedBMK(mBoth)
    tBoth.m.parameters['Nv'] = 0
    tBoth.m.DR.ndparameters[0] = 0
    tBoth.m.parameters['Nsea'] = 0
    tBoth.m.DR.ndparameters[6] = 0
    pt.W = 82.
    pt.Q2 = 1.
    pt.t = 0.0
    pt.xi = pt.Q2 / ( 2.0 * pt.W * pt.W + pt.Q2)
    pt.xB = 2*pt.xi/(1.+pt.xi)
    tBoth.m.g.parint.p = 0
    tBoth.m.g.parint.pid = 1
    tBoth.m.g.init()
    aux = tBoth.X(pt)/1.e4
    assert_almost_equal(aux, 0.56075998541187819, 3)

