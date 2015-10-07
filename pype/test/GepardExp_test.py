
import shutil, copy, math
from nose.tools import *

import numpy as np

import Model, Approach, Fitter, Data

# Gepard only
mGepard = Model.ComptonGepard(ansatz='FITEXP', scheme='CSBAR')
tGepard = Approach.hotfixedBMK(mGepard)

# Setting gepard to test values
tGepard.m.g.parint.p = 1
tGepard.m.g.parint.nf = 3
tGepard.m.g.parint.pid = 1
tGepard.m.g.astrong.asp = np.array([0.05, 0.05, 0.05])
tGepard.m.g.parflt.q02 = 1.0
#tGepard.m.g.parchr.ansatz = np.array([c for c in 'TEST  ']) # array(6)

# Shortcut for seting gepard model parameters
def setpar(i, val):
    mGepard.parameters[mGepard.parameters_index[i]] = val


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

def test_gepardXDVCStexp():
    """Calculate LO DVCS partial cross section using exp-t gepard (no evolution).
    This is equal to fotran test.F, with ANSATZ=FITEXP"""
    pt.W = 82.
    pt.Q2 = 1.
    pt.t = 0.0
    pt.xi = pt.Q2 / ( 2.0 * pt.W * pt.W + pt.Q2)
    pt.xB = 2.*pt.xi/(1+pt.xi)
    tGepard.m.g.parint.p = 0
    tGepard.m.g.init()
    aux = tGepard.X(pt)
    assert_almost_equal(aux/1e4, 6941.1676907979203/1e4, 5)

test_gepardXDVCStexp.gepardsuite = 1

def test_gepardXDVCSevolexp():
    """Calculate NLO DVCS cross section using exp-t gepard (+evolution).
    This is equal to fotran test.F, with ANSATZ=FITEXP, and W=3.5"""
    pt.W = 55          
    pt.Q2 = 3.          
    pt.xi = pt.Q2 / ( 2.0 * pt.W * pt.W + pt.Q2)
    pt.xB = 2.*pt.xi/(1+pt.xi)
    del pt.t
    tGepard.m.g.parint.p = 1
    tGepard.m.g.init()
    aux = tGepard.X(pt)
    assert_almost_equal(aux, 98.538728060134133, 5)

test_gepardXDVCSevolexp.gepardsuite = 1
