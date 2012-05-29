
import shutil, copy, math, sys
from nose.tools import *
import numpy as np

import utils, Model, Approach, Data, Fitter

from constants import Mp, Mp2

m = Model.ComptonGepard(ansatz='HOUCHE', q02=2.0)
t = Approach.hotfixedBMK(m)

# Setting gepard to values as in radNNLONS.F
# (Fig. 12 of NPB08)

t.m.g.parint.nf = 4
asp0 = 0.35 / 2./ np.pi
t.m.g.astrong.asp = np.array([asp0, asp0, asp0])
t.m.g.astrong.mu02 = 2.0
t.m.g.mbcont.phi = 1.9
t.m.g.parchr.fftype = np.array([c for c in 'SINGLET   ']) # array(10)
t.m.g.parchr.process = np.array([c for c in 'DVCSZG'])  # array(6)

# Seting model parameters to be as in test.F


pt = Data.DummyPoint()

Q2 = 10000.
xi = 0.1

def test_DISNLO():
    """Les Houches gluon PDF at NLO"""
    pt.Q2 = Q2
    pt.t = 0
    pt.xi = xi
    pt.xB = 2*pt.xi/(1.+pt.xi)
    #FIXME: how to avoid this:
    try:
        del t.m.qdict[pt.Q2]
    except:
        pass
    t.m.g.parint.p = 1
    t.m.g.newcall = 1
    t.m.g.init()
    assert_almost_equal(pt.xi*m.ImH(pt)/np.pi, 0.9028145873)

test_DISNLO.gepardsuite = 1

def test_DISNNLO():
    """Les Houches gluon PDF at NNLO"""
    pt.Q2 = Q2
    pt.t = 0
    pt.xi = xi
    pt.xB = 2*pt.xi/(1.+pt.xi)
    #FIXME: how to avoid this:
    try:
        del t.m.qdict[pt.Q2]
    except:
        pass
    t.m.g.parint.p = 2
    t.m.g.newcall = 1
    t.m.g.init()
    assert_almost_equal(pt.xi*m.ImH(pt)/np.pi, 0.9068212068)

test_DISNNLO.gepardsuite = 1
