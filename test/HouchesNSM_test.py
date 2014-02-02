
import shutil, copy, math, sys
from nose.tools import *
import numpy as np

import utils, Model, Approach, Data, Fitter

from constants import Mp, Mp2

m = Model.ComptonGepard(ansatz='NSMHOU', fftype='NONSINGLET', q02=2.0)
t = Approach.hotfixedBMK(m)

t.m.g.parint.nf = 4
asp0 = 0.35 / 2./ np.pi
t.m.g.astrong.asp = np.array([asp0, asp0, asp0])
t.m.g.astrong.mu02 = 2.0
t.m.g.mbcont.phi = 1.9
t.m.g.parint.pid = -2

# Seting model parameters to be as in test.F


pt = Data.DummyPoint()

Q2 = 10000.
xi = 0.1  # note that this is actually xB

def test_PDFevolINIT_NSM():
    """Les Houches (u_v-d_v) PDF at input scale"""
    pt.Q2 = 2.0
    pt.t = 0
    pt.xi = xi
    t.m.g.parint.p = 0
    t.m.g.newcall = 1
    t.m.g.init()
    assert_almost_equal(pt.xi*m.gpdHzeroQ(pt), 0.271436487, 6)

test_PDFevolINIT_NSM.gepardsuite = 1

def test_PDFevolLO_NSM():
    """Les Houches (u_v-d_v) PDF at LO"""
    pt.Q2 = Q2
    pt.t = 0
    pt.xi = xi
    t.m.g.parint.p = 0
    t.m.g.newcall = 1
    t.m.g.init()
    assert_almost_equal(pt.xi*m.gpdHzeroQ(pt), 0.288538040, 6)

test_PDFevolLO_NSM.gepardsuite = 1

def test_PDFevolNLO_NSM():
    """Les Houches (u_v-d_v) PDF at NLO"""
    pt.Q2 = Q2
    pt.t = 0
    pt.xi = xi
    t.m.g.parint.p = 1
    t.m.g.newcall = 1
    t.m.g.init()
    # FIXME: We have 4th decimal error here (0.06%), with
    # K-factor (NLO-LO)/LO equal 2.7%. Error is not
    # phenomenologically important but should be tracked down.
    assert_almost_equal(pt.xi*m.gpdHzeroQ(pt), 0.280377064, 6)

test_PDFevolNLO_NSM.gepardsuite = 1

#def test_PDFevolNLO_NSM_dbg():
    #"""Les Houches (u_v-d_v) PDF at NLO, for debugging with Pegasus"""
    #pt.Q2 = 1.e2
    #pt.t = 0
    #pt.xi = 1.e-2
    #t.m.g.parint.nf = 3
    ##t.m.g.mbcont.phi = np.pi/2.
    ##t.m.g.mbcont.c = 0.3
    #t.m.g.parint.p = 1
    #t.m.g.newcall = 1
    #t.m.g.init()
    #assert_almost_equal(pt.xi*m.gpdHzeroQ(pt), 0.0815943813, 3)

#test_PDFevolNLO_NSM_dbg.gepardsuite = 1

def test_PDFevolNNLO_NSM():
    """Les Houches u_v PDF at NNLO (not implemented in gepard)"""
    pt.Q2 = Q2
    pt.t = 0
    pt.xi = xi
    t.m.g.parint.p = 2
    t.m.g.newcall = 1
    t.m.g.init()
    assert_almost_equal(pt.xi*m.gpdHzeroQ(pt), 0.279974155E-01, 7)

test_PDFevolNNLO_NSM.gepardsuite = 1
test_PDFevolNNLO_NSM.newfeature = 1
