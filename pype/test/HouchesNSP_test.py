
import shutil, copy, math, sys
from nose.tools import *
import numpy as np

import utils, Model, Approach, Data, Fitter

from consts import Mp, Mp2

m = Model.ComptonGepard(ansatz='NSPHOU', q02=2.0)
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

def test_PDFevolINIT_NSP():
    """Les Houches (u+ubar-d-dbar) PDF at input scale"""
    pt.Q2 = 2.0
    pt.t = 0
    pt.xi = xi
    t.m.g.parint.p = 0
    t.m.g.newcall = 1
    t.m.g.init()
    assert_almost_equal(pt.xi*m.gpdHzeroQ(pt), 0.245479296, 6)

test_PDFevolINIT_NSP.gepardsuite = 1

def test_PDFevolLO_NSP():
    """Les Houches (u+ubar-d-dbar) PDF  LO evolution"""
    pt.Q2 = Q2
    pt.t = 0
    pt.xi = xi
    t.m.g.parint.p = 0
    t.m.g.newcall = 1
    t.m.g.init()
    assert_almost_equal(pt.xi*m.gpdHzeroQ(pt), 0.267597475, 6)

test_PDFevolLO_NSP.gepardsuite = 1

def test_PDFevolNLO_NSP():
    """Les Houches (u+ubar-d-dbar) PDF NLO evolution"""
    pt.Q2 = Q2
    pt.t = 0
    pt.xi = xi
    t.m.g.parint.p = 1
    t.m.g.newcall = 1
    t.m.g.init()
    assert_almost_equal(pt.xi*m.gpdHzeroQ(pt), 2.60426989E-01, 6)

test_PDFevolNLO_NSP.gepardsuite = 1

#def test_PDFevolNLO_NSP_dbg():
    #"""Les Houches (u+ubar-d-dbar) PDF NLO evolution"""
    #pt.Q2 = 1.e2
    #pt.t = 0
    #pt.xi = 1.e-2
    #t.m.g.parint.nf = 3
    #t.m.g.mbcont.phi = np.pi/2.
    #t.m.g.mbcont.c = 0.3
    #t.m.g.parint.p = 1
    #t.m.g.newcall = 1
    #t.m.g.init()
    #assert_almost_equal(pt.xi*m.gpdHzeroQ(pt), 7.26024364E-02, 3)


def test_PDFevolNNLO_NSP():
    """Les Houches (u+ubar-d-dbar) PDF NNLO evolution"""
    pt.Q2 = Q2
    pt.t = 0
    pt.xi = xi
    t.m.g.parint.p = 2
    t.m.g.newcall = 1
    t.m.g.init()
    assert_almost_equal(pt.xi*m.gpdHzeroQ(pt), 2.59924771E-01, 3)

test_PDFevolNNLO_NSP.gepardsuite = 1
# Should work but error is larger than expected
test_PDFevolNNLO_NSP.newfeature = 1
