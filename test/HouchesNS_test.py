
import shutil, copy, math, sys
from nose.tools import *
import numpy as np

import utils, Model, Approach, Data, Fitter

from constants import Mp, Mp2

m = Model.ComptonGepard(ansatz='NSHOUC', fftype='NONSINGLET', q02=2.0)
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

def test_PDFevolINIT_NS():
    """Les Houches u_v PDF at input scale"""
    pt.Q2 = 2.0
    pt.t = 0
    pt.xi = xi
    t.m.g.parint.p = 0
    t.m.g.newcall = 1
    t.m.g.init()
    assert_almost_equal(pt.xi*m.gpdHzeroQ(pt), 0.590079319, 6)

test_PDFevolINIT_NS.gepardsuite = 1

def test_PDFevolLO_NS():
    """Les Houches u_v PDF at LO"""
    pt.Q2 = Q2
    pt.t = 0
    pt.xi = xi
    t.m.g.parint.p = 0
    t.m.g.newcall = 1
    t.m.g.init()
    assert_almost_equal(pt.xi*m.gpdHzeroQ(pt), 0.572672517, 6)

test_PDFevolLO_NS.gepardsuite = 1

def test_PDFevolNLO_NS():
    """Les Houches u_v PDF at NLO"""
    pt.Q2 = Q2
    pt.t = 0
    pt.xi = xi
    t.m.g.parint.p = 1
    t.m.g.newcall = 1
    t.m.g.init()
    assert_almost_equal(pt.xi*m.gpdHzeroQ(pt), 0.553322489, 2)

test_PDFevolNLO_NS.gepardsuite = 1

def test_PDFevolNNLO_NS():
    """Les Houches u_v PDF at NNLO"""
    pt.Q2 = Q2
    pt.t = 0
    pt.xi = xi
    t.m.g.parint.p = 2
    t.m.g.newcall = 1
    t.m.g.init()
    assert_almost_equal(pt.xi*m.gpdHzeroQ(pt), 0.551455516, 4)

test_PDFevolNNLO_NS.gepardsuite = 1
test_PDFevolNNLO_NS.extendedtesting = 1
