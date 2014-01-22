
import shutil, copy, math, sys
from nose.tools import *
import numpy as np

import utils, Model, Approach, Data, Fitter

from constants import Mp, Mp2

m = Model.ComptonGepard(ansatz='HOUCHE', fftype='SINGLET', q02=2.0)
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

def test_PDFevolINIT():
    """Les Houches gluon PDF at input scale"""
    pt.Q2 = 2.0
    pt.t = 0
    pt.xi = xi
    t.m.g.parint.p = 0
    t.m.g.newcall = 1
    t.m.g.init()
    assert_almost_equal(m.gpdHzeroG(pt), 1.26375109, 6)

test_PDFevolINIT.gepardsuite = 1
test_PDFevolINIT.extendedtesting = 1


def test_PDFevolLO():
    """Les Houches gluon PDF at LO"""
    pt.Q2 = Q2
    pt.t = 0
    pt.xi = xi
    t.m.g.parint.p = 0
    t.m.g.newcall = 1
    t.m.g.init()
    assert_almost_equal(m.gpdHzeroG(pt), 0.887657279, 6)

test_PDFevolLO.gepardsuite = 1

def test_PDFevolNLO():
    """Les Houches gluon PDF at NLO"""
    pt.Q2 = Q2
    pt.t = 0
    pt.xi = xi
    t.m.g.parint.p = 1
    t.m.g.newcall = 1
    t.m.g.init()
    assert_almost_equal(m.gpdHzeroG(pt), 0.9028145873)

test_PDFevolNLO.gepardsuite = 1

def test_PDFevolNNLO():
    """Les Houches gluon PDF at NNLO"""
    pt.Q2 = Q2
    pt.t = 0
    pt.xi = xi
    t.m.g.parint.p = 2
    t.m.g.newcall = 1
    t.m.g.init()
    assert_almost_equal(m.gpdHzeroG(pt), 0.9068212068)

test_PDFevolNNLO.gepardsuite = 1
test_PDFevolNNLO.extendedtesting = 1

def test_PDFevolNNLO_Q10():
    """Les Houches gluon PDF at NNLO. Short evolution."""
    pt.Q2 = 10.
    pt.t = 0
    pt.xi = xi
    t.m.g.parint.p = 2
    t.m.g.newcall = 1
    t.m.g.init()
    assert_almost_equal(m.gpdHzeroG(pt)/10, 1.26362026/10, 5)

test_PDFevolNNLO_Q10.gepardsuite = 1
test_PDFevolNNLO_Q10.extendedtesting = 1

