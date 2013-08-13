
import shutil, copy, math
from nose.tools import *
import numpy as np

import utils, Model, Approach, Data, Fitter

from constants import Mp, Mp2
from results import KM10b

m = Model.ComptonGepard(process='DVMP')
t = Approach.BMK(m)
t.m.parameters.update(KM10b)

pt = Data.DummyPoint()

def test_gepardTFFs():
    """Calculate LO DVMP TFFs for rho production."""
    pt.W = 75.
    pt.Q2 = 6.6
    pt.t = -0.025
    pt.xi = pt.Q2 / ( 2.0 * pt.W * pt.W + pt.Q2)
    pt.xB = 2*pt.xi/(1.+pt.xi)
    t.m.g.init()
    t.m.g.newcall = 1
    aux = t.m.ReHrho(pt)
    assert_almost_equal(aux/100, 276.163/100, 2)
    aux = t.m.ImHrho(pt)
    assert_almost_equal(aux/100, 495.104/100, 2)

test_gepardTFFs.gepardsuite = 1

def test_gepardXrhot():
    """Calculate LO DVMP cross section d sigma / dt"""
    pt.W = 75.
    pt.Q2 = 6.6
    pt.t = -0.025
    pt.xi = pt.Q2 / ( 2.0 * pt.W * pt.W + pt.Q2)
    pt.xB = 2*pt.xi/(1.+pt.xi)
    t.m.g.init()
    t.m.g.newcall = 1
    aux = t.Xrhot(pt)
    assert_almost_equal(aux/1000., 1137.124/1000., 2)

test_gepardXrhot.gepardsuite = 1

