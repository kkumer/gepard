
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

def test_c1():
    """Calculate NLO DVMP coef. functions"""
    t.m.g.init()
    t.m.g.newcall = 1
    aux = t.m.g.cdvemf((0.5+1.j), 2)
    # comparing to DM's DVEM-c1_forKreso.nb
    # quark part
    assert_almost_equal(aux[0]/100, (30.3836+9.2856j)/100, 3)
    # "pure singlet" part
    assert_almost_equal(aux[1]/100, (-0.496221+3.85004j)/100, 3)
    # gluon part
    assert_almost_equal(aux[2]/100, (35.2725+34.0699j)/100, 2)

test_c1.gepardsuite = 1
