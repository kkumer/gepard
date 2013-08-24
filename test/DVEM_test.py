
import shutil, copy, math
from nose.tools import *
import numpy as np

import utils, Model, Approach, Data, Fitter

from constants import Mp, Mp2
from results import KM10b, dvmppars

# KM10b model
m = Model.ComptonGepard(p=0)
t = Approach.BMK(m)
t.m.parameters.update(KM10b)

# generic LO model from big DVMP draft
mlo = Model.ComptonGepard(p=0)
tlo = Approach.BMK(mlo)
tlo.m.parameters.update(dvmppars)


def test_gepardTFFs():
    """Calculate LO DVMP TFFs for rho production at input scale."""
    xB = 1e-4
    pt = Data.DummyPoint(init = {'Q2':4., 't':0, 'xi':xB/2., 'xB':xB})
    utils.fill_kinematics(pt)
    re, im =  (tlo.m.ReHrho(pt), tlo.m.ImHrho(pt))
    # following agrees with DM to best than percent
    assert_almost_equal(re/1e4, 4766.8993/1e4, 3)
    assert_almost_equal(im/1e4, 12395.53/1e4, 3)

test_gepardTFFs.gepardsuite = 1

def test_gepardTFFsEvol():
    """Calculate LO DVMP TFFs for rho production + evolution."""
    pt = Data.DummyPoint()
    pt.W = 75.
    pt.Q2 = 6.6
    pt.t = -0.025
    pt.xi = pt.Q2 / ( 2.0 * pt.W * pt.W + pt.Q2)
    pt.xB = 2*pt.xi/(1.+pt.xi)
    t.m.g.init()
    t.m.g.newcall = 1
    aux = t.m.ReHrho(pt)
    assert_almost_equal(aux/100, 285.15512/100, 2)
    aux = t.m.ImHrho(pt)
    assert_almost_equal(aux/100, 511.39622404/100, 2)

test_gepardTFFsEvol.gepardsuite = 1

def test_gepardXrhot():
    """Calculate LO DVMP cross section d sigma / dt"""
    pt = Data.DummyPoint()
    pt.process = 'gammastarp2Mp'
    pt.W = 75.
    pt.Q2 = 6.6
    pt.t = -0.025
    pt.xi = pt.Q2 / ( 2.0 * pt.W * pt.W + pt.Q2)
    pt.xB = 2*pt.xi/(1.+pt.xi)
    t.m.g.init()
    t.m.g.newcall = 1
    aux = t.X(pt)
    assert_almost_equal(aux/1000., 1212.62165/1000., 2)

test_gepardXrhot.gepardsuite = 1

