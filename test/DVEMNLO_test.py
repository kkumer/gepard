
import shutil, copy, math
from nose.tools import *
import numpy as np

import utils, Model, Approach, Data, Fitter

from constants import Mp, Mp2
from results import KM10b,dvmppars


def test_c1():
    """Calculate NLO DVMP coef. functions"""
	# KM10b model
	# For speed, we can stay with p=0 because this
	# doesn't influence c1 coefficients
    m = Model.ComptonGepard(p=0)
    t = Approach.BMK(m)
    t.m.parameters.update(KM10b)
    t.m.g.init()
    t.m.g.newcall = 1
    t.m.g.parint.pid = 3
    aux = t.m.g.c1dvmp(1, (0.5+1.j), 2)
    # comparing to DM's DVEM-c1_forKreso.nb
    # quark part
    assert_almost_equal(aux[0]/100, (30.3836+9.2856j)/100, 3)
    # "pure singlet" part
    assert_almost_equal(aux[1]/100, (-0.496221+3.85004j)/100, 3)
    # gluon part
    assert_almost_equal(aux[2]/100, (35.2725+34.0699j)/100, 2)

test_c1.gepardsuite = 1

def test_gepardTFFsNLO():
    """Calculate NLO DVMP TFFs for rho production at input scale."""
	# NLO model
    mnlo = Model.ComptonGepard(p=1, speed=1)
    tnlo = Approach.BMK(mnlo)
    tnlo.m.parameters.update(dvmppars)
    tnlo.m.g.parint.pid = 3
    xB = 1e-4
    pt = Data.DummyPoint(init = {'Q2':4., 't':0, 'xi':xB/(2.-xB), 'xB':xB})
    utils.fill_kinematics(pt)
    re, im =  (tnlo.m.ReHrho(pt), tnlo.m.ImHrho(pt))
    # following agrees with DM to better than percent
    assert_almost_equal(re/1e4,  5328.3678/1e4, 3)
    assert_almost_equal(im/1e4, 22676.063/1e4, 3)

test_gepardTFFsNLO.gepardsuite = 1
test_gepardTFFsNLO.long = 1
