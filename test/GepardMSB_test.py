
# Testing gepard MSBAR NLO code - singlet case

# Numbers below are produced by gepard's ex/evolut


import shutil, copy, math
from nose.tools import *
import numpy as np

import utils, Model, Approach, Data, Fitter

from constants import Mp, Mp2

m = Model.ComptonGepard(process='DVCS', ansatz='FITBP', q02=2.5)
t = Approach.hotfixedBMK(m)
# Setting gepard to values as in radNLONS.F
# (Fig. 7 of NPB08)
t.m.g.parint.nf = 4
t.m.g.astrong.asp = np.array([0.05, 0.05, 0.05])
t.m.g.parchr.fftype = np.array([c for c in 'SINGLET   ']) # array(10)
t.m.g.parchr.scheme = np.array([c for c in 'MSBAR']) # array(?)
t.m.g.mbcont.phi = 1.9

# Seting model parameters to be as in test.F
def setpar(i, val):
    t.m.parameters[t.m.parameters_index[i]] = val

MassP = 0.938272 # proton mass
#setpar(11, 0.2)
setpar(11, 0.0) # pure valence
setpar(12, 1.1)
setpar(13, 0.15)
setpar(14, (2.0 * MassP)**2)
setpar(15, MassP**2)
setpar(16, 3.0)
setpar(21, 0.5)
setpar(22, 1.0)
setpar(23, 0.15)
setpar(24, (2.0 * MassP)**2)
setpar(25, MassP**2)
setpar(26, 2.0)
setpar(31, 2.0)
setpar(32, 0.5)
setpar(33, 1.0)
setpar(34, (2.0 * MassP)**2)
setpar(35, MassP**2)
setpar(36, 1.0)
setpar(41, 1.0)
setpar(42, 0.5)
setpar(43, 1.0)
setpar(44, (2.0 * MassP)**2)
setpar(45, MassP**2)
setpar(46, 1.0)

#  'hard' ansatz
setpar(21, 0.4)
setpar(22, 1.1 + 0.05)
setpar(11, 2./3. - 0.4)

pt = Data.DummyPoint()
pt.Q2 = 25.
pt.t = -0.25
pt.xi = 1.e-5
pt.xB = 2*pt.xi/(1.+pt.xi)


def test_radMSBARLO():
    """Singlet LO CFF H evolved"""
    #FIXME: how to avoid this:
    try:
        del t.m.qdict[pt.Q2]
    except:
        pass
    t.m.g.parint.p = 0
    t.m.g.newcall = 1
    t.m.g.init()
    assert_almost_equal(m.ReH(pt)/1e5, 251460.03959908773/1e5)
    assert_almost_equal(m.ImH(pt)/1e5, 1015357.1865059549/1e5)


test_radMSBARLO.gepardsuite = 1


def test_radMSBARNLO():
    """Singlet NLO MSBAR CFF H evolved"""
    #FIXME: how to avoid this:
    try:
        del t.m.qdict[pt.Q2]
    except:
        pass
    t.m.g.parint.p = 1
    t.m.g.newcall = 1
    t.m.g.init()
    assert_almost_equal(m.ReH(pt)/1e5, 142867.21556625995/1e5)
    assert_almost_equal(m.ImH(pt)/1e5, 653095.26655367797/1e5)

test_radMSBARNLO.gepardsuite = 1
