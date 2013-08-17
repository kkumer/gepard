
# Testing gepard MSBAR NLO code

# Numbers below are produced by gepard's redNLONS
# When ratios of absolute values are taken 
#        (NLO/LO - 1)*100 = -16.0881222
# they correspond to thick green dashed line
# of panel 2 of Fig. 7 in NPB 08
# (leftmost point, also visible in
# radNLONS1.dat, at
# the start of 2nd and 4th block of numbers.


import shutil, copy, math
from nose.tools import *
import numpy as np

import utils, Model, Approach, Data, Fitter

from constants import Mp, Mp2

m = Model.ComptonGepard(process='DVCS', ansatz='NSFIT', q02=2.5)
t = Approach.hotfixedBMK(m)

# Setting gepard to values as in radNLONS.F
# (Fig. 7 of NPB08)

t.m.g.parint.nf = 4
t.m.g.astrong.asp = np.array([0.05, 0.05, 0.05])
t.m.g.parchr.fftype = np.array([c for c in 'NONSINGLET']) # array(10)
t.m.g.parchr.scheme = np.array([c for c in 'MSBDI']) # array(10)
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
setpar(11, 4./15.)

pt = Data.DummyPoint()


def test_radMSBARLONS():
    """Non-singlet LO CFF H"""
    pt.Q2 = 2.5
    pt.t = -1.
    pt.xi = 0.01
    pt.xB = 2*pt.xi/(1.+pt.xi)
    #FIXME: how to avoid this:
    try:
        del t.m.qdict[pt.Q2]
    except:
        pass
    t.m.g.parint.p = 0
    t.m.g.newcall = 1
    t.m.g.init()
    assert_almost_equal(m.ReH(pt), -1.6620214256774486)
    assert_almost_equal(m.ImH(pt), -12.604465892107740)


test_radMSBARLONS.gepardsuite = 1


def test_radMSBARNLONS():
    """Non-singlet NLO MSBAR CFF H"""
    pt.Q2 = 2.5
    pt.t = -1.
    pt.xi = 0.01
    pt.xB = 2*pt.xi/(1.+pt.xi)
    #FIXME: how to avoid this:
    try:
        del t.m.qdict[pt.Q2]
    except:
        pass
    t.m.g.parint.p = 1
    t.m.g.newcall = 1
    t.m.g.init()
    assert_almost_equal(m.ReH(pt), -1.2740328118218003)
    assert_almost_equal(m.ImH(pt), -10.591847869828122)

test_radMSBARNLONS.gepardsuite = 1

## relo = -1.6620214256774486
## imlo = -12.604465892107740
## renlo =  -1.2740328118218003
## imnlo = -10.591847869828122
## print (sqrt(renlo**2+imnlo**2)/sqrt(relo**2+imlo**2) - 1)*100
## -16.08812221458558

