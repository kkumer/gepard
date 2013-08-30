# Three tests below, LO, NLO and NNLO, correspond to
# left-most points, x=0.01, of panel 2 of Fig. 12 in NPB 08
# paper. Their ratios are plotted there, for xi=0.01, as
#
# (NLO/LO - 1)*100 = -9.627907  (thin blue dashed line)
#
# (NNLO/NLO - 1)*100 = -1.575703 (thin red solid line)
#
# Actually, plot starts for somewhat larger xi, but
# numbers above can be seen in radNNLONS1.dat, at
# the start of 4th and 2nd block of numbers, respectively.

import shutil, copy, math
from nose.tools import *
import numpy as np

import utils, Model, Approach, Data, Fitter

from constants import Mp, Mp2

m = Model.ComptonGepard(ansatz='NSFIT', scheme='CSBAR', q02=2.5)
t = Approach.hotfixedBMK(m)

# Setting gepard to values as in radNNLONS.F
# (Fig. 12 of NPB08)

t.m.g.parint.nf = 4
t.m.g.parint.pid = 1
t.m.g.astrong.asp = np.array([0.05, 0.05, 0.05])
t.m.g.parchr.fftype = np.array([c for c in 'NONSINGLET']) # array(10)
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

pt = Data.DummyPoint()


def test_radLONS():
    """Non-singlet LO CFF H"""
    pt.Q2 = 10.
    pt.t = -0.25
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
    assert_almost_equal(m.ReH(pt), -3.3434664508535783)
    assert_almost_equal(m.ImH(pt), 3.274603763172275)
   
test_radLONS.gepardsuite = 1

def test_radNLONS():
    """Non-singlet NLO CFF H"""
    pt.Q2 = 10.
    pt.t = -0.25
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
    assert_almost_equal(m.ReH(pt), -2.8809502819615411)
    assert_almost_equal(m.ImH(pt), 3.096381114053143)

test_radNLONS.gepardsuite = 1

def test_radNNLONS():
    """Non-singlet NNLO CFF H"""
    pt.Q2 = 10.
    pt.t = -0.25
    pt.xi = 0.01
    pt.xB = 2*pt.xi/(1.+pt.xi)
    try:
        del t.m.qdict[pt.Q2]
    except:
        pass
    t.m.g.parint.p = 2
    t.m.g.newcall = 1
    t.m.g.init()
    assert_almost_equal(m.ReH(pt), -2.8495779492285682)
    assert_almost_equal(m.ImH(pt), 3.0344836334101961)

test_radNNLONS.gepardsuite = 1

## FORTRAN:
## lo = (-3.3434664508535783     ,  3.2746037631722751  ) 
## nlo = ( -2.8809502819615411     ,  3.0963811140531430  ) 
## nnlo = ( -2.8495779492285682     ,  3.0344836334101961  )
## python:
## relo = -3.3434664508535783
## imlo = 3.274603763172275
## renlo = -2.8809502819615411
## imnlo = 3.096381114053143
## rennlo = -2.8495779492285682
## imnnlo = 3.0344836334101961
## print (sqrt(renlo**2+imnlo**2)/sqrt(relo**2+imlo**2) - 1)*100
## print (sqrt(rennlo**2+imnnlo**2)/sqrt(renlo**2+imnlo**2) - 1)*100
## 
## -9.62790688897865
## -1.57570328069094

