# Three tests below, LO, NLO and NNLO, correspond to
# xi=1e-4?? points of panel 2 of Fig. 13 in NPB 08
# paper. Their ratios are plotted there, as
#
# (NLO/LO - 1)*100 = -9.627907  (thin blue dashed line)
#
# (NNLO/NLO - 1)*100 = -1.575703 (thin red solid line)
#
# Actually, plot starts for somewhat larger xi, but
# numbers can be seen in radNNLONS1.dat, at
# the start of 2nd and 4th block of numbers. ???

import shutil, copy, math
from nose.tools import *
import numpy as np

import utils, Model, Approach, Data, Fitter

from constants import Mp, Mp2

m = Model.ComptonGepard(ansatz='FITBP', q02=2.5)
t = Approach.hotfixedBMK(m)

# Setting gepard to values as in radNNLONS.F
# (Fig. 12 of NPB08)

t.m.g.parint.nf = 4
t.m.g.astrong.asp = np.array([0.05, 0.05, 0.05])
t.m.g.mbcont.phi = 1.9
t.m.g.parchr.fftype = np.array([c for c in 'SINGLET   ']) # array(10)
#t.m.g.parchr.fftype = np.array([c for c in 'NONSINGLET']) # array(10)

# Seting model parameters to be as in test.F
def setpar(i, val):
    t.m.parameters[t.m.parameters_index[i]] = val

MassP = 0.938272 # proton mass
setpar(11, 0.2)
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

Q2 = 2.5
xi = 1.e-5

def test_radLO():
    """Singlet LO CFF H"""
    pt.Q2 = Q2
    pt.t = -0.25
    pt.xi = xi
    pt.xB = 2*pt.xi/(1.+pt.xi)
    #FIXME: how to avoid this:
    try:
        del t.m.qdict[pt.Q2]
    except:
        pass
    t.m.g.parint.p = 0
    t.m.g.newcall = 1
    t.m.g.init()
    #print (m.ReH(pt), m.ImH(pt))
    #print np.sqrt(m.ImH(pt)**2 + m.ReH(pt)**2)
    assert_almost_equal(m.ReH(pt), 39544.823112887607)
    assert_almost_equal(m.ImH(pt), 402367.23596533033)
   
test_radLO.gepardsuite = 1
test_radLO.newfeature = 0

def test_radNLO():
    """Singlet NLO CFF H"""
    pt.Q2 = Q2
    pt.t = -0.25
    pt.xi = xi
    pt.xB = 2*pt.xi/(1.+pt.xi)
    #FIXME: how to avoid this:
    try:
        del t.m.qdict[pt.Q2]
    except:
        pass
    t.m.g.parint.p = 1
    t.m.g.newcall = 1
    t.m.g.init()
    #print (m.ReH(pt), m.ImH(pt))
    #print np.sqrt(m.ImH(pt)**2 + m.ReH(pt)**2)
    assert_almost_equal(m.ReH(pt), 5747.0424614455933)
    assert_almost_equal(m.ImH(pt), 201256.45352582674)

test_radNLO.gepardsuite = 1
test_radNLO.newfeature = 0

def test_radNNLO():
    """Singlet NNLO CFF H"""
    pt.Q2 = Q2
    pt.t = -0.25
    pt.xi = xi
    pt.xB = 2*pt.xi/(1.+pt.xi)
    try:
        del t.m.qdict[pt.Q2]
    except:
        pass
    t.m.g.parint.p = 2
    t.m.g.newcall = 1
    t.m.g.init()
    #print (m.ReH(pt), m.ImH(pt))
    #print np.sqrt(m.ImH(pt)**2 + m.ReH(pt)**2)
    assert_almost_equal(m.ReH(pt), 13475.876522596294)
    assert_almost_equal(m.ImH(pt), 197331.78427187083)

test_radNNLO.gepardsuite = 1
test_radNNLO.newfeature = 0
