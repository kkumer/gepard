"""
Testing new flexible E ansaetze EFL(EXP) and EPH(EXP)
"""

import copy
from nose.tools import *
import numpy as np

import utils, Model, Approach, Data
from math import exp

m = Model.ComptonGepard(ansatz='EPH')
t = Approach.hotfixedBMK(m)

# 11 : 'NS',        112 : 'EAL0S',
# 12 : 'AL0S',      113 : 'EALPS',
# 13 : 'ALPS',      114 : 'EM02S',
# 14 : 'M02S',      115 : 'EDELM2S',
# 15 : 'DELM2S',    116 : 'EPS',
# 16 : 'PS',        117 : 'ESECS',
# 17 : 'SECS',      119 : 'ESKEWS',
# 18 : 'KAPS',      122 : 'EAL0G',
# 19 : 'SKEWS',     123 : 'EALPG',
# 21 : 'NG',        124 : 'EM02G',
# 22 : 'AL0G',      125 : 'EDELM2G',
# 23 : 'ALPG',      126 : 'EPG',
# 24 : 'M02G',      127 : 'ESECG',
# 25 : 'DELM2G',    129 : 'ESKEWG',
# 26 : 'PG',        137 : 'ETHIS',
# 27 : 'SECG',      147 : 'ETHIG' }
# 28 : 'KAPG',
# 29 : 'SKEWG',
# 37 : 'THIS',
# 47 : 'THIG',
# 48 : 'DELB',

# nloLOParameters
m.parameters['NS']     =  0.152039
m.parameters['AL0S']   =  1.15751
m.parameters['M02S']   =  0.478391
m.parameters['SECS']   = -0.15152
m.parameters['THIS']   =  0.0
m.parameters['KAPS']   =  0.7
m.parameters['SKEWS']  =  0.0

m.parameters['AL0G']   =  1.24732
m.parameters['M02G']   =  0.7
m.parameters['SECG']   = -0.81217
m.parameters['THIG']   =  0.0
m.parameters['KAPG']   = -0.2
m.parameters['SKEWG']  =  0.0

m.parameters['DELB']   =  1.0

m.parameters['EAL0S']   =  1.15751
m.parameters['EM02S']   =  0.478391
m.parameters['ESECS']   = -0.15152
m.parameters['ETHIS']   =  0.0
m.parameters['EKAPS']   =  0.7
m.parameters['ESKEWS']  =  0.0

m.parameters['EAL0G']   =  1.24732
m.parameters['EM02G']   =  0.7
m.parameters['ESECG']   = -0.81217
m.parameters['ETHIG']   =  0.0
#m.parameters['EKAPG']   = -0.2  #ignored par
m.parameters['ESKEWG']  =  0.0

pt = Data.DummyPoint()

def test_CFF():
    """Calculate CFF H for EPH ansatz (Equal to FIT)."""
    pt.Q2 = 8.
    pt.t = -0.2
    pt.xi = 0.01
    t.m.g.parint.p = 0
    t.m.g.init()
    assert_almost_equal(m.ImH(pt)/100, 66.5255/100, 4)
    assert_almost_equal(m.ReH(pt)/100, 26.9984/100, 4)

test_CFF.gepardsuite = 1

def test_CFFE():
    """Calculate CFF E for EPH ansatz."""
    pt.Q2 = 8.
    pt.t = -0.2
    pt.xi = 0.01
    t.m.g.parint.p = 0
    t.m.g.init()
    assert_almost_equal(m.ImE(pt)/100, 35.5974571/100, 4)
    assert_almost_equal(m.ReE(pt)/100, 10.8939245/100, 4)

test_CFFE.gepardsuite = 1
