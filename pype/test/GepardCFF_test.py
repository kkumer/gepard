
import copy

import Approach
import Data
import Model
import numpy as np
import utils
from nose.tools import *

data = utils.loaddata('data/ep2epgamma', approach=Approach.hotfixedBMK)  

m = Model.ComptonGepard(scheme='CSBAR')
t = Approach.hotfixedBMK(m)

t.m.g.parint.pid = 1

# Seting model parameters to be as in test.F
def setpar(i, val):
    t.m.parameters[t.m.parameters_index[i]] = val

# nloLOParameters

setpar(11, 0.152039)
setpar(12, 1.15751)
setpar(13, 0.15)
setpar(14, 0.478391)
setpar(17, -0.15152)
setpar(18, 0.7)
setpar(19, 0.0)

#setpar(21, 0.4)
setpar(22, 1.24732)
setpar(23, 0.15)
setpar(24, 0.7)
setpar(27, -0.81217)
setpar(28, -0.2)
setpar(29, 0.0)

#setpar(111, 0.0)



pt = Data.DummyPoint()


def test_CFF():
    """Calculate CFF H."""
    pt.Q2 = 4.
    pt.t = -0.2
    pt.xi = 0.01
    t.m.g.parint.p = 0
    t.m.g.init()
    print((m.ReH(pt), m.ImH(pt)))
    # assert_almost_equal(m.ImH(pt)/100, 66.5255/100, 4)
    assert_almost_equal(m.ImH(pt)/100, 60.1622/100, 4)   # for Q2=4.
    # assert_almost_equal(m.ReH(pt)/100, 26.9984/100, 4)
    assert_almost_equal(m.ReH(pt)/100, 18.9540/100, 4)  # for Q2=4.

test_CFF.gepardsuite = 1

def test_CFFE():
    """Calculate CFF E."""
    pt.Q2 = 8.
    pt.t = -0.2
    pt.xi = 0.01
    t.m.g.parint.p = 0
    t.m.g.init()
    assert_almost_equal(m.ImE(pt)/100, 43.6024/100, 4)
    # assert_almost_equal(m.ImE(pt)/100, 42.11355/100, 4)  # for Q2=4.
    assert_almost_equal(m.ReE(pt)/100, 13.5297/100, 4)
    # assert_almost_equal(m.ReE(pt)/100, 13.2678/100, 4)  # for Q2=4.

test_CFF.gepardsuite = 1
