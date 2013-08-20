
from nose.tools import *
import numpy as np

import utils, Model, Approach
from results import DM12

data = utils.loaddata('data/gammastarp2gammap', approach=Approach.hotfixedBMK)  

m = Model.Gepard(ansatz='EFLEXP')
th = Approach.BMK(m)
m.parameters.update(DM12)

pt0 = data[39][0]
pt0.t = -0.1
pt0.Q2 = 4.
pt0.xi = 0.001
pt0.xB = 2.*pt0.xi/(1.+pt0.xi)



def test_GPDtraj():
    """Calculate GPDs on trajectory xi=x"""
    assert_almost_equal(m.gpdHtrajQ(pt0)/1000., 1.34402, 5)
    assert_almost_equal(m.gpdHtrajG(pt0), 2.6866, 4)
    assert_almost_equal(m.gpdEtrajQ(pt0)/1000., 1.98069, 5)
    assert_almost_equal(m.gpdEtrajG(pt0)*100, 3.58365, 5)

test_GPDtraj.gepardsuite = 1

def test_GPDzero():
    """Calculate GPDs on trajectory xi=0"""
    assert_almost_equal(m.gpdHzeroQ(pt0)/1000., 1.83647, 5)
    assert_almost_equal(m.gpdHzeroG(pt0), 6.9642, 4)
    assert_almost_equal(m.gpdEzeroQ(pt0)/1000., 3.09861, 5)
    assert_almost_equal(m.gpdEzeroG(pt0)*100, 7.56646, 5)

test_GPDzero.gepardsuite = 1
