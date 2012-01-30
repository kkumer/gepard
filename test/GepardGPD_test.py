
from nose.tools import *
import numpy as np

import utils, Model, Approach
from results import KM10

data = utils.loaddata('data/gammastarp2gammap', approach=Approach.hotfixedBMK)  

m = Model.Gepard(ansatz='FIT')
m.parameters.update(KM10)
m.parameters['KAPS'] = 1.0
#m.parameters['KAPS'] = 0
th = Approach.BMK(m)

pt0 = data[39][0]


def test_GPDtraj():
    """Calculate GPDs on trajectory xi=x"""
    assert_almost_equal(m.gpdHtrajQ(pt0)/1000., 4.222761776)
    assert_almost_equal(m.gpdHtrajG(pt0), -1.266816602)
    assert_almost_equal(m.gpdEtrajQ(pt0)/1000., 4.222761776)
    assert_almost_equal(m.gpdEtrajG(pt0), 0)
