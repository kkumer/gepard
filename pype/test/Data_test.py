
import copy
from nose.tools import *
import numpy as np

import utils, Approach

data = utils.loaddata('data/ep2epgamma', approach=Approach.hotfixedBMK)  

pt0 = data[34][-1]

def test_errs_addition():
    """Testing addition of stat and syst errs."""
    assert_almost_equal(pt0.errstat**2+pt0.errsyst**2, 5.2602024250699046e-05)

def test_errs_HallA06():
    """Hall A 2006 errsyst=errnorm"""
    assert_almost_equal(pt0.errnorm, pt0.errsyst)

test_errs_HallA06.extendedtesting = 1
