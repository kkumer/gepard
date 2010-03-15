"""
This is collection of tests. 
Just run 'nosetest' in the pype directory.

"""

from nose.tools import *
import numpy as np

import utils, models, Approach
from fits import DMGLO1  #use some testpars here?

ff = models.FormFactors()
b = Approach.hotfixedBMK(ff, optimization = False)

data = utils.loaddata()  #FIXME: write tests without dependence on data
pt0 = data[1][0]
pt0.prepare(b)

def test_CFF():
    assert_almost_equal(ff.ImH(pt0), 17.681527454585797)
    assert_almost_equal(ff.ReH(pt0), -2.471928504474433)

def test_Xunp():
    assert_almost_equal(b.Xunp(pt0, 1, -1,  DMGLO1, {'phi':1.}), 1.8934179005138172)
    ar = b.Xunp(pt0, 1, -1,  DMGLO1, {'phi':np.array([0.3, 1.1])})
    assert isinstance(ar, np.ndarray)
    assert_almost_equal(ar[0], 1.4431498922585635)
    assert_almost_equal(ar[1], 1.9830803062066602)

def test_Xunp2():
    """New feature: any kinematic variable could be in vars.

    """
    assert_almost_equal(b.Xunp(pt0, 1, -1,  DMGLO1, 
        {'phi':1., 'xB':0.07}), 3.0240991086297577)

def test_Xunp3():
    """New feature: ndarray of any kinematical variable could be in vars.

    """
    ar = b.Xunp(pt0, 1, -1,  DMGLO1, {'phi':1, 'xB':np.array([0.07, 0.11])})
    assert isinstance(ar, np.ndarray)
    assert_almost_equal(ar[1], 0.68881435298587579)

def test_fit():
    """Testing set of fitting observables."""

    fitpoints = data[1] + data[8] + data[29] + data[30]  # DM's GLO1 set
    [pt.prepare(b) for pt in fitpoints]
    chisq = 0.
    for pt in fitpoints:
        chisq = chisq + (
                (getattr(b, pt.yaxis)(pt, DMGLO1) - pt.val)**2 / pt.err**2 )
    assert_almost_equal(chisq, 24.974520510810581)

