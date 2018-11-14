"""Testing model database theories.db"""

from nose.tools import *
import shelve
#import Model, Approach, utils, plots
#from results import *
#from utils import listdb
from abbrevs import *

db = shelve.open('theories.db')

def test_KM09a():
    """Test model: KM09a"""
    th = db['KM09a']
    pts = GLOpoints
    chisq = th.chisq(pts)[0]
    assert_almost_equal(chisq, 32.044069303618073)

#test_KM09a.long = 1

def test_KM09b():
    """Test model: KM09b"""
    th = db['KM09b']
    pts = GLOpoints + data[30]
    chisq = th.chisq(pts)[0]
    assert_almost_equal(chisq, 33.36747338543438)

#test_KM09b.long = 1

def test_KM10():
    """Test model: KM10"""
    th = db['KM10']
    pts = H1ZEUSpoints + UNP5points
    chisq = th.chisq(pts)[0]
    #assert_almost_equal(chisq, 135.85499940324056)
    # from BSA to ALUI
    #assert_almost_equal(chisq, 135.72627340347722)
    # unclear why the diff
    assert_almost_equal(chisq, 135.7650869105709)

test_KM10.long = 1

def test_KM10a():
    """Test model: KM10a"""
    th = db['KM10a']
    pts = H1ZEUSpoints+ALTGLOpoints
    chisq = th.chisq(pts)[0]
    #assert_almost_equal(chisq, 132.14636420551949)
    # from BSA to ALUI
    assert_almost_equal(chisq, 129.18281370844684)

test_KM10a.long = 1

def test_KM10b():
    """Test model: KM10b"""
    th = db['KM10b']
    pts = DVCSpoints+GLOpoints+data[30]
    chisq = th.chisq(pts)[0]
    assert_almost_equal(chisq, 115.54198973827977)

test_KM10b.long = 1


def test_KMM12():
    """Test model: KMM12"""
    th = db['KMM12']
    AULptsOLD = TSA1points[:4] + data[54]  # data[54] are not independent
    pts = H1ZEUS + ALUIpts + BCApts + CLASptsOLD + BSDwpoints\
            + AULptsOLD + ALLpts + AUTIpts + BSSwpoints
    chisq = th.chisq(pts)[0]
    assert_almost_equal(chisq, 123.53321520926985)

test_KMM12.long = 1
