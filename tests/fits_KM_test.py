"""Testing published KM fits."""

import gepard as g
from gepard.fits import *
from pytest import approx, fixture, mark


def test_KM09a():
    """Test model: KM09a."""
    chisq = th_KM09a.chisq(GLOpoints)
    assert chisq == approx(32.044069303618073)


def test_KM09b():
    """Test model: KM09b."""
    pts = GLOpoints + g.data.dset[30]
    chisq = th_KM09b.chisq(pts)
    assert chisq == approx(33.36747338543438)


@mark.slow
def test_KM10():
    """Test model: KM10."""
    pts = H1ZEUSpoints + UNP5points
    chisq = th_KM10.chisq(pts)
    # assert_almost_equal(chisq, 135.85499940324056)
    # # from BSA to ALUI
    # assert_almost_equal(chisq, 135.72627340347722)
    # # unclear why the diff
    assert chisq == approx(135.7650869105709)


@mark.skip('KM10a not yet transferred')
def test_KM10a():
    """Test model: KM10a."""
    pts = H1ZEUSpoints+ALTGLOpoints
    chisq = th_KM10a.chisq(pts)[0]
    #assert_almost_equal(chisq, 132.14636420551949)
    # from BSA to ALUI
    assert chisq == approx(129.18281370844684)


@mark.slow
def test_KM10b():
    """Test model: KM10b."""
    pts = DVCSpoints+GLOpoints+g.data.dset[30]
    chisq = th_KM10b.chisq(pts)
    assert chisq == approx(115.54198973827977)


@mark.skip('KMM12 not yet transferred')
def test_KMM12():
    """Test model: KMM12."""
    AULptsOLD = TSA1points[:4] + data[54]  # data[54] are not independent
    pts = H1ZEUS + ALUIpts + BCApts + CLASptsOLD + BSDwpoints\
            + AULptsOLD + ALLpts + AUTIpts + BSSwpoints
    chisq = th_KMM12.chisq(pts)[0]
    # For difference to EPJA review, see comment for KM15 below
    assert chisq == approx(123.53321520926985)


@mark.slow
def test_KM15():
    """Test model: KM15."""
    pts = GLO15b
    #H1ZEUS + ALUIpts + BCApts + CLAS08pts + AULpts + ALLpts + AUTIpts + GLO15new
    chisq = th_KM15.chisq(pts)
    # Number is slightly different from EPJA review due to change in
    # treatment of uncertainties of Hall A 2015 data 
    # in commit 4033d7accb5bea7c371e4145343ef650cb38b6b9
    # on Feb 14 2017 
    assert chisq == approx(245.82778194537178)
