"""Testing the Goloskokov-Kroll GPD/CFF model."""

import gepard as g
from pytest import approx, fixture, mark

gpdargs = (0.21, 0.13, -0.2, 7.)

pt = g.DataPoint({'xB': 0.1, 't': -0.28, 'Q2': 2.})

@fixture
def thGK():
    th = g.cff.GoloskokovKrollCFF()
    return th


def test_GKGPDHuval(thGK):
    """Test GK model: H^u_val."""
    res = thGK.Huval(*gpdargs)
    assert res == approx(2.249834779299078)


def test_GKGPDHudsea(thGK):
    """Test GK model: H^u_sea."""
    res = thGK.Hudsea(*gpdargs)
    assert res == approx(0.136954792369421)


def test_GKImH(thGK):
    """Test GK model: Im(CFFH)."""
    res = thGK.ImH(pt)
    assert res == approx(13.1057763522748)


def test_GKReH(thGK):
    """Test GK model: Re(CFFH)."""
    res = thGK.ReH(pt)
    assert res == approx(-1.94670099930845)


def test_GKImHt(thGK):
    """Test GK model: Im(CFFHt)."""
    res = thGK.ImHt(pt)
    assert res == approx(2.3675474323)


def test_GKReHt(thGK):
    """Test GK model: Re(CFFHt)."""
    res = thGK.ReHt(pt)
    assert res == approx(1.0193853587666957)


def test_GKImE(thGK):
    """Test GK model: Im(CFFE)."""
    res = thGK.ImE(pt)
    assert res == approx(-3.3357207038)


def test_GKReE(thGK):
    """Test GK model: Re(CFFE)."""
    res = thGK.ReE(pt)
    assert res == approx(-5.1576181607)


def test_GKImEt(thGK):
    """Test GK model: Im(CFFEt)."""
    res = thGK.ImEt(pt)
    #assert res == approx(1.960464892) # with wrong al0 from KMS paper
    assert res == approx(39.522321114885344)


def test_GKReEt(thGK):
    """Test GK model: Re(CFFEt)."""
    res = thGK.ReEt(pt)
    #assert res == approx(63.030901624) # with wrong al0 from KMS paper
    assert res == approx(100.94485648382198)
