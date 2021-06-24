"""Testing code for special functions of complex argument."""

import gepard as g
import numpy as np
from pytest import approx

z_test = 1.7 + 3.7j
zar_test = np.array([z_test, z_test + 10, z_test + 10j])

#  -- Numbers below are produced by Wolfram Mathematica


def test_pochhammer():
    """Test Pochammer symbol."""
    assert g.special.pochhammer(z_test, 3) == approx(
            -93.906 + 26.566j, rel=1.e-14)


def test_S1():
    """Test Harmonic sum."""
    assert g.special.S1(z_test) == approx(
            2.03584940426394444 + 1.03236755885145415j, rel=1.e-14)


def test_S2():
    """Test Harmonic sum."""
    assert g.special.S2(z_test) == approx(
            1.525155519155933 + 0.199617618107181j, rel=1.e-14)


def test_S2_array():
    """Test Harmonic sum of numpy array."""
    assert g.special.S2(zar_test) == approx(
          np.array([1.525155519155933 + 0.199617618107181j,
                    1.569896512574346 + 0.0227340297888166j,
                    1.633492895143775 + 0.071185483777968j]), rel=1.e-14)


def test_S3():
    """Test Harmonic sum."""
    assert g.special.S3(z_test) == approx(
             1.2147345977685899+0.0240107948528736j, rel=1.e-14)


def test_S4():
    """Test Harmonic sum."""
    assert g.special.S4(z_test) == approx(
            1.08654557333671745 + 0.00026621041228996j, rel=1.e-14)


def test_SB3():
    """Test function from Eq. (4.44e) of arXiv:1310.5394."""
    # The number is from fortran-gepard, not Mathematica
    assert g.special.SB3(z_test) == approx(
           0.00106358734374320529-0.00111600322992239098j, rel=1.e-5)
