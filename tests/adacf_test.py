"""Testing code for anomalous dimensions and Wilson coefficients."""

from pytest import approx

import gepard as g

n_test = 1.11 + 1.7j

#  -- Testing anomalous dimensions - LO


def test_wgammaNS():
    """Test NS+ at LO."""
    assert g.adim.NS(n_test, 0, 3, 1) == approx(
            3.663142570252238+4.775593439086258j)


