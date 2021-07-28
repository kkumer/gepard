"""Testing code for anomalous dimensions and Wilson coefficients."""

import gepard as g
import numpy as np
from pytest import approx

n_test = 1.11 + 1.7j

#  -- Testing anomalous dimensions --

# Numbers are from MMA adacf.m

def test_wgamma_NS():
    """Test NS(+/-) an. dim. at LO and NLO."""
    assert g.adim.non_singlet(n_test, 3, 1) == approx(
            np.array([3.66314 +4.77559j, 23.0105 +18.9183j]), 6)
    assert g.adim.non_singlet(n_test, 3, -1) == approx(
            np.array([3.66314 +4.77559j, 23.0185 +18.8685j]), 6)


# def test_wgamma_singlet():
    # """Test singlet an. dim. at LO and NLO."""
    # assert g.adim.singlet(n_test, 3, 1) == approx(np.array([
            # [[3.66314+4.77559j, -1.13792+1.3199j],
             # [0.467652+1.54209j, 10.4322 +12.8949j]],
            # [[24.2238 +17.5493j,  4.38729 +9.94242j],
            # [-2.90091+14.6529j, 29.6072 +53.1367j]]]), 6)
