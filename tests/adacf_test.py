"""Testing code for anomalous dimensions and Wilson coefficients."""

import gepard as g
import numpy as np
from pytest import approx

n_test = 1.11 + 1.7j

#  -- Testing anomalous dimensions --

# Numbers are from MMA adacf.m


def test_adim_block():
    """Test block an. dim. at LO and NLO."""
    res = g.adim.block(n_test, 3)
    assert res.shape == (1, 2, 4, 4)
    assert res[0, :, :, :] == approx(np.array([
            # LO
            [[3.66314+4.77559j, -1.13792+1.3199j, 0, 0],
             [0.467652+1.54209j, 10.4322+12.8949j, 0, 0],
             [0, 0, 3.66314+4.77559j, 0],
             [0, 0, 0, 3.66314+4.77559j]],
           # NLO
            [[24.2238+17.5493j,  4.38729+9.94242j, 0, 0],
             [-2.90091+14.6529j, 29.6072+53.1367j, 0, 0],
             [0, 0, 23.0105+18.9183j, 0],
             [0, 0, 0, 23.0185+18.8685j]]]), rel=1.e-4)
