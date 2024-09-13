"""Testing code for DA models."""

import gepard as g
import numpy as np
from pytest import approx, mark

CZ = g.GegenbauerDA(p=0, ngegens=2)
CZ.parameters['a2'] = 2/3

def test_DA():
    """Evaluate all DA Gegenbauers."""
    assert 3*CZ.gegenbauers().sum() == approx(5)

def test_DA_xspace():
    """Evaluate DA at x=0.5."""
    assert CZ.x(0.5) == approx(0)

