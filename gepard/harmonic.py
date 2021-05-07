"""Harmonic sums and related functions."""

# from cmath import exp
from scipy.special import psi

from gepard.constants import EMC

def S1(z: complex) -> complex:
    """Harmonic sum S_1."""

    S1 = EMC + psi(z+1)

    return S1
