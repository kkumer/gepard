"""QCD anomalous dimensions."""

# from cmath import exp
# from scipy.special import loggamma as clngamma

from gepard.constants import CF
from gepard.harmonic import S1

def NS(n: complex, p: int, nf: int, prty: int) -> complex:
    """Non-singlet anomalous dimension.

    Args:
        n (complex): index of moment (= Mellin moment for integer n)
        p (int): order of perturbatio series (0=LO, 1=NLO, ...)
        nf (int): number of active quark flavors
        prty (int): 1 for NS^{+}, -1 for NS^{-}

    """
    # Make sure we only ask for what's implemented
    assert prty == 1
    assert p == 0

    lo = CF*(-3.0-2.0/(n*(1.0+n))+4.0*S1(n))

    return lo
