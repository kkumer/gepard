"""QCD anomalous dimensions."""

# from cmath import exp
# from scipy.special import loggamma as clngamma

import numpy as np

from gepard.constants import CA, CF, TF
from gepard.harmonic import S1


def non_singlet(n: complex, p: int, nf: int, prty: int) -> complex:
    """Non-singlet anomalous dimension.

    Args:
        n (complex): which moment (= Mellin moment for integer n)
        p (int): order of pQCD series (0=LO, 1=NLO, ...)
        nf (int): number of active quark flavors
        prty (int): 1 for NS^{+}, -1 for NS^{-}

    """
    # Make sure we only ask for what's implemented
    assert prty == 1
    assert p == 0

    lo = CF*(-3.0-2.0/(n*(1.0+n))+4.0*S1(n))

    return lo


def singlet(n: complex, p: int, nf: int, prty: int) -> np.ndarray:
    """Singlet anomalous dimensions.

    Args:
        n (complex): which moment (= Mellin moment for integer n)
        p (int): order of pQCD series (0=LO, 1=NLO, ...)
        nf (int): number of active quark flavors
        prty (int): 1 for NS^{+}, -1 for NS^{-}

    Returns:
        2x2 complex matrix ((QQ, QG),
                            (GQ, GG))

    """
    # Make sure we only ask for what's implemented
    assert prty == 1
    assert p == 0

    qq = CF*(-3.0-2.0/(n*(1.0+n))+4.0*S1(n))
    qg = (-4.0*nf*TF*(2.0+n+n*n))/(n*(1.0+n)*(2.0+n))
    gq = (-2.0*CF*(2.0+n+n*n))/((-1.0+n)*n*(1.0+n))
    gg = (-22*CA/3.-8.0*CA*(1/((-1.0+n)*n)+1/((1.0+n)*(2.0+n))-S1(n))+8*nf*TF/3.)/2.

    return np.array([[qq, qg], [gq, gg]])
