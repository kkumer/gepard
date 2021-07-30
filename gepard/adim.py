"""QCD anomalous dimensions.

Note:
    Functions in this module have as first argument Mellin moment
    n=j+1, where j is conformal moment used everywhere else.
    Thus, they should always be called as f(j+1, ...).
"""
# from cmath import exp
# from scipy.special import loggamma as clngamma

import numpy as np

from gepard.constants import CA, CF, CG, TF
from gepard.special import (S1, S2, S3, MellinF2, S2_prime, S2_tilde, S3_prime,
                            poch, psi, zeta)


def non_singlet_LO(n: complex, nf: int, prty: int=1) -> complex:
    """Non-singlet LO anomalous dimension.

    Args:
        n (complex): which moment (= Mellin moment for integer n)
        nf (int): number of active quark flavors
        prty (int): 1 for NS^{+}, -1 for NS^{-}, irrelevant at LO

    """
    return CF*(-3.0-2.0/(n*(1.0+n))+4.0*S1(n))


def singlet_LO(n: complex, nf: int, prty: int=1) -> np.ndarray:
    """Singlet LO anomalous dimensions.

    Args:
        n (complex): which moment (= Mellin moment for integer n)
        nf (int): number of active quark flavors
        prty (int): C parity, irrelevant at LO

    Returns:
        2x2 complex matrix ((QQ, QG),
                            (GQ, GG))

    """
    qq0 = CF*(-3.0-2.0/(n*(1.0+n))+4.0*S1(n))
    qg0 = (-4.0*nf*TF*(2.0+n+n*n))/(n*(1.0+n)*(2.0+n))
    gq0 = (-2.0*CF*(2.0+n+n*n))/((-1.0+n)*n*(1.0+n))
    gg0 = (-22*CA/3.-8.0*CA*(1/((-1.0+n)*n)+1/((1.0+n)*(2.0+n))-S1(n))+8*nf*TF/3.)/2.

    return np.array([[qq0, qg0], 
                     [gq0, gg0]])



def non_singlet_NLO(n: complex, nf: int, prty: int) -> complex:
    """Non-singlet anomalous dimension.

    Args:
        n (complex): which moment (= Mellin moment for integer n)
        nf (int): number of active quark flavors
        prty (int): 1 for NS^{+}, -1 for NS^{-}

    """
    # From Curci et al.
    nlo = (CF * CG * (
            16*S1(n)*(2*n+1)/poch(n, 2)**2 +
            16*(2*S1(n) - 1/poch(n, 2)) * (S2(n)-S2_prime(n/2, prty)) +
            64 * S2_tilde(n, prty) + 24*S2(n) - 3 - 8*S3_prime(n/2, prty) -
            8*(3*n**3 + n**2 - 1)/poch(n, 2)**3 -
            16*prty*(2*n**2 + 2*n + 1)/poch(n, 2)**3) +
           CF * CA * (S1(n)*(536/9 + 8*(2*n+1)/poch(n, 2)**2) - 16*S1(n)*S2(n) +
            S2(n)*(-52/3 + 8/poch(n, 2)) - 43/6 -
            4*(151*n**4 + 263*n**3 + 97*n**2 + 3*n + 9)/9/poch(n, 2)**3) +
           CF * nf * TF * (-(160/9)*S1(n) + (32/3)*S2(n) + 4/3 +
               16*(11*n**2 + 5*n - 3)/9/poch(n, 2)**2)) / 4

    return nlo


def singlet_NLO(n: complex, nf: int, prty: int=1) -> np.ndarray:
    """Singlet NLO anomalous dimensions matrix.

    Args:
        n (complex): which moment (= Mellin moment for integer n)
        nf (int): number of active quark flavors
        prty (int): C parity

    Returns:
        Matrix (LO, NLO) where each is in turn
        2x2 complex matrix ((QQ, QG),
                            (GQ, GG))

    """
    qq1 = non_singlet_NLO(n, nf, 1) - 4*CF*TF*nf*(5*n**5+32*n**4+49*n**3+38*n**2 +
            28*n+8)/((n-1)*n**3*(n+1)**3*(n+2)**2)

    qg1 = 0.25*(-8*CF*nf*TF*((-4*S1(n))/n**2+(4+8*n +26*n**3+11*n**4 +
              15*(n*n))/(n**3*(1+n)**3*(2+n))+((2+n+n*n)*(5-2*S2(n) +
               2*(S1(n)*S1(n))))/(n*(1+n)*(2+n))) -
              8*CA*nf*TF*((8*(3+2*n)*S1(n))/((1+n)**2*(2+n)**2) +
               (2*(16+64*n+128*n**3+85*n**4+36*n**5+25*n**6 +
                14*n**7+6*n**8+n**9+104*(n*n)))/(
                (-1+n)*n**3*(1+n)**3*(2+n)**3)+((2+n+n*n)*(2*S2(n) -
                2*(S1(n)*S1(n))-2*S2(n/2)))/(n*(1+n)*(2+n))))
 
 
    gq1 = 0.25*(-10.666666666666667*CF*nf*TF*((1+n)**(-2) +
              ((-2.6666666666666665+S1(n))*(2+n+n*n))/((-1+n)*n*(1+n))) -
              4*(CF*CF)*((-4*S1(n))/(1+n)**2-(-4-12*n+28*n**3+43*n**4 +
                30*n**5+12*n**6-n*n)/((-1+n)*n**3*(1+n)**3) +
                ((2+n+n*n)*(10*S1(n)-2*S2(n)-2*(S1(n)*S1(n))))/((-1+n)*n*(1+n))) -
              8*CF*CA*((0.11111111111111112*(144+432*n-1304*n**3-1031*n**4 +
                  695*n**5+1678*n**6+1400*n**7+621*n**8+109*n**9 -
                  152*(n*n)))/((-1+n)**2*n**3*(1+n)**3*(2+n)**2) -
                  (0.3333333333333333*S1(n)*(-12-22*n+17*n**4 +
                41*(n*n)))/((-1+n)**2*n**2*(1+n))+((2+n+n*n)*(S2(n) +
                 S1(n)*S1(n)-S2(n/2)))/((-1+n)*n*(1+n))))
 

    gg1 = 0.25*(CF*nf*TF*(8+(16*(-4-4*n-10*n**3+n**4+4*n**5+2*n**6 -
              5*(n*n)))/((-1+n)*n**3*(1+n)**3*(2+n)))+CA*nf*TF*(10.666666666666667 -
               17.77777777777778*S1(n)+(1.7777777777777777*(12+56*n +
                 76*n**3+38*n**4+94*(n*n)))/((-1+n)*n**2*(1+n)**2*(2+n))) +
               CA*CA*(-21.333333333333333+59.55555555555556*S1(n)+(64*S1(n)*(-2 -
                2*n+8*n**3+5*n**4+2*n**5+7*(n*n)))/((-1+n)**2*n**2*(1+n)**2*(2+n)**2) -
                (0.4444444444444445*(576+1488*n-1632*n**3-2344*n**4+1567*n**5 +
                6098*n**6+6040*n**7+2742*n**8+457*n**9+560*(n*n)))/((-1+n)**2*n**3*(1 +
                n)**3*(2+n)**3)-16*S1(n)*S2(n/2)+(32*(1+n+n*n)*S2(n/2))/((-1 +
                n)*n*(1+n)*(2+n))-4*S3(n/2)+32*(S1(n)/n**2-0.625*zeta(3)+MellinF2(n) -
                    0.5*zeta(2)*(-psi(n/2)+psi((1+n)/2)))))


    return np.array([[qq1, qg1], 
                     [gq1, gg1]])


def block(n: complex, nf: int) -> np.ndarray:
    """LO+NLO block anomalous dimensions matrix.

    Args:
        n (complex): which moment (= Mellin moment for integer n)
        nf (int): number of active quark flavors

    Returns:
        2x4x4 array (LO, NLO) where each is in turn
        4x4 complex matrix ((QQ, QG,   0,   0),
                            (GQ, GG,   0,   0),
                            ( 0,  0, NS+,   0),
                            ( 0,  0,   0, NS-))
    """
    lo_SI = singlet_LO(n, nf)
    lo_NS = non_singlet_LO(n, nf)

    nlo_SI = singlet_NLO(n, nf)
    nlo_plus_NS = non_singlet_NLO(n, nf, prty=1)
    nlo_minus_NS = non_singlet_NLO(n, nf, prty=-1)

    return np.block([
                     [[lo_SI, np.zeros((2, 2))],
                      [np.zeros((2, 2)), np.diagflat([lo_NS, lo_NS])]],
                     # NLO
                     [[nlo_SI, np.zeros((2, 2))],
                      [np.zeros((2, 2)), np.diagflat([nlo_plus_NS, nlo_minus_NS])]]]) 
