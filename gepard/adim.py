"""QCD anomalous dimensions."""

# from cmath import exp
# from scipy.special import loggamma as clngamma

import numpy as np

from gepard.constants import CA, CF, CG, TF
from gepard.special import (S1, S2, S3, MellinF2, S2_prime, S2_tilde, S3_prime,
                            poch, psi, zeta)


def non_singlet_LO(n: complex, nf: int, prty: int) -> complex:
    """Non-singlet LO anomalous dimension.

    Args:
        n (complex): which moment (= Mellin moment for integer n)
        nf (int): number of active quark flavors
        prty (int): 1 for NS^{+}, -1 for NS^{-}

    """

    return CF*(-3.0-2.0/(n*(1.0+n))+4.0*S1(n))


def non_singlet(n: complex, nf: int, prty: int) -> complex:
    """Non-singlet anomalous dimension.

    Args:
        n (complex): which moment (= Mellin moment for integer n)
        nf (int): number of active quark flavors
        prty (int): 1 for NS^{+}, -1 for NS^{-}

    """
    lo = non_singlet_LO(n, nf, prty)

#     # Old, valid only for prty=1
#     nlo = 0.25*(CF*nf*TF*(1.3333333333333333-17.77777777777778*S1(n) +
#             10.666666666666667*S2(n)+(1.7777777777777777*(-3+5*n +
#             11*(n*n)))/(n**2*(1+n)**2))+CF*CA*(-
#           7.166666666666667+(59.55555555555556+(8*(1 +
#           2*n))/(n**2*(1+n)**2))*S1(n)+(-17.333333333333332 +
#           8/(n*(1+n)))*S2(n)-16*S1(n)*S2(n)-(0.4444444444444445*(9 +
#           3*n+263*n**3+151*n**4+97*(n*n)))/(n**3*(1+n)**3)) +
#           CF*CG*(-3+(16*(1+2*n)*S1(n))/(n**2*(1+n)**2) +
#               24*S2(n)-(8*(-1+3*n**3+n*n))/(n**3*(1+n)**3) -
#               (16*(1+2*n+2*(n*n)))/(n**3*(1+n)**3)+16*(-(1/(n*(1+n))) +
#             2*S1(n))*(S2(n)-S2(n/2))-8*S3(n/2)+64*(S1(n)/n**2 -
#             0.625*zeta(3)+MellinF2(n)-0.5*zeta(2)*(-psi(n/2)+psi((1+n)/2))))) 

#     j = n - 1
#     # MC, TF-prop part is wront
#     nlo = (CF*nf*TF/3*(1+4/3*(13+27*j+11*j**2)/poch(j+1, 2)**2 -
#             40/3*S1(j+1)+8*S2(j+1)) +
#           CF*CA/4*(-43/6-4/9*(523+1590*j+1792*j**2+867*j**3 +
#             151*j**4)/poch(j+1,2)**3 +
#             (536/9 + 8*(2*j+3)/poch(j+1, 2)**2)*S1(j+1) +
#             (8/poch(j+1, 2) - 52/3)*S2(j+1) - 16*S1(j+1)*S2(j+1)) +
#           CF*CG*(-3/4-2*(3+11*j+10*j**2+3*j**3)/poch(j+1, 2)**3 +
#           4*(-prty)*(5+6*j+2*j**2)/poch(j+1, 2)**3 + 
#           4*(2*j+3)*S1(j+1)/poch(j+1, 2)**2 +
#           16*MellinF2(j+1) + 6*S2(j+1) + 4*(2*S1(j+1)-1/poch(j+1, 2))*(S2(j+1) -
#          (1-(-prty))/2*S2((j+1)/2)-(1+(-prty))/2*S2(j/2)) -
#           2*((1-(-prty))/2*S3((j+1)/2)+(1+(-prty))/2*S3(j/2))))


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



    return np.array([lo, nlo])


def singlet_LO(n: complex, nf: int, prty: int=1) -> np.ndarray:
    """Singlet LO anomalous dimensions.

    Args:
        n (complex): which moment (= Mellin moment for integer n)
        nf (int): number of active quark flavors
        prty (int): C parity

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


def singlet(n: complex, nf: int, prty: int=1) -> np.ndarray:
    """Singlet anomalous dimensions.

    Args:
        n (complex): which moment (= Mellin moment for integer n)
        p (int): order of pQCD series (0=LO, 1=NLO, ...)
        nf (int): number of active quark flavors
        prty (int): C parity

    Returns:
        Matrix (LO, NLO) where each is in turn
        2x2 complex matrix ((QQ, QG),
                            (GQ, GG))

    """


    j = n - 1
    qq1 = non_singlet(n, nf, 1)[1] - 4*CF*nf*TF*(160+404*j+427*j**2+227*j**3 +
                             57*j**4+5*j**5)/(j*poch(j+1,2)*poch(j+1,3)**2)

#     # MC WRONG 
#     qg1 = (-2*CA*nf*TF*((-2*S1(j+1)**2+2*S2(j+1)-(1-(-prty))*S2((j+1)/2) -
#           (1+(-prty))*S2(j/2))*(j**2+3*j+4)/poch(j+1,3) +
#           8*S1(j+1)*(2*j+5)/poch(j+2,2)**2+2*(480+1488*j+2252*j**2+2273*j**3 +
#           1711*j**4+963*j**5+382*j**6+99*j**7+15*j**8+j**9)/(j*poch(j+1,3)**3)) -
#           2*CF*nf*TF*((2*S1(j+1)**2-2*S2(j+1)+5)*(j**2+3*j+4)/poch(j+1,3) -
#          4*S1(j+1)**2/(j+1)**2+(64+160*j+159*j**2+70*j**3+11*j**4)/poch(j+1,3)**3 ))
#         
#     # MC OK
#     gq1 = (-CF**2*((-2*S1(j+1)**2+10*S1(j+1)-2*S2(j+1))*(4+3*j+j**2)/(j*(j+1)*(j+2)) -
#           4*S1(j+1)/(j+2)**2-(96+464*j+821*j**2+740*j**3+373*j**4+102*j**5 +
#           12*j**6)/(j*(j+1)**3*(j+2)**3)) -
#           2*CF*CA*((S1(j+1)**2+S2(j+1)-(1-(-prty))/2*S2((j+1)/2) -
#             (1+(-prty))/2*S2(j/2))*(4+3*j+j**2)/(j*(j+1)*(j+2))-S1(j+1)*(24 +
#             128*j+143*j**2+68*j**3+17*j**4)/(3*j**2*(j+1)**2*(j+2)) +
#             (2592+21384*j+72582*j**2+128014*j**3+133818*j**4+88673*j**5 +
#             38022*j**6+10292*j**7+1602*j**8+109*j**9)/(9*j**2*(j+1)**3*(j +
#             2)**3*(j+3)**2))-8/3*CF*nf*TF*((S1(j+1) -
#             8/3)*(4+3*j+j**2)/(j*(j+1)*(j+2))+1/(j+2)**2))
# 
#     # MC WRONG 
#     gg1 = (CA*nf*TF*(-40/9*S1(j+1)+8/3+8/9*(138+312*j+275*j**2+114*j**3 +
#            19*j**4)/(j*(j+1)**2*(j+2)**2*(j+3)))+CF*nf*TF*(2+4*(-16-8*j +
#             41*j**2+74*j**3+51*j**4+16*j**5+2*j**6)/(j*(j+1)**3*(j+2)**3*(j+3))) +
#            CA**2*(134/9*S1(j+1)+16*S1(j+1)*(18+66*j+81*j**2+48*j**3+15*j**4 +
#             2*j**5)/(j**2*(j+1)**2*(j+2)**2*(j+3)**2) -
#             16/3+8*((1-(-prty))/2*S2((j+1)/2)+(1+(-prty))/2*S2(j/2))*(3 +
#             3*j+j**2)/(j*(j+1)*(j+2)*(j+3))-4*S1(j+1)*((1-(-prty))/2*S2((j+1)/2) +
#             (1+(-prty))/2*S2(j/2))+8*MellinF2(j+1)-(1-(-prty))/2*S3((j+1)/2) -
#             (1+(-prty))/2*S3(j/2)-1/9*(15552 + 101088*j + 308808*j**2 + 
#             529962*j**3 + 557883*j**4 + 376129*j**5 + 163542*j**6 + 
#             44428*j**7 + 6855*j**8 + 457*j**9)/(j**2*(j+1)**3*(j+2)**3*(j+3)**3)))
# 
 
    qg1 = 0.25*(-8*CF*nf*TF*((-4*S1(n))/n**2+(4+8*n +26*n**3+11*n**4 +
              15*(n*n))/(n**3*(1+n)**3*(2+n))+((2+n+n*n)*(5-2*S2(n) +
               2*(S1(n)*S1(n))))/(n*(1+n)*(2+n))) -
              8*CA*nf*TF*((8*(3+2*n)*S1(n))/((1+n)**2*(2+n)**2) +
               (2*(16+64*n+128*n**3+85*n**4+36*n**5+25*n**6 +
                15*n**7+6*n**8+n**9+104*(n*n)))/(
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


    return np.array([
                     singlet_LO(n, nf, prty),
                     [[qq1, qg1], [gq1, gg1]]
                    ])
