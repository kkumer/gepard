"""NLO part of DVMP hard scattering (Wilson) coefficients."""

import math
from typing import Tuple

from . import adim, constants, qcd
from .special import S1, S2, S3, SB3, deldelS2, delS2, parity, poch


def c1dvmp(m, sgntr: int, j: complex, k: int) -> Tuple[complex, complex, complex]:
    """NLO DVMP hard scattering coefficients.

    Args:
        m: Theory instance
        sgntr: signature (+1 or -1)
        j: conformal moment
        k: moment

    Returns:
        (quark, pure singlet, gluon) NLO coefficients.

    Todo:
        There is an unresolved discrepancy in the code implementing
        Eq. (4.53b) of Towards DVMP paper term. Dieter's notebook
        seems to have one different sign. See comments in the code
        for MCG1CF_new.

    """
    LRR2 = math.log(m.rr2)
    LRGPDF2 = math.log(m.rf2)
    LRDAF2 = math.log(m.rdaf2)
    b0 = qcd.beta(0, m.nf)

    ptyk = parity(k)

    gamGQ = adim.singlet_LO(j+1, m.nf, ptyk)[1, 0]

    gamGGCA = 4*S1(j+1) - 12/j/(j+3) + 4/(j+1)/(j+2)
    gamQQCF = 4*S1(k+1) - 3 - 2/poch(k+1, 2)
    gamQGNF = -(4+2*(j+1)*(j+2))/(j+1)/(j+2)/(j+3)

    # Spliced from Mathematica:

    #  ... quark part

    MCQ1CF = -23/3+(0.5*(1.+3.*(1.+j)*(2.+j)))/((
             1 + j)**2*(2.+j)**2)+(0.5*(1.+3.*(1.+k)*(2.+k)))/((1.+k)**2
             *(2.+k)**2)+0.5*(-3.-2./((1.+j)*(2.+j))+4.*S1(1.+j))*((-
             0.5*(1.+(1.+j)*(2.+j)))/((1.+j)*(2.+j))-(0.5*(1.+(1.+k)*(
             2.+k)))/((1.+k)*(2.+k))-LRGPDF2+S1(1.+j)+S1(1.+k))+0.5 * \
             ((-0.5*(1.+(1.+j)*(2.+j)))/((1.+j)*(2.+j))-(0.5*(1.+(1.+k
             )*(2.+k)))/((1.+k)*(2.+k))-LRDAF2+S1(1.+j)+S1(1.+k))*(-
             3.-2./((1.+k)*(2.+k))+4.*S1(1.+k))

    MCQ1BET0 = -5/6+0.5/((1.+j)*(2.+j))+0.5/((
               1 + k)*(2.+k))+0.5*LRR2-S1(1.+j)-S1(1.+k)

    SUMA = 0j
    SUMB = 0j

    for LI in range(0, k+1):
        SUMANDB = (0.125*(1.+2.*LI)*(S2(0.5*(1.+j))-S2(-0.5+0.5*(
           1.+j))+S2(-0.5+0.5*LI)-S2(0.5*LI)))/((0.5*(1.+j)-0.5*LI)*(2.+j+LI))
        SUMB += SUMANDB
        SUMA += parity(LI)*SUMANDB

    DELc1aGJK = (-2.*ptyk)/((1.+j)*(2.+j)*(1.+k)*(2.+k))+(0.5
       *(-1.+(1.+j)*(2.+j))*(-2.+(1.+k)*(2.+k)*(S2(0.5*(1.+k)) -
       S2(-0.5+0.5*(1.+k))))*ptyk)/((1.+j)*(2.+j))+(1.+k)*(2.
       +k)*(2.+0.5*(1.+k)*(2.+k))*((0.25*k*(2.+(1.+k)**2)*(S2(
       0.5*(1.+j))-S2(-0.5+0.5*(1.+j))+S2(-0.5+0.5*k)-S2(0.5*k
       )))/((0.5*(1.+j)-0.5*k)*(2.+j+k)*(3.+2.*k)*(4.+(1.+k)*(2.
       +k)))-(0.25*(S2(0.5*(1.+j))-S2(-0.5+0.5*(1.+j))-S2(0.5
       *(1.+k))+S2(-0.5+0.5*(1.+k))))/((0.5*(1.+j)+0.5*(-1.-k))
       *(3.+j+k))+(0.25*(3.+k)*(2.+(2.+k)**2)*(S2(0.5*(1.+j)) -
       S2(-0.5+0.5*(1.+j))-S2(0.5*(2.+k))+S2(-0.5+0.5*(2.+k)))
       )/((0.5*(1.+j)+0.5*(-2.-k))*(4.+j+k)*(3.+2.*k)*(4.+(1.+k)
       *(2.+k))))*ptyk+2.*(j-k)*(3.+j+k)*(-SUMA-constants.ZETA3+S3(1.+j
       )+(0.125*(1.+k)*(S2(0.5*(1.+j))-S2(-0.5+0.5*(1.+j))-S2
       (0.5*(1.+k))+S2(-0.5+0.5*(1.+k)))*ptyk)/((0.5*(1.+j)+0.5
       *(-1.-k))*(3.+j+k)))

    DELc1bGKJ = 1/((1.+k)*(2.+k))+0.5*(-2.-(1.+k)**2-((1.+j)*(2
       +j))/((1.+k)*(2.+k)))*(S2(0.5*(1.+j))-S2(-0.5+0.5*(1. +
       j)))-0.5*(1.+k)*(S2(0.5*(1.+k))-S2(-0.5+0.5*(1.+k)))-(0.125
       *(1.+k)*(2.+k)*(4.+(1.+k)*(2.+k))*(S2(0.5*(1.+j)) -
       S2(-0.5+0.5*(1.+j))-S2(0.5*(1.+k))+S2(-0.5+0.5*(1.+k)))) \
       / ((0.5*(1.+j)+0.5*(-1.-k))*(3.+j+k))-(0.5*(1.+k)*(2.+k)*(
          (0.25*k*(2.+(1.+k)**2)*(S2(0.5*(1.+j))-S2(-0.5+0.5*(1. +
       j))+S2(-0.5+0.5*k)-S2(0.5*k)))/((0.5*(1.+j)-0.5*k)*(2. +
       j+k))+(0.25*(3.+k)*(2.+(2.+k)**2)*(S2(0.5*(1.+j))-S2(-0.5
       +0.5*(1.+j))-S2(0.5*(2.+k))+S2(-0.5+0.5*(2.+k))))/((0.5
       *(1.+j)+0.5*(-2.-k))*(4.+j+k))))/(3.+2.*k)+2.*(-j+k)*(3
       +j+k)*(-SUMB-0.5*S1(1.+k)*(S2(0.5*(1.+j))-S2(-0.5+0.5
       *(1.+j)))+SB3(1+j))

    MCQ1CG = 0.9565348003631189+DELc1aGJK-(2.*(1.+(1.+j)*(2.+j)
       )*(1.-sgntr))/((1.+j)**2*(2.+j)**2)-DELc1bGKJ*sgntr+(-(1 /
       ((1.+k)*(2.+k)))+2.*S1(1.+k))*(1.-sgntr+0.5*(1.+j)*(2.+j
       )*sgntr*(-S2(0.5*j)+S2(0.5*(1.+j)))) + (2.*sgntr*ptyk)/(
               (1.+j)*(2.+j)*(1.+k)*(2.+k))-(
               2.*(1.+(1.+k)*(2.+k))*(1.+
       ptyk))/((1.+k)**2*(2.+k)**2)+(-(1/((1.+j)*(2.+j))) +
       2*S1(1.+j))*(1.+ptyk-0.5*(1.+k)*(2.+k)*(-S2(0.5*k) +
       S2(0.5*(1.+k)))*ptyk)+2.*(1.+j)*(2.+j)*((-0.5*(-1.+(1.+j)*(
       2.+j)))/((1.+j)**2*(2.+j)**2)+constants.ZETA3-(0.5*sgntr*(-S2(0.5 *
       j)+S2(0.5*(1.+j))))/((1.+j)*(2.+j))-S3(1.+j)-sgntr*SB3(
       1+j))+2.*(1.+k)*(2.+k)*((-0.5*(-1.+(1.+k)*(2.+k)))/((1. +
       k)**2*(2.+k)**2)+constants.ZETA3-S3(1.+k)+(0.5*(-S2(0.5*k) +
           S2(0.5*(1.+k)))*ptyk)/((1.+k)*(2.+k))+ptyk*SB3(1+k))

    MCQ1 = constants.CF * MCQ1CF + constants.CG * MCQ1CG + b0 * MCQ1BET0


    #  ... pure singlet quark part starting at NLO

    # MCPS1 = (-2.*(0.5+1/((1.+j)*(2.+j))+1/((1.+k)*(2.+k))))/((1
       # +j)*(2.+j))-(2.*(2.+(1.+j)*(2.+j))*(-1.-LRGPDF2+2.*S1(1
       # +j)+2.*S1(1.+k)))/(j*(1.+j)*(2.+j)*(3.+j))+(0.5*k*(1.+k
       # )*(2.+k)*(3.+k)*((0.25*(S2(0.5*(1.+j))-S2(-0.5+0.5*(1. +
       # j))+S2(-0.5+0.5*k)-S2(0.5*k)))/((0.5*(1.+j)-0.5*k)*(2. +
       # j+k))-(0.25*(S2(0.5*(1.+j))-S2(-0.5+0.5*(1.+j))-S2(0.5
       # *(2.+k))+S2(-0.5+0.5*(2.+k))))/((0.5*(1.+j)+0.5*(-2.-k))
       # *(4.+j+k))))/(3.+2.*k)

    # The commented-out code above may be slightly faster.
    # Here correction from 1612.01937 is implemented
    MCPS1 = (-LRGPDF2-m.corr_c1dvmp_one+2*S1(j+1)+2*S1(k+1)-1)*gamGQ/(j+
             3)/constants.CF - (1/2 + 1/poch(j+1, 2) +
             1/poch(k+1, 2))*2/poch(j+1,2) + poch(k,4)*(
             deldelS2((j+1)/2, k/2) -
             deldelS2((j+1)/2, (k+2)/2)) / 2 / (2*k+3)

    #  ... gluon part

    # MCG1CF_old = (-0.5*(2.+(1.+k)*(2.+k)))/((1.+j)*(2.+j)*(1.+k)*(2
       # +k))+S1(1.+j)/((1.+j)*(2.+j))-(-1.5-3./((1.+j)*(2.+j)) +
       # 0.5/((1.+k)*(2.+k))+S1(1.+j))/((1.+k)*(2.+k))-(0.5*(2.+(
       # 1.+j)*(2.+j))*(-1.5-1/((1.+j)*(2.+j))+2./((1.+k)*(2.+k)) -
       # LRGPDF2+3.*S1(1.+j)))/((1.+j)*(2.+j))+0.5*(-0.75-1/((1. +
       # j)*(2.+j))-0.5/((1.+k)*(2.+k))-LRDAF2+S1(1.+j)+S1(1.+k)
       # )*(-3.-2./((1.+k)*(2.+k))+4.*S1(1.+k))+0.125*(-39.+(2.+(
       # 1.+k)*(2.+k))*(-S2(0.5*k)+S2(0.5*(1.+k))))+0.25*(1.+k) * \
       # (2.+k)*(2.+(1.+k)*(2.+k))*((0.5*(-S2(0.5*j)+S2(0.5*(1. +
       # j))))/((1.+k)*(2.+k))-(0.5*((0.25*(-1.+k)*k*(S2(0.5*(1. +
       # j))-S2(-0.5+0.5*(1.+j))+S2(-0.5+0.5*k)-S2(0.5*k)))/((0.5
       # *(1.+j)-0.5*k)*(2.+j+k))-(0.25*(3.+k)*(4.+k)*(S2(0.5*(
       # 1.+j))-S2(-0.5+0.5*(1.+j))-S2(0.5*(2.+k))+S2(-0.5+0.5 *
       # (2.+k))))/((0.5*(1.+j)+0.5*(-2.-k))*(4.+j+k))))/(3.+2.*k)
       # )
    
    DELC1FG = ((delS2((j+1)/2)/(2*poch(k+1,2))-(poch(k-1,2)*\
             deldelS2((j+1)/2,k/2)-poch(k+3,2)*deldelS2((j+1)/2,
             (k+2)/2))/(2*(2*k+3)))*poch(k+1,2)*(poch(k+1,2)+2)/4
             -(poch(k+1,2)+2)/(2*poch(j+1,2)*poch(k+1,2)))

    # NOTE: There is an unresolved discrepancy in the next term! In the
    # Eq. (4.53b) of Towards DVMP paper term -4/poch(k+1,2)^2
    # in the second row, should have the opposite sign to agree
    # with the expression in Dieter's notebook!
    # Here, in _new, correction from 1612.01937 is implemented

    MCG1CF_new = ((-LRDAF2+S1(j+1)+S1(k+1)-3/4-
            1/(2*poch(k+1,2))-1/poch(j+1,2))*gamQQCF/2+
            (-LRGPDF2+m.corr_c1dvmp_one+3*S1(j+1)-1/2+(2*S1(j+1)-
             1)/poch(k+1,2)-1/poch(j+1,2))*(j+3)/2*gamQGNF/2
            -(35-(poch(k+1,2)+2)*delS2((k+1)/2) -   # this sign ?!
             m.corr_c1dvmp_sgn*4/poch(k+1,2)**2)/8+((poch(k+1,2)+2)*S1(j+
                 1)/poch(k+1,2)+1)/poch(j+1,2) + DELC1FG)

    # From DM's notebook:
#     MCG1CF_DM = ((-LRDAF2+S1(j+1)+S1(k+1)-3/4-
#             1/(2*poch(k+1,2))-1/poch(j+1,2))*gamQQCF/2+
#             (-LRGPDF2+3*S1(j+1)-3/2+2/poch(k+1,2)-
#                 1/poch(j+1,2))*(j+3)/2*gamQGNF/2
#             -(39-(poch(k+1,2)+2)*delS2((k+1)/2))/8-
#             (S1(j+1)-3/2-3/poch(j+1,2) +  # this sign ?!
#             1/(2*poch(k+1,2)))/poch(k+1,2) + S1(j+1)/poch(j+1,2) + DELC1FG)

    MCG1CF = MCG1CF_new

    # MCG1CA = 1.572467033424113-(4.+10.*(1.+j)*(2.+j))/((1.+j) **
       # 2*(2.+j)**2)-(3.*(-6.+2.*S1(1.+j)+S1(1.+k)))/(j*(3.+j)) \
       # +0.5*(4./((1.+j)*(2.+j))-12./(j*(3.+j))+4.*S1(1.+j))*((0.5
       # *(2.+(1.+j)*(2.+j)))/((1.+j)*(2.+j))-LRGPDF2+S1(1.+j) +
       # 1.5*S1(1.+k))+(-2.+(1.+k)*(2.+k)*S1(1.+k))/((1.+j)*(2. +
       # j)*(1.+k)*(2.+k))+0.5*(S2(0.5*j)-S2(0.5*(1.+j)))+0.125 * \
       # (2.-(1.+k)*(2.+k)*(-S2(0.5*k)+S2(0.5*(1.+k))))+0.25*k*(
       # 1.+k)*(2.+k)*(3.+k)*((-0.5*(-S2(0.5*j)+S2(0.5*(1.+j))))
       # /((1.+k)*(2.+k))-((0.25*(-1.+k)*(S2(0.5*(1.+j))-S2(-0.5
       # +0.5*(1.+j))+S2(-0.5+0.5*k)-S2(0.5*k)))/((0.5*(1.+j)-0.5
       # *k)*(2.+j+k))+(0.25*(4.+k)*(S2(0.5*(1.+j))-S2(-0.5+0.5
       # *(1.+j))-S2(0.5*(2.+k))+S2(-0.5+0.5*(2.+k))))/((0.5*(1.
       # +j)+0.5*(-2.-k))*(4.+j+k)))/(3.+2.*k))

    # Eq. (4.53e)
    DELC1AG = (-(delS2((j+1)/2)/2/poch(k+1,2) + ((k-1)*deldelS2((j+1)/2,k/2) +
              (k+4)*deldelS2((j+1)/2,(k+2)/2))/(2*k+3))*poch(k,4)/4 +
           (poch(k+1,2)*S1(k+1)-2)/(poch(j+1,2)*poch(k+1,2)))

    # Eq. (4.53a)
    MCG1CA = ((-LRGPDF2+S1(j+1)+3*S1(k+1)/2+1/2 +
              1/poch(j+1,2))*gamGGCA/2 -
              3*(2*S1(j+1)+S1(k+1)-6)/(j*(j+3)) +
              (8+4*constants.ZETA2-poch(k+1,2)*delS2((k+1)/2))/8 -
              delS2((j+1)/2)/2-(10*poch(j+1,2)+4)/(poch(j+1,2)**2) +
              DELC1AG)

    MCG1BET0 = -0.5*LRGPDF2+0.5*LRR2

    MCG1 = constants.CF * MCG1CF + constants.CA * MCG1CA + b0 * MCG1BET0

    return (MCQ1, MCPS1, MCG1)
