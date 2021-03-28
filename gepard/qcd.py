"""Perturbative QCD functions.

Coefficients of beta function of QCD up to order alpha_s^5 (N^3LO),
normalized by

     d a_s / d ln mu_r^2  =   BETA0 a_s^2 + BETA1 a_s^3 + ...

with  a_s = alpha_s/(4*pi).

Beyond NLO the QCD colour factors are hard-wired in this routine,
and the numerical coefficients are truncated to six digits.

QCD coupling alpha strong is obtained by integrating the evolution
equation for a fixed number of massless flavours  NF.  Except at
leading order (LO), result is obtained using a fourth-order
Runge-Kutta integration.

Todo:
    * Implement alpha strong running beyond LO
"""

from math import log


def beta(p: int, nf: int) -> float:
    """QCD beta function coefficient.

    Args:
        p (int): pQCD order, 0=LO, 1=NLO, ...
       nf (int): number of active quark flavors

    Returns:
        float: QCD beta function coefficient

    Examples:
        >>> beta(0, 3)
        -9.0

    """
    # Colour factors

    CF = 4.0/3.0
    CA = 3.0
    TF = 0.5

    B00 = 11./3. * CA
    B01 = -4./3. * TF
    B10 = 34./3. * CA**2
    B11 = -20./3. * CA*TF - 4. * CF*TF

    if p == 0:
        beta = - B00 - B01 * nf
    elif p == 1:
        beta = - B10 - B11 * nf
    elif p == 2:
        beta = - 1428.50 + 279.611 * nf - 6.01852 * nf**2
    elif p == 3:
        beta = - 29243.0 + 6946.30 * nf - 405.089 * nf**2 - 1.49931 * nf**3
    else:
        raise ValueError('N^5LO not yet implemented :-)')

    return beta


def as2pf(p: int, nf: int,  r2: float, as0: float, r20: float) -> float:
    """QCD beta function coefficient.

    Args:
        p: pQCD order, 0=LO, 1=NLO, ...
       nf: number of active quark flavors
       r2: final momentum scale squared
      as0: initial value for a_strong/(2 Pi)
      r20: initial momentum scale squared

    Returns:
        float: final value for a_strong/(2 Pi)

    Examples:
        >>> as2pf(0, 3, 8., 0.3, 4.)
        0.15497879500975464

    """
    # a below is as defined in 4pi expansion and is returned to
    # 2pi expansion convention just before return

    if p != 0:
        raise ValueError('Only LO implemented!')

    a = 0.5 * as0
    lrrat = log(r2/r20)

    a = 0.5 * as0 / (1. - 0.5 * beta(0, nf) * as0 * lrrat)

    # Return to .../(2pi)  expansion
    a = 2 * a

    return a
