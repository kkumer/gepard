r"""Perturbative QCD functions.

Coefficients of beta function of QCD up to order :math:`\alpha_s^5` (NNNLO),
normalized by

.. math::

     \frac{d a_s}{d \ln \mu_r^2}  =   \beta_0 a_s^2 + \beta_1 a_s^3 + \cdots

with

.. math::

    a_s = \frac{\alpha_{s}}{4\pi}.

Beyond NLO the QCD colour factors are hard-wired in this routine,
and the numerical coefficients are truncated to six digits.

QCD coupling alpha strong is obtained by integrating the evolution
equation for a fixed number of massless flavours  NF.  Except at
leading order (LO), result is obtained using a fourth-order
Runge-Kutta integration.
"""

from math import log

from . import constants


def beta(p: int, nf: int) -> float:
    """QCD beta function coefficient.

    Args:
        p (int): pQCD order, 0=LO, 1=NLO, 2=NNLO
        nf (int): number of active quark flavors

    Returns:
        float: QCD beta function coefficient

    Examples:
        >>> beta(0, 3)
        -9.0

    """
    B00 = 11./3. * constants.CA
    B01 = -4./3. * constants.TF
    B10 = 34./3. * constants.CA**2
    B11 = -20./3. * constants.CA*constants.TF - 4. * constants.CF*constants.TF

    if p == 0:
        beta = - B00 - B01 * nf
    elif p == 1:
        beta = - B10 - B11 * nf
    elif p == 2:
        beta = - 1428.50 + 279.611 * nf - 6.01852 * nf**2
    elif p == 3:
        beta = - 29243.0 + 6946.30 * nf - 405.089 * nf**2 - 1.49931 * nf**3
    else:
        raise ValueError('NNNNLO not yet implemented :-)')

    return beta


def _fbeta1(a: float, nf: int) -> float:
    return a**2 * (beta(0, nf) + a * beta(1, nf))


def as2pf(p: int, nf: int,  r2: float, as0: float, r20: float) -> float:
    """QCD beta function coefficient.

    Args:
        p: pQCD order, 0=LO, 1=NLO, 2=NNLO
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
    # a below is as defined in 1/4pi expansion and is returned to
    # 1/2pi expansion convention just before return
    NASTPS = 20

    a = 0.5 * as0
    lrrat = log(r2/r20)
    dlr = lrrat / NASTPS

    if p == 0:
        a = 0.5 * as0 / (1. - 0.5 * beta(0, nf) * as0 * lrrat)
    elif p == 1:
        for k in range(1, NASTPS+1):
            xk0 = dlr * _fbeta1(a, nf)
            xk1 = dlr * _fbeta1(a + 0.5 * xk0, nf)
            xk2 = dlr * _fbeta1(a + 0.5 * xk1, nf)
            xk3 = dlr * _fbeta1(a + xk2, nf)
            a = a + (xk0 + 2 * xk1 + 2 * xk2 + xk3) / 6
    else:
        raise ValueError('Only LO and NLO implemented!')

    # Return to .../(2pi)  expansion
    a = 2 * a

    return a
