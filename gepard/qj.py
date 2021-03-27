"""Building blocks for GPD models.

FIXME: ptj ModuleType --> DataPoint
"""

import gepard.special


def qj(j: complex, t: float, poch: int, norm: float, al0: float, alp: float) -> complex:
    r"""GPD building block Q_j with reggeized t-dependence.

    Args:
        j: complex conformal moment
        t: momentum transfer
        poch: pochhammer (beta?)
        al0: Regge intercept
        alp: Regge slope

    Returns:
        conformal moment of GPD with Reggeized t-dependence,

    Examples:
        >>> gepard.qj.qj(0.5+0.3j, 0.,  4, 3, 0.5, -0.8)  #doctest: +ELLIPSIS
        (5.668107685...-4.002820...j)

    Notes:
        From Eq. (40) or (41) of arXiv:0904.0458
        but without residual beta(t) from Eq. (19):

        .. math::

            Q_j = N \frac{B(1-\alpha_0+j, \beta+1)}{B(2-\alpha_0, \beta+1)}
                  \frac{1+j-\alpha_0}{1+j-\alpha(t)}

    """
    alpt = al0 + alp * t
    qj = (
        norm
        * gepard.special.pochhammer(2 - al0, poch)
        / gepard.special.pochhammer(1 - al0 + j, poch)
        * (1 + j - al0)
        / (1 + j - alpt)
    )
    return qj


def betadip(j: complex, t: float, m02: float, delm2: float, pp: int) -> complex:
    r"""GPD residual dipole t-dependence function beta from Eq. (19) of NPB.

    Args:
        j: conformal moment
        t: momentum transfer
        m02: mass param
        delm2: mass param - interplay with j
        pp: exponent (=2 for dipole)

    Returns:
        Dipole residual t-dependence
    """
    return 1. / (1. - t / (m02 + delm2*j))**pp
