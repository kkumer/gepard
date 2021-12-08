"""Building blocks for GPD models.

FIXME: ptj ModuleType --> DataPoint
"""

import numpy as np

from . import special


def qj(j: np.ndarray, t: float, poch: int, norm: float, al0: float,
        alp: float, alpf: float = 0, val: int = 0) -> np.ndarray:
    r"""GPD building block Q_j with reggeized t-dependence.

    Args:
        j: complex conformal moment
        t: momentum transfer
        poch: pochhammer (beta+1)
        al0: Regge intercept
        alp: Regge slope (leading pole)
        alpf: Regge slope (full)
        val: 0 for sea, 1 for valence partons

    Returns:
        conformal moment of GPD with Reggeized t-dependence,

    Examples:
        >>> gepard.qj.qj(0.5+0.3j, 0.,  4, 3, 0.5, -0.8)  #doctest: +ELLIPSIS
        (5.668107685...-4.002820...j)

    Notes:
        There are two implementations of Regge t dependence, one uses
        `alpf` :math:`\alpha'` parameter (`alp` should be set to zero)

        .. math::

            Q_j = N \frac{B(1-\alpha(t)+j, \beta+1)}{B(2-\alpha_0, \beta+1)}

        This is then equal to Eq. (41) of arXiv:hep-ph/0605237
        without residual :math:`F(\Delta^2)`.

        The default uses `alp` parameter, and expression is
        from Eq. (40) or (41) of arXiv:0904.0458, which takes into account
        only the leading pole:

        .. math::

            Q_j = N \frac{B(1-\alpha_0+j, \beta+1)}{B(2-\alpha_0, \beta+1)}
                  \frac{1+j-\alpha_0}{1+j-\alpha(t)}

        and is again defined without residual beta(t) from Eq. (19).


    """
    assert alp*alpf == 0   # only one version of alpha' can be used
    alpt = al0 + alp * t
    qj = (
        norm
        * special.pochhammer(2 - val - al0 - alpf * t, poch)
        / special.pochhammer(1 - al0 + j, poch)
        * (1 + j - al0)
        / (1 + j - alpt)
    )
    return qj


def betadip(j: np.ndarray, t: float, m02: float, delm2: float, pp: int) -> np.ndarray:
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
