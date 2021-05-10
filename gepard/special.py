"""Special functions of complex arguments.

poch: Pochhammer symbol
  S1: harmonic sum
"""

from scipy.special import psi   # type: ignore

from gepard.constants import EMC


def pochhammer(z: complex, m: int) -> complex:
    """Pochhammer symbol.

    Args:
        z (complex):
        m (int):

    Returns:
        complex: pochhammer(z,m)

    Examples:
        >>> pochhammer(3.7+2j, 3)
        (42.723000000000006+122.54j)

    """
    p = z
    for k in range(1, m):
        p = p * (z + k)
    return p


def S1(z: complex) -> complex:
    """Harmonic sum S_1."""

    S1 = EMC + psi(z+1)

    return S1
