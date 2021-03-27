"""Special functions of complex arguments.

poch: Pochhammer symbol
"""

# from cmath import exp
# from scipy.special import loggamma as clngamma


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
