"""Special functions of complex arguments.

poch: Pochhammer symbol
  S1: harmonic sum
"""

from typing import Union

import numpy as np
from scipy.special import psi  # type: ignore

from gepard.constants import EMC


def pochhammer(z: Union[complex, np.ndarray], m: int) -> Union[complex, np.ndarray]:
    """Pochhammer symbol.

    Args:
        z (complex):
        m (int):

    Returns:
        complex: pochhammer(z,m)
    """
    p = z
    for k in range(1, m):
        p = p * (z + k)
    return p


def S1(z: Union[complex, np.ndarray]) -> Union[complex, np.ndarray]:
    """Harmonic sum S_1."""
    S1 = EMC + psi(z+1)
    return S1
