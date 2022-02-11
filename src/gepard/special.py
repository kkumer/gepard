"""Special functions of complex arguments.

poch: Pochhammer symbol
dpsi: polygamma function
Sn: harmonic sums
"""

from math import factorial, log
from typing import Union

import numpy as np
from scipy.special import psi, zeta  # type: ignore


def parity(n: int) -> int:
    """Parity of n, i. e. (-1)^n."""
    return 1 - 2 * (n % 2)


def pochhammer(z: Union[complex, np.ndarray], m: int) -> Union[complex, np.ndarray]:
    """Pochhammer symbol.

    Args:
        z: complex argument
        m: integer index

    Returns:
        complex: pochhammer(z,m)

    """
    p = z
    for k in range(1, m):
        p = p * (z + k)
    return p


poch = pochhammer  # just an abbreviation


def dpsi_one(z: complex, m: int) -> complex:
    """Polygamma - m'th derivative of Euler gamma at z."""
    # Algorithm from Vogt, cf. julia's implementation
    sub = 0j

    if z.imag < 10:
        subm = (-1/z)**(m+1) * factorial(m)
        while z.real < 10:
            sub += subm
            z += 1
            subm = (-1/z)**(m+1) * factorial(m)

    a1 = 1.
    a2 = 1./2.
    a3 = 1./6.
    a4 = -1./30.
    a5 = 1./42.
    a6 = -1./30.
    a7 = 5./66.

    if m != 1:
        for k2 in range(2, m+1):
            a1 = a1 * (k2-1)
            a2 = a2 * k2
            a3 = a3 * (k2+1)
            a4 = a4 * (k2+3)
            a5 = a5 * (k2+5)
            a6 = a6 * (k2+7)
            a7 = a7 * (k2+9)

    rz = 1. / z
    dz = rz * rz
    res = (sub + (-1)**(m+1) * rz**m *
           (a1 + rz * (a2 + rz * (a3 + dz *
            (a4 + dz * (a5 + dz * (a6 + a7 * dz)))))))
    return res


dpsi = np.vectorize(dpsi_one)


def S1(z: Union[complex, np.ndarray]) -> Union[complex, np.ndarray]:
    """Harmonic sum S_1."""
    return np.euler_gamma + psi(z+1)


def S2(z: Union[complex, np.ndarray]) -> Union[complex, np.ndarray]:
    """Harmonic sum S_2."""
    return zeta(2) - dpsi(z+1, 1)


def S3(z: Union[complex, np.ndarray]) -> Union[complex, np.ndarray]:
    """Harmonic sum S_3."""
    return zeta(3) + dpsi(z+1, 2) / 2


def S4(z: Union[complex, np.ndarray]) -> Union[complex, np.ndarray]:
    """Harmonic sum S_4."""
    return zeta(4) - dpsi(z+1, 3) / 6


def S2_prime(z: Union[complex, np.ndarray], prty: int) -> Union[complex, np.ndarray]:
    """Curci et al Eq. (5.25)."""
    # note this is related to delS2
    return (1+prty)*S2(z)/2 + (1-prty)*S2(z-1/2)/2


def S3_prime(z: Union[complex, np.ndarray], prty: int) -> Union[complex, np.ndarray]:
    """Curci et al Eq. (5.25)."""
    return (1+prty)*S3(z)/2 + (1-prty)*S3(z-1/2)/2


def delS2(z: Union[complex, np.ndarray]) -> Union[complex, np.ndarray]:
    """Return Harmonic sum S_2 difference.

    Args:
        z: complex argument

    Returns:
        delS2((z+1)/2) From Eq. (4.13) of 1310.5394.
        Note halving of the argument.

    """
    return S2(z) - S2(z - 1/2)


def deldelS2(j: Union[complex, np.ndarray], k: int) -> Union[complex, np.ndarray]:
    """Return diference of harmonic sum S_2 differences.

    Args:
        j: complex argument
        k: integer index

    Returns:
        Equal to delS2((j+1)/2, (k+1)/2) From Eq. (4.38) of 1310.5394.
        Note halving of the argument.

    """
    return (delS2(j) - delS2(k)) / (4*(j-k)*(2*j+2*k+1))


def Sm1(z: Union[complex, np.ndarray], k: int) -> Union[complex, np.ndarray]:
    """Aux fun. FIXME: not tested."""
    return - log(2) + 0.5 * (1-2*(k % 2)) * (psi((z+2)/2) - psi((z+1)/2))


def MellinF2(n: Union[complex, np.ndarray]) -> Union[complex, np.ndarray]:
    """Return Mellin transform  i.e. x^(N-1) moment of Li2(x)/(1+x).

    Args:
        n: complex argument

    Returns:
        According to Eq. (33) in Bluemlein and Kurth, hep-ph/9708388

    """
    abk = np.array([0.9999964239, -0.4998741238,
                    0.3317990258, -0.2407338084, 0.1676540711,
                   -0.0953293897, 0.0360884937, -0.0064535442])
    psitmp = psi(n)
    mf2 = 0

    for k in range(1, 9):
        psitmp = psitmp + 1 / (n + k - 1)
        mf2 += (abk[k-1] *
                ((n - 1) * (zeta(2) / (n + k - 1) -
                 (psitmp + np.euler_gamma) / (n + k - 1)**2) +
                 (psitmp + np.euler_gamma) / (n + k - 1)))

    return zeta(2) * log(2) - mf2


def SB3(j: Union[complex, np.ndarray]) -> Union[complex, np.ndarray]:
    """Eq. (4.44e) of arXiv:1310.5394."""
    return 0.5*S1(j)*(-S2(-0.5+0.5*j)+S2(0.5*j))+0.125*(-S3(
             - 0.5 + 0.5 * j) + S3(0.5 * j)) - 2 * (0.8224670334241131 * (
                 -S1(0.5 * (-1 + j)) + S1(0.5 * j)) - MellinF2(1 + j))


def S2_tilde(n: Union[complex, np.ndarray], prty: int) -> Union[complex, np.ndarray]:
    """Eq. (30) of  Bluemlein and Kurth, hep-ph/9708388."""
    G = psi((n+1)/2) - psi(n/2)
    return -(5/8)*zeta(3) + prty*(S1(n)/n**2 - (zeta(2)/2)*G + MellinF2(n))
