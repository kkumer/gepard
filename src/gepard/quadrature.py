"""Quadrature formulas, abscissas and weights.

Also numerical derivative routine.
"""

from typing import Tuple

import numpy as np
from scipy.special.orthogonal import p_roots

roots4, weights4 = p_roots(4)  # roots and weigths for 4-th order quadrature

roots5, weights5 = p_roots(5)

roots10, weights10 = p_roots(10)

roots18, weights18 = p_roots(18)

roots35, weights35 = p_roots(35)

roots81, weights81 = p_roots(81)


def mellin_barnes(c, phi, accuracy: int = 3) -> Tuple[np.ndarray, np.ndarray]:
    """Construct basic MB array.

    Args:
        c: point of crossing the real axis
        phi: angle with the real axis
        accuracy: 2*accuracy points per interval

    Returns:
        (complex coordinates of MB contour, integration weighths)

    """
    c = 0.35
    phij = phi*1j
    roots, weights = p_roots(2**accuracy)
    division = np.array([0., 0.01, 0.08, 0.15, 0.3, 0.5, 1.0,
                         1.5, 2.0, 4.0, 6.0, 8.0, 10.0])
    summ = division[1:] + division[:-1]
    diff = division[1:] - division[:-1]
    x = []
    wg = []
    dv = 0
    for dv in range(len(division)-1):
        x.append((diff[dv]*roots + summ[dv])/2)
        wg.append(diff[dv]*weights/2)
        x_array = np.array(x).flatten()
        wg_array = np.array(wg).flatten()
    n_array = c + 1 + x_array * np.exp(phij)
    return n_array, wg_array


def nd_mellin_barnes(accuracy: int = 3) -> Tuple[np.ndarray, np.ndarray]:
    """Construct MB array for non-diagonal evolution operator.

    Args:
        accuracy: 2*accuracy points per interval

    Returns:
        (complex coordinates of MB contour, integration weighths)

    Note:
        Code should be merged with standard mellin_barnes.

    """
    cnd = -0.25
    phij = 1.57079632j
    roots, weights = p_roots(2**accuracy)
    division = np.array([0., 0.01, 0.025, 0.067, 0.18, 0.5, 1.3,
                         3.7, 10, 25, 67, 180, 500])
    summ = division[1:] + division[:-1]
    diff = division[1:] - division[:-1]
    x = []
    wg = []
    dv = 0
    for dv in range(len(division)-1):
        x.append((diff[dv]*roots + summ[dv])/2)
        wg.append(diff[dv]*weights/2)
        x_array = np.array(x).flatten()
        wg_array = np.array(wg).flatten()
    znd_array = cnd + x_array * np.exp(phij)
    return znd_array, wg_array


def quadSciPy81(func, a, b, args=()):
    """Compute a definite integral using 81-order Gaussian quadrature."""
    y = (b-a)*(roots81+1)/2.0 + a
    return (b-a)/2.0*np.sum(weights81*func(y, *args), 0)


def quadSciPy35(func, a, b, args=()):
    """Compute a definite integral using 35-order Gaussian quadrature."""
    y = (b-a)*(roots35+1)/2.0 + a
    return (b-a)/2.0*np.sum(weights35*func(y, *args), 0)


def quadSciPy18transposed(func, a, b, args=()):
    """Compute a definite integral using 18-order Gaussian quadrature."""
    y = (b-a)*(roots18+1)/2.0 + a
    return (b-a)/2.0*np.sum((weights18*func(y, *args)).transpose(), 0)


def quadSciPy10(func, a, b, args=()):
    """Compute a definite integral using tenth-order Gaussian quadrature."""
    y = (b-a)*(roots10+1)/2.0 + a
    return (b-a)/2.0*np.sum(weights10*func(y, *args), 0)


def quadSciPy10transposed(func, a, b, args=()):
    """Compute a definite integral using fifth-order Gaussian quadrature."""
    y = (b-a)*(roots10+1)/2.0 + a
    return (b-a)/2.0*np.sum((weights10*func(y, *args)).transpose(), 0)


def quad10(func, a, b, args=()):
    """Compute a definite integral using tenth-order Gaussian quadrature."""
    # Slower, but doesn't require ndarray from our code.
    int = 0.0
    for i in range(len(roots10)):
        y = (b-a)*(roots10[i]+1)/2.0 + a
        int = int + (b-a)/2.0*(weights10[i]*func(y, *args))
    return int


def quadSciPy5(func, a, b, args=()):
    """Compute a definite integral using fifth-order Gaussian quadrature."""
    y = (b-a)*(roots5+1)/2.0 + a
    return (b-a)/2.0*sum(weights5*func(y, *args), 0)


def quadSciPy5transposed(func, a, b, args=()):
    """Compute a definite integral using fifth-order Gaussian quadrature."""
    y = (b-a)*(roots5+1)/2.0 + a
    return (b-a)/2.0*sum((weights5*func(y, *args)).transpose(), 0)


def quadSciPy4(func, a, b, args=()):
    """Compute a definite integral using fifth-order Gaussian quadrature."""
    y = (b-a)*(roots4+1)/2.0 + a
    return (b-a)/2.0*sum(weights4*func(y, *args), 0)


# Choice of routine used for P.V. integration
PVquadrature = quadSciPy18transposed

# Choice of routine used for harmonic projection
Hquadrature = quadSciPy10transposed

# Choice of routine used for t-integration
tquadrature = quadSciPy5

# Choice of routine used for t->b Fourier transform
bquadrature = quadSciPy10

# Choice of routine used for <r(theta)> integrals
rthtquadrature = quadSciPy10


def deriv(func, x, h, nevals):
    """Return derivative func'(x) using Ridders-Neville algorithm."""
    # Adapted from Press et al., Numerical Recipes
    # in a blind, non-pythonic way
    con = 1.4  # scale decrease per step
    safe = 2   # return when error is safe worse than the best so far
    big = 1.e3
    hh = h
    a = np.zeros((nevals+1, nevals+1))
    a[1, 1] = (func(x+hh) - func(x-hh))/(2*hh)
    err = big
    for i in range(2, nevals+1):
        hh = hh / con
        a[1, i] = (func(x+hh) - func(x-hh))/(2*hh)
        fac = con**2
        for j in range(2, i+1):
            a[j, i] = (a[j-1, i]*fac - a[j-1, i-1])/(fac-1)
            fac = con**2 * fac
            errt = max(abs(a[j, i]-a[j-1, i]),
                       abs(a[j, i]-a[j-1, i-1]))
            if (errt <= err):
                err = errt
                drv = a[j, i]
        if abs(a[i, i]-a[i-1, i-1]) >= safe*err:
            return drv, err
