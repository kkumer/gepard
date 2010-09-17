""" Quadrature formulas """

from numpy import sum, ones, add, diff, isinf, isscalar, \
     asarray, real, trapz, arange, empty, ndarray
from scipy.special.orthogonal import p_roots

[roots4, weights4] = p_roots(4)  # roots and weigths for 4-th order quadrature
roots4 = real(roots4)

[roots5, weights5] = p_roots(5)  # roots and weigths for 5-th order quadrature
roots5 = real(roots5)

[roots10, weights10] = p_roots(10)  # roots and weigths for 10-th order quadrature
roots10 = real(roots10)

[roots18, weights18] = p_roots(18)  # roots and weigths for 10-th order quadrature
roots18 = real(roots18)

[roots35, weights35] = p_roots(35)  # roots and weigths for 35-th order quadrature
roots35 = real(roots35)

def quadSciPy35(func,a,b,args=()):
    """Compute a definite integral using 35-order Gaussian quadrature.
    Adapted from scipy."""
    y = (b-a)*(roots35+1)/2.0 + a
    return (b-a)/2.0*sum(weights35*func(y,*args),0)

def quadSciPy18(func,a,b,args=()):
    """Compute a definite integral using 18-order Gaussian quadrature.
    Adapted from scipy."""
    y = (b-a)*(roots18+1)/2.0 + a
    return (b-a)/2.0*sum(weights18*func(y,*args),0)

def quadSciPy10(func,a,b,args=()):
    """Compute a definite integral using tenth-order Gaussian quadrature.
    Adapted from scipy."""
    y = (b-a)*(roots10+1)/2.0 + a
    return (b-a)/2.0*sum(weights10*func(y,*args),0)

def quadSciPy10transposed(func,a,b,args=()):
    """Compute a definite integral using fifth-order Gaussian quadrature.
    Adapted from scipy."""
    y = (b-a)*(roots10+1)/2.0 + a
    return (b-a)/2.0*sum((weights10*func(y,*args)).transpose(),0)

def quad10(func,a,b,args=()):
    """Compute a definite integral using tenth-order Gaussian quadrature.
    Adapted from scipy. Slower, but doesn't require ndarray from our code."""
    int = 0.0
    for i in range(len(roots10)):
        y = (b-a)*(roots10[i]+1)/2.0 + a
        int = int + (b-a)/2.0*(weights10[i]*func(y,*args))
    return int

def quadSciPy5(func,a,b,args=()):
    """Compute a definite integral using fifth-order Gaussian quadrature.
    Adapted from scipy."""
    y = (b-a)*(roots5+1)/2.0 + a
    return (b-a)/2.0*sum(weights5*func(y,*args),0)

def quadSciPy5transposed(func,a,b,args=()):
    """Compute a definite integral using fifth-order Gaussian quadrature.
    Adapted from scipy."""
    y = (b-a)*(roots5+1)/2.0 + a
    return (b-a)/2.0*sum((weights5*func(y,*args)).transpose(),0)

def quadSciPy4(func,a,b,args=()):
    """Compute a definite integral using fifth-order Gaussian quadrature.
    Adapted from scipy."""
    y = (b-a)*(roots4+1)/2.0 + a
    return (b-a)/2.0*sum(weights4*func(y,*args),0)

# Choice of routine used for P.V. integration
PVquadrature = quadSciPy18

# Choice of routine used for harmonic projection
Hquadrature = quadSciPy10transposed

# Choice of routine used for t-integration
tquadrature = quadSciPy5

