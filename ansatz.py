""" Here actual ansatz for Compton Form Factor (CFF) H is defined, as
well as ansaetze for Dirac and Pauli elastic form factors F_1 and F_2"""

from numpy import log, pi
from numpy import ndarray

from quadrature import PVquadrature


### Dirac and Pauli proton form factors ###

# Dipole approximation

def F1(t):
    return (1.41 * (1.26 - t))/((0.71 - t)**2 * (3.53 - t))

def F2(t):
    return 3.2 / ((0.71 - t)**2 * (3.53 - t))


### CFFs ###

def ImcffH(pt, pars, xi=0):
    """Imaginary part of CFF H."""
    (NS, alS, alpS, M2S, sS, bS, Nv, alv, alpv, M2v, sv, bv, C, M2C) = pars
    if isinstance(xi, ndarray):
        # function was called with third argument that is xi nd array
        x = xi
    elif xi != 0:
        # function was called with third argument that is xi number
        x = xi
    else:
        # xi should be taken from pt object
        x = pt.xi
    t = pt.t
    twox = 2.*x / (1.+x)
    onex = (1.-x) / (1.+x)
    val = ( (2.*4./9. + 1./9.) * Nv * (1.+sv) * twox**(-alv-alpv*t) *
             onex**bv / (1. - onex*t/M2v)  )
    sea = ( (2./9.) * NS * (1.+sS) * twox**(-alS-alpS*t) *
             onex**bS / (1. - onex*t/M2S)**2 )
    return pi * (val + sea)

def disparg(x, pt, pars):
    """ Integrand of the dispersion integral with variable change x->x^(1/(1-ga))=u in
    order to tame the singularity at x=0. """
    ga = 0.9  # Nice value obtained by experimentation in Mathematica
    u = x**(1./(1.-ga))
    res = u**ga * ( ImcffH(pt, pars, u) - ImcffH(pt, pars) )
    return (2.*u) / (pt.xi**2 - u**2) * res / (1.-ga)

def RecffH(pt, pars):
    """ Real part of CFF H, obtained from imaginary part using dispersion relation."""
    res = PVquadrature(disparg, 0, 1, (pt, pars))
    pv = res + log(pt.xi**2 / (1.-pt.xi**2)) * ImcffH(pt, pars)
    return pv/pi - pars[-2]/(1.-pt.t/pars[-1])**2  # this is P.V./pi - subtraction constant C/(1-t/M2c)^2

## Simpler integration but less accurate
# def dispargOLD(x, pt, pars):
#     """ Integrand of the dispersion integral. """
#     res = ImcffH(pt, pars, x) - ImcffH(pt, pars)
#     return (2.*x)/(pt.xi**2 - x**2) * res
# 
# def RecffHOLD(pt, pars):
#     """ Real part using dispersion relation """
#     res = PVquadrature(dispargOLD, 0, 1, (pt, pars))
#     pv = res + log(pt.xi**2 / (1.-pt.xi**2)) * ImcffH(pt, pars)
#     return pv/pi - pars[-1]  # this is P.V./pi - subtraction constant C
