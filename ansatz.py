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
    (NS, alS, alpS, M2S, rS, bS, Nv, alv, alpv, M2v, rv, bv, C, M2C, 
            tNv, tM2v, trv, tbv) = pars
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
    val = ( (2.*4./9. + 1./9.) * Nv * rv * twox**(-alv-alpv*t) *
             onex**bv / (1. - onex*t/M2v)  )
    sea = ( (2./9.) * NS * rS * twox**(-alS-alpS*t) *
             onex**bS / (1. - onex*t/M2S)**2 )
    return pi * (val + sea) / (1.+x)

def dispargV(x, fun, pt, pars):
    """ Integrand of the dispersion integral with variable change x->x^(1/(1-ga))=u in
    order to tame the singularity at x=0. """
    ga = 0.9  # Nice value obtained by experimentation in Mathematica
    u = x**(1./(1.-ga))
    res = u**ga * ( fun(pt, pars, u) - fun(pt, pars) )
    return (2.*u) / (pt.xi**2 - u**2) * res / (1.-ga)

def RecffH(pt, pars):
    """ Real part of CFF H, obtained from imaginary part using dispersion relation."""
    res = PVquadrature(dispargV, 0, 1, (ImcffH, pt, pars))
    pv = res + log(pt.xi**2 / (1.-pt.xi**2)) * ImcffH(pt, pars)
    return pv/pi - pars[-6]/(1.-pt.t/pars[-5])**2  # this is P.V./pi - subtraction constant C/(1-t/M2c)^2

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

def ImcffHt(pt, pars, xi=0):
    """Imaginary part of CFF Ht i.e. \tilde{H}."""

    (NS, alS, alpS, M2S, rS, bS, Nv, alv, alpv, M2v, rv, bv, C, M2C, 
            tNv, tM2v, trv, tbv) = pars
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
    val = ( (2.*4./9. + 1./9.) * tNv * trv * twox**(-alv-alpv*t) *
             onex**tbv / (1. - onex*t/tM2v)  )
    return pi * val / (1.+x)

def dispargA(x, fun, pt, pars):
    """ Integrand of the dispersion integral with variable change x->x^(1/(1-ga))=u in
    order to tame the singularity at x=0. """
    ga = 0.9  # Nice value obtained by experimentation in Mathematica
    u = x**(1./(1.-ga))
    res = u**ga * ( fun(pt, pars, u) - fun(pt, pars) )
    return (2.* pt.xi) / (pt.xi**2 - u**2) * res / (1.-ga)

def RecffHt(pt, pars):
    """ Real part of CFF Ht, obtained from imaginary part using dispersion relation."""

    res = PVquadrature(dispargA, 0, 1, (ImcffHt, pt, pars))
    pv = res + log((1.+pt.xi)/(1.-pt.xi)) * ImcffHt(pt, pars)
    return pv/pi   # this is P.V./pi 

def ImcffE(pt, pars):
        return 0.

def ImcffEt(pt, pars):
        return 0.

def RecffE(pt, pars):
    return pars[-6]/(1.-pt.t/pars[-5])**2  # this is just subtraction constant C/(1-t/M2c)^2

def RecffEt(pt, pars):
        return (2.2390424 * (1. - (1.7*(0.0196 - pt.t))/(1. - pt.t/2.)**2)
                  )/((0.0196 - pt.t)*pt.xi)


