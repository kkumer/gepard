""" Implementation of dispersion relations

FIXME: Duplicate DR code should be removed from Model.py

"""

from numpy import log, pi

from quadrature import PVquadrature

def dispV(x, fun, pt):
    """ Integrand of the dispersion integral (vector case) 
    
      x --  integration variable
    fun --  spectral function, i.e., Im(CFF)
     pt --  DataPoint instance specifying kinematics
    
    """
    # We use variable change x->x^(1/(1-ga))=u in
    # order to tame the singularity at x=0. 
    ga = 0.9  # Nice value obtained by experimentation
    u = x**(1./(1.-ga))
    res = u**ga * ( fun(pt, u) - fun(pt) )
    #sys.stderr.write('u = %f, fun = %f' % (u, fun(pt,u)))
    return (2.*u) / (pt.xi**2 - u**2) * res / (1.-ga)

def dispA(x, fun, pt):
    """ Integrand of the dispersion integral (axial-vector case)
    
      x --  integration variable
    fun --  spectral function, i.e., Im(CFF)
     pt --  DataPoint instance specifying kinematics
    
    """
    ga = 0.9  # Value same as for V-case (FIXME: is this the best choice?)
    u = x**(1./(1.-ga))
    res = u**ga * ( fun(pt, u) - fun(pt) )
    return (2.* pt.xi) / (pt.xi**2 - u**2) * res / (1.-ga)


def intV(fun, pt):
    """ Vector P.V. dispersion integral of fun divided by pi (without subtraction!)"""
    res = PVquadrature(dispV, 0, 1, (fun, pt))
    pvpi = (res + log(pt.xi**2 / (1.-pt.xi**2)) * fun(pt)) / pi
    return pvpi 

def intVNN(fun, pt):
    """ Vector P.V. dispersion integral of fun divided by pi (without subtraction!)"""
    res = PVquadrature(dispV, 0, 1, (fun, pt))
    pvpi = (res.reshape(res.size,1) + log(pt.xi**2 / (1.-pt.xi**2)) * fun(pt)) / pi
    return pvpi 


def intA(fun, pt):
    """ Axial-vector P.V. dispersion integral of fun divided by pi"""
    res = PVquadrature(dispA, 0, 1, (fun, pt))
    pvpi = (res + log((1.+pt.xi)/(1.-pt.xi)) * fun(pt))/pi
    return pvpi 

