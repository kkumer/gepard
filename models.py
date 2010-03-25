""" Here are definitions of models. "Model" is a set of functions
(typically CFFs or GPDs) which depend on parameters (some of which
can be provided by minimization routine in the fitting procedure).
Theoretical "approach", when given an instance of a model and
parameter values can calculate observables.
Classes:
    Model -- base class
    FormFactors(Model) -- dipole FFs and dispersion-relation CFFs
                          with ansatz as in arXiv:0904.0458

"""
import pickle, sys

from numpy import log, pi
from numpy import ndarray, array

from quadrature import PVquadrature
from utils import AttrDict, flatten


def dispargV(x, fun, pt, pars):
    """ Integrand of the dispersion integral (vector case) 
    
    fun -- Im(CFF)
    With variable change x->x^(1/(1-ga))=u in
    order to tame the singularity at x=0. 
    
    """
    ga = 0.9  # Nice value obtained by experimentation in Mathematica
    u = x**(1./(1.-ga))
    res = u**ga * ( fun(pt, pars, u) - fun(pt, pars) )
    return (2.*u) / (pt.xi**2 - u**2) * res / (1.-ga)

def dispargA(x, fun, pt, pars):
    """ Integrand of the dispersion integral (axial-vector case)
    
    fun -- Im(CFF)
    With variable change x->x^(1/(1-ga))=u in
    order to tame the singularity at x=0. 
    
    """
    ga = 0.9  # Value same as for V-case (FIXME: is this the best choice?)
    u = x**(1./(1.-ga))
    res = u**ga * ( fun(pt, pars, u) - fun(pt, pars) )
    return (2.* pt.xi) / (pt.xi**2 - u**2) * res / (1.-ga)


class Model(object):
    """Later some methods or attributes may be added here."""

    def printCFFs(self, pt, pars={}):
        cffs = ['ImH', 'ReH', 'ImE', 'ReE', 'ImHt', 'ReHt', 'ImEt', 'ReEt']
        vals = map(lambda cff: str(getattr(self, cff)(pt, pars)), cffs)
        s = "{" + 8*"%s -> %s, "
        s = s[:-2] + "}"
        return s % flatten(tuple(zip(cffs, vals)))



class FormFactors(Model):
    """Compton and elastic Form Factors.

    methods: ImH, ReH, ImE, ReE, ImHt, ReHt, ImEt, ReEt, F1, F2
    attributes:  H, E, Ht, Et
    
    """

    # Attribute-access parameter dictionaries with default values
    H = AttrDict({
          'NS' : 1.5,       
         'alS' : 1.13,      
        'alpS' : 0.15,      
          'MS' : 0.707,       
          'rS' : 1.0,       
          'bS' : 2.0,   
          'Nv' : 1.35,      
         'alv' : 0.43,      
        'alpv' : 0.85,      
          'Mv' : 1.011,
          'rv' : 0.496383,  
          'bv' : 2.15682,   
           'C' : 6.90484,   
          'MC' : 1.33924
          })

    Ht = AttrDict({
          'Nv' : 0.6,      
          'Mv' : 2.69667,
          'rv' : 5.97923,  
          'bv' : 3.25607
          })
    #parsE = AttrDict() # not needed
    #parsEt = AttrDict() # not needed
    pars = AttrDict({'H':H, 'Ht':Ht}) # this one holds all parameter dicts

    
    def ImH(self, pt, pars={}, xi=0):
        """Imaginary part of CFF H."""

        # override defaults (FIXME: should this be permanent?)
        # FIXME: Is it costly to do this update for every call of ImH?
        self.H.update(pars) 
        p = self.H # shortcut

        # FIXME: The following solution is not elegant
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
        val = ( (2.*4./9. + 1./9.) * p.Nv * p.rv * twox**(-p.alv-p.alpv*t) *
                 onex**p.bv / (1. - onex*t/(p.Mv**2))  )
        sea = ( (2./9.) * p.NS * p.rS * twox**(-p.alS-p.alpS*t) *
                 onex**p.bS / (1. - onex*t/(p.MS**2))**2 )
        return pi * (val + sea) / (1.+x)

    def ReH(self, pt, pars={}):
        """ Real part of CFF H, 
        
        obtained from imaginary part using dispersion relation.
        
        """
        # override defaults (FIXME: should this be permanent?)
        self.H.update(pars) 
        p = self.H # shortcut

        res = PVquadrature(dispargV, 0, 1, (self.ImH, pt, pars))
        pv = res + log(pt.xi**2 / (1.-pt.xi**2)) * self.ImH(pt, pars)
        # P.V./pi - subtraction constant C/(1-t/MC^2)^2
        return pv/pi - p.C/(1.-pt.t/p.MC**2)**2  

    def ImHt(self, pt, pars={}, xi=0):
        """Imaginary part of CFF Ht i.e. \tilde{H}."""

        # override defaults 
        self.Ht.update(pars) 
        p = self.Ht # shortcut

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
        val = ( (2.*4./9. + 1./9.) * p.tNv * p.trv * 
            # Regge trajectory params taken from H:
            twox**(-self.H.alv-self.H.alpv*t) *
                 onex**p.tbv / (1. - onex*t/(p.tMv**2))  )
        return pi * val / (1.+x)

    def ReHt(self, pt, pars={}):
        """ Real part of CFF Ht. 
        
        Obtained from imaginary part using dispersion relation.

        """
        res = PVquadrature(dispargA, 0, 1, (self.ImHt, pt, pars))
        pv = res + log((1.+pt.xi)/(1.-pt.xi)) * self.ImHt(pt, pars)
        return pv/pi   # this is P.V./pi 

    def ImE(self, pt, pars={}):
        return 0

    def ReE(self, pt, pars={}):
        """Real part of CFF {\cal E}."""

        # This is same subtraction constant C/(1-t/M2c)^2
        # as for H, but with opposite sign
        return self.H.C/(1.-pt.t/self.H.MC**2)**2 

    def ImEt(self, pt, pars={}):
        return 0

    def ReEt(self, pt, pars={}):
            return (2.2390424 * (1. - (1.7*(0.0196 - pt.t))/(1. 
                - pt.t/2.)**2))/((0.0196 - pt.t)*pt.xi)

    def F1(self, t):
        """Dirac elastic form factor. 
        
        Dipole approximation from DM's notebook.
        
        """
        return (1.41 * (1.26 - t))/((0.71 - t)**2 * (3.53 - t))

    def F2(self, t):
        """Pauli elastic form factor.
        
        Dipole approximation from DM's notebook.

        """
        return 3.2 / ((0.71 - t)**2 * (3.53 - t))


class NNFormFactors(FormFactors):
    """Neural network CFF H.
    F1 and F2 are taken from FormFactors. Im(CFF H) is given by neural
    nets, while Re(CFF H) and other GPDs are zero.
    
    """
    def __init__(self):
        global nets
        nets = pickle.load(open('nets.pkl', 'r'))
        sys.stderr.write('Neural nets loaded from nets.pkl')
    
    def ImH(self, pt, pars={}, xi=0):
        ar = []
        for net in nets:
            ar.append(net.activate([pt.xB, pt.t]))
        return array(ar).flatten()


    def ReH(self, pt, pars={}):
        return 0

    def ImHt(self, pt, pars={}, xi=0):
        return 0

    def ReHt(self, pt, pars={}):
        return 0

    def ImE(self, pt, pars={}):
        return 0

    def ReE(self, pt, pars={}):
        return 0

    def ImEt(self, pt, pars={}):
        return 0

    def ReEt(self, pt, pars={}):
        return 0

class Hdominance(FormFactors):

    def ReE(self, pt, pars={}):
        return 0

    def ReEt(self, pt, pars={}):
        return 0
