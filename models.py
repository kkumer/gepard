""" Definitions of models. 

"Model" is a set of functions (typically CFFs or GPDs) which depend on
parameters (some of which can be provided by minimization routine in the
fitting procedure).  Theoretical "approach", when given an instance of a model
and parameter values can calculate observables.

"""
import pickle, sys

from numpy import log, pi
from numpy import ndarray, array

from quadrature import PVquadrature
from utils import AttrDict, flatten


class Model(object):
    """Later some methods or attributes may be added here."""


class ElasticFormFactors(Model):
    """Dirac and Pauli elastic form factors F_1 and F_2."""


class ElasticDipole(ElasticFormFactors):
    """Dipole approximation from DM's notebook."""

    def F1(self, t):
        """Dirac elastic form factor."""
        return (1.41 * (1.26 - t))/((0.71 - t)**2 * (3.53 - t))

    def F2(self, t):
        """Pauli elastic form factor."""
        return 3.2 / ((0.71 - t)**2 * (3.53 - t))


class ComptonFormFactors(Model):
    """Twist-two, no-transversity set of 4 CFFs.

    They are set to be zero here. Actual models are built by subclassing this.

    """
    funcnames = ['ImH', 'ReH', 'ImE', 'ReE', 'ImHt', 'ReHt', 'ImEt', 'ReEt']
    # Initial definition of CFFs. All just return zero.
    for name in funcnames:
        exec('def %s(self, pt, pars={}): return 0.' % name)


    def values(self, pt, pars={}):
        """Print values of CFFs. Pastable into Mathematica."""
        vals = map(lambda cff: str(getattr(self, cff)(pt, pars)), self.funcnames)
        s = "{" + 8*"%s -> %s, "
        s = s[:-2] + "}"
        return s % flatten(tuple(zip(self.funcnames, vals)))


class ComptonDispersionRelations(ComptonFormFactors):
    """Use dispersion relations for ReH and ReE

    methods: ReH, ReE, ReHt, ReEt, subtraction
    Subclass should implement ansaetze for ImH, ImE, ImHt, ImEt 
    and subtraction. This class implements just dispersion integrals.
    
    """

    def dispargV(self, x, fun, pt, pars):
        """ Integrand of the dispersion integral (vector case) 
        
        fun -- Im(CFF)
        With variable change x->x^(1/(1-ga))=u in
        order to tame the singularity at x=0. 
        
        """
        ga = 0.9  # Nice value obtained by experimentation in Mathematica
        u = x**(1./(1.-ga))
        res = u**ga * ( fun(pt, pars, u) - fun(pt, pars) )
        return (2.*u) / (pt.xi**2 - u**2) * res / (1.-ga)

    def dispargA(self, x, fun, pt, pars):
        """ Integrand of the dispersion integral (axial-vector case)
        
        fun -- Im(CFF)
        With variable change x->x^(1/(1-ga))=u in
        order to tame the singularity at x=0. 
        
        """
        ga = 0.9  # Value same as for V-case (FIXME: is this the best choice?)
        u = x**(1./(1.-ga))
        res = u**ga * ( fun(pt, pars, u) - fun(pt, pars) )
        return (2.* pt.xi) / (pt.xi**2 - u**2) * res / (1.-ga)

    def subtraction(self, pt, pars={}):
        return 0  # default

    def ReH(self, pt, pars={}):
        """ Real part of CFF H, 
        
        Given by dispersion integral over ImH - subtraction constant.
        
        """
        # override defaults (FIXME: should this be permanent?)
        self.pars.update(pars) 
        p = self.pars # shortcut

        res = PVquadrature(self.dispargV, 0, 1, (self.ImH, pt, p))
        pv = res + log(pt.xi**2 / (1.-pt.xi**2)) * self.ImH(pt, p)
        # P.V./pi - subtraction constant C/(1-t/MC^2)^2
        return pv/pi - self.subtraction(pt, p)

    def ReHt(self, pt, pars={}):
        """ Real part of CFF Ht. 
        
        Given by dispersion integral over ImHt

        """
        # override defaults 
        self.pars.update(pars) 
        p = self.pars # shortcut

        res = PVquadrature(self.dispargA, 0, 1, (self.ImHt, pt, p))
        pv = res + log((1.+pt.xi)/(1.-pt.xi)) * self.ImHt(pt, p)
        return pv/pi   # this is P.V./pi 

    def ReE(self, pt, pars={}):
        """Real part of CFF E.
        
        Given by dispersion integral over ImE + subtraction constant.
        
        """
        # override defaults 
        self.pars.update(pars) 
        p = self.pars # shortcut

        res = PVquadrature(self.dispargV, 0, 1, (self.ImE, pt, p))
        pv = res + log(pt.xi**2 / (1.-pt.xi**2)) * self.ImE(pt, p)
        # This is same subtraction constant
        # as for H, but with opposite sign
        return pv/pi + self.subtraction(pt, p)

    def ReEt(self, pt, pars={}):
        """ Real part of CFF Et. 
        
        Given by dispersion integral over ImEt

        """
        # override defaults 
        self.pars.update(pars) 
        p = self.pars # shortcut

        res = PVquadrature(self.dispargA, 0, 1, (self.ImEt, pt, p))
        pv = res + log((1.+pt.xi)/(1.-pt.xi)) * self.ImEt(pt, p)
        return pv/pi   # this is P.V./pi 


class ComptonModelDR(ComptonDispersionRelations):
    """Model for CFFs as in arXiv:0904.0458."""

    # Attribute-access parameter dictionary with default values
    pars = AttrDict({
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
          'MC' : 1.33924,
          'tNv' : 0.6,      
          'tMv' : 2.69667,
          'trv' : 5.97923,  
          'tbv' : 3.25607
          })

    def subtraction(self, pt, pars={}):
        return pars.C/(1.-pt.t/pars.MC**2)**2

    def ImH(self, pt, pars={}, xi=0):
        """Imaginary part of CFF H."""

        # override defaults (FIXME: should this be permanent?)
        # FIXME: Is it costly to do this update for every call of ImH?
        self.pars.update(pars) 
        p = self.pars # shortcut

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

    def ImHt(self, pt, pars={}, xi=0):
        """Imaginary part of CFF Ht i.e. \tilde{H}."""

        # override defaults 
        self.pars.update(pars) 
        p = self.pars # shortcut

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
            twox**(-p.alv-p.alpv*t) *
                 onex**p.tbv / (1. - onex*t/(p.tMv**2))  )
        return pi * val / (1.+x)

    def ImE(self, pt, pars={}, xi=0):
        """Imaginary part of CFF E."""
        # Just changing function signature w.r.t. ComptonFormFactors
        # to make it compatible for dispersion integral
        return 0

    def ReEt(self, pt, pars={}):
        """Instead of disp. rel. use pole formula."""
        return (2.2390424 * (1. - (1.7*(0.0196 - pt.t))/(1. 
            - pt.t/2.)**2))/((0.0196 - pt.t)*pt.xi)


class ComptonNNH(ComptonFormFactors):
    """Neural network CFF H.

    Im(CFF H) is given by neural nets in file 'nets.pkl', while 
    Re(CFF H) and other GPDs are zero.
    
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


##  --- Complete models built from the above components ---

class ModelDR(ComptonModelDR, ElasticDipole):
    """Complete model as in arXiv:0904.0458.."""


class ModelNN(ComptonNNH, ElasticDipole):
    """Complete model."""

