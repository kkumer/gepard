""" Definitions of models. 

"Model" is a set of functions (typically CFFs or GPDs) which depend on
parameters (some of which can be provided by minimization routine in the
fitting procedure).  Theoretical "approach", when given an instance of a model
and parameter values can calculate observables.

"""
#from IPython.Debugger import Tracer; debug_here = Tracer()
import pickle, sys

from numpy import log, pi, imag, real
from numpy import ndarray, array
from termcolor import colored

from quadrature import PVquadrature
from utils import flatten, hubDict

import pygepard as g
import optModel

class Model(object):
    """Base class for all models."""

    def __init__(self, optimization=False):
        self.optimization = optimization
        # Intially all parameters are fixed and should be released by user
        exec('fixed = {' + ", ".join(map(lambda x: "'fix_%s': %s" % x, 
                    zip(self.parameter_names, len(self.parameter_names)*['True']))) + '}')
        self.parameters.update(fixed)
        self.ndparameters = array([self.parameters[name] for name in self.parameter_names])
        #self.res = array([0. for k in range(18)])

    def release_parameters(self, *args):
        """Release parameters for fitting.

        If allowed ranges have to be changed from default, user needs
        to modify parameters_dict dictionary directly.

        """
        for par in args:
            if par not in self.parameter_names:
                raise ValueError('Parameter "%s" is not defined in model %s' 
                        % (par, self))
            self.parameters['fix_'+par] = False

    def fix_parameters(self, *args):
        """Fix parameters so they are not fitting variables."""
        if args[0] == 'ALL':
            # fix 'em all
            for par in self.parameter_names:
                self.parameters['fix_'+par] = True
        else:
            for par in args:
                if par not in self.parameter_names:
                    raise ValueError('Parameter "%s" is not defined in model %s' 
                            % (par, self))
                self.parameters['fix_'+par] = True

    def print_parameters(self, compare_with=[]):
        """Pretty-print parameters and their values.

        Variable parameters are printed green, while parameters with values
        at the limits of their range are printed red.
        If additional models are given in compare_with list, their parameter
        values are also printed and differences larger than 5% are denoted
        by blue and red coloring.

        """
        s = ""
        for name in self.parameter_names:
            value = self.parameters[name]
            row = '%4s -> %-5.3g' % (name, value)
            if self.parameters.has_key('limit_'+name):
                lo, hi = self.parameters['limit_'+name]
                if (abs((lo-value)*(hi-value)) < 0.001):
                    row = colored(row, 'red')
            if self.parameters['fix_'+name] == False:
                row = colored(row, 'green')
            for model in compare_with:
                try:
                    value2 =  model.parameters[name]
                except KeyError:
                    # compared model doesnt' have this parameter
                    value2 = 0
                app = '   %-5.3g' % value2
                # calculate relative diff, or absolute if value is zero
                diff = value - value2
                if value != 0:
                    diff = diff/abs(value)
                if diff > 0.05:
                    app = colored(app, 'red')
                elif diff < -0.05:
                    app = colored(app, 'blue')
                row += app
            row += '\n'
            s += row
        print s


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

    allCFFs = ['ImH', 'ReH', 'ImE', 'ReE', 'ImHt', 'ReHt', 'ImEt', 'ReEt']
    allCFFeffs = ['ImHeff', 'ReHeff', 'ImEeff', 'ReEeff', 
                     'ImHteff', 'ReHteff', 'ImEteff', 'ReEteff']

    def CFFvalues(self, pt):
        """Print values of CFFs. Pastable into Mathematica."""
        vals = map(lambda cff: str(getattr(self, cff)(pt)), ComptonFormFactors.allCFFs)
        s = "{" + 8*"%s -> %s, "
        s = s[:-2] + "}"
        return s % flatten(tuple(zip(ComptonFormFactors.allCFFs, vals)))

    # Initial definition of all CFFs. All just return zero.
    for name in allCFFs:
        exec('def %s(self, pt): return 0.' % name)

    for name in allCFFeffs:
        exec('def %s(self, pt): return 0.' % name)

class ComptonDispersionRelations(ComptonFormFactors):
    """Use dispersion relations for ReH and ReE

    methods: ReH, ReE, ReHt, ReEt, subtraction
    Subclass should implement ansaetze for ImH, ImE, ImHt, ImEt 
    and subtraction. This class implements just dispersion integrals.
    
    """

    def dispargV(self, x, fun, pt):
        """ Integrand of the dispersion integral (vector case) 
        
        fun -- Im(CFF)
        With variable change x->x^(1/(1-ga))=u in
        order to tame the singularity at x=0. 
        
        """
        ga = 0.9  # Nice value obtained by experimentation in Mathematica
        u = x**(1./(1.-ga))
        res = u**ga * ( fun(pt, u) - fun(pt) )
        return (2.*u) / (pt.xi**2 - u**2) * res / (1.-ga)

    def dispargA(self, x, fun, pt):
        """ Integrand of the dispersion integral (axial-vector case)
        
        fun -- Im(CFF)
        With variable change x->x^(1/(1-ga))=u in
        order to tame the singularity at x=0. 
        
        """
        ga = 0.9  # Value same as for V-case (FIXME: is this the best choice?)
        u = x**(1./(1.-ga))
        res = u**ga * ( fun(pt, u) - fun(pt) )
        return (2.* pt.xi) / (pt.xi**2 - u**2) * res / (1.-ga)

    def subtraction(self, pt):
        return 0  # default

    def ReH(self, pt):
        """ Real part of CFF H, 
        
        Given by dispersion integral over ImH - subtraction constant.
        
        """
        if self.optimization:
            pvpi = optModel.pvquadrature('H', self.ndparameters, pt.t, pt.xi)
            return pvpi - optModel.subtractionr(self.ndparameters, pt.t)
        else:
            res = PVquadrature(self.dispargV, 0, 1, (self.ImH, pt))
            pvpi = (res + log(pt.xi**2 / (1.-pt.xi**2)) * self.ImH(pt)) / pi
            # P.V./pi - subtraction constant C/(1-t/MC^2)^2
            return pvpi - self.subtraction(pt)

    def ReHt(self, pt):
        """ Real part of CFF Ht. 
        
        Given by dispersion integral over ImHt

        """
        if self.optimization:
            pvpi = optModel.pvquadrature('Ht', self.ndparameters, pt.t, pt.xi)
        else:
            res = PVquadrature(self.dispargA, 0, 1, (self.ImHt, pt))
            pvpi = (res + log((1.+pt.xi)/(1.-pt.xi)) * self.ImHt(pt))/pi
        return pvpi   # this is P.V./pi 

    def ReE(self, pt):
        """Real part of CFF E.
        
        Given by dispersion integral over ImE + subtraction constant.
        
        """
        if self.optimization:
            pvpi = optModel.pvquadrature('E', self.ndparameters, pt.t, pt.xi)
            return pvpi + optModel.subtractionr(self.ndparameters, pt.t)
        else:
            res = PVquadrature(self.dispargV, 0, 1, (self.ImE, pt))
            pvpi = (res + log(pt.xi**2 / (1.-pt.xi**2)) * self.ImE(pt)) / pi
            # This is same subtraction constant
            # as for H, but with opposite sign
            return pvpi + self.subtraction(pt)

    def ReEt(self, pt):
        """ Real part of CFF Et. 
        
        Given by dispersion integral over ImEt

        """
        if self.optimization:
            pvpi = optModel.pvquadrature('Et', self.ndparameters, pt.t, pt.xi)
        else:
            res = PVquadrature(self.dispargA, 0, 1, (self.ImEt, pt))
            pvpi = (res + log((1.+pt.xi)/(1.-pt.xi)) * self.ImEt(pt))
        return pvpi   # this is P.V./pi 


class ComptonModelDR(ComptonDispersionRelations):
    """Model for CFFs as in arXiv:0904.0458."""

    def __init__(self, **kwargs):
        # initial values of parameters and limits on their values
        self.parameters = {
              'NS' : 1.5,                                 
             'alS' : 1.13,                              
            'alpS' : 0.15,                              
              'MS' : 0.707,                               
              'rS' : 1.0,                               
              'bS' : 2.0,     'limit_bS' : (0.4, 5.0),
              'Nv' : 1.35,                              
             'alv' : 0.43,                              
            'alpv' : 0.85,                              
              'Mv' : 1.0,     'limit_Mv' : (0.4, 1.5),
              'rv' : 0.5,     'limit_rv' : (0., 8.),
              'bv' : 2.2,     'limit_bv' : (0.4, 5.),
               'C' : 7.0,      'limit_C' : (-10., 10.),
              'MC' : 1.3,     'limit_MC' : (0.4, 2.),
             'tNv' : 0.0,                             
             'tMv' : 2.7,    'limit_tMv' : (0.4, 2.),
             'trv' : 6.0,    'limit_trv' : (0., 8.),
             'tbv' : 3.0,    'limit_tbv' : (0.4, 5.)   }

        # order matters to fit.MinuitFitter, so it is defined by:
        self.parameter_names = ['NS', 'alS', 'alpS', 'MS', 'rS', 'bS',
                                'Nv', 'alv', 'alpv', 'Mv', 'rv', 'bv',
                                'C', 'MC',
                                'tNv', 'tMv', 'trv', 'tbv']

        # now do whatever else is necessary
        ComptonFormFactors.__init__(self, **kwargs)

    def subtraction(self, pt):
        return self.parameters['C']/(1.-pt.t/self.parameters['MC']**2)**2

    def ImH(self, pt, xi=0):
        """Imaginary part of CFF H."""
        p = self.parameters # just a shortcut
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
        val = ( (2.*4./9. + 1./9.) * p['Nv'] * p['rv'] * twox**(-p['alv']-p['alpv']*t) *
                 onex**p['bv'] / (1. - onex*t/(p['Mv']**2))  )
        sea = ( (2./9.) * p['NS'] * p['rS'] * twox**(-p['alS']-p['alpS']*t) *
                 onex**p['bS'] / (1. - onex*t/(p['MS']**2))**2 )
        return pi * (val + sea) / (1.+x)


    def ImHt(self, pt, xi=0):
        """Imaginary part of CFF Ht i.e. \tilde{H}."""
        p = self.parameters # just a shortcut
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
        val = ( (2.*4./9. + 1./9.) * p['tNv'] * p['trv'] * 
            # Regge trajectory params taken from H:
            twox**(-p['alv']-p['alpv']*t) *
                 onex**p['tbv'] / (1. - onex*t/(p['tMv']**2))  )
        return pi * val / (1.+x)

    def ImE(self, pt, xi=0):
        """Imaginary part of CFF E."""
        # Just changing function signature w.r.t. ComptonFormFactors
        # to make it compatible for dispersion integral
        return 0

    def ReEt(self, pt):
        """Instead of disp. rel. use pole formula."""
        return (2.2390424 * (1. - (1.7*(0.0196 - pt.t))/(1. 
            - pt.t/2.)**2))/((0.0196 - pt.t)*pt.xi)


class ComptonModelDRPP(ComptonModelDR):
    """Model for CFFs as in arXiv:0904.0458. + free pion pole normalization"""

    def __init__(self, **kwargs):
        # initial values of parameters and limits on their values
        self.parameters = {
              'NS' : 1.5,                                 
             'alS' : 1.13,                              
            'alpS' : 0.15,                              
              'MS' : 0.707,                               
              'rS' : 1.0,                               
              'bS' : 2.0,     'limit_bS' : (0.4, 5.0),
              'Nv' : 1.35,                              
             'alv' : 0.43,                              
            'alpv' : 0.85,                              
              'Mv' : 1.0,     'limit_Mv' : (0.4, 1.5),
              'rv' : 0.5,     'limit_rv' : (0., 8.),
              'bv' : 2.2,     'limit_bv' : (0.4, 5.),
               'C' : 7.0,      'limit_C' : (-10., 10.),
              'MC' : 1.3,     'limit_MC' : (0.4, 2.),
             'tNv' : 0.0,                             
             'tMv' : 2.7,    'limit_tMv' : (0.4, 2.),
             'trv' : 6.0,    'limit_trv' : (0., 8.),
             'tbv' : 3.0,    'limit_tbv' : (0.4, 5.),
             'NPP' : 1.0,    'limit_tPP' : (-8, 8.)   }

        # order matters to fit.MinuitFitter, so it is defined by:
        self.parameter_names = ['NS', 'alS', 'alpS', 'MS', 'rS', 'bS',
                                'Nv', 'alv', 'alpv', 'Mv', 'rv', 'bv',
                                'C', 'MC',
                                'tNv', 'tMv', 'trv', 'tbv', 'NPP']

        # now do whatever else is necessary
        ComptonFormFactors.__init__(self, **kwargs)

    def ReEt(self, pt):
        """Instead of disp. rel. use pole formula * NPP."""
        return self.parameters['NPP']*(2.2390424 * (1. - (1.7*(0.0196 - pt.t))/(1. 
            - pt.t/2.)**2))/((0.0196 - pt.t)*pt.xi)

class ComptonModelDRsea(ComptonDispersionRelations):
    """DR Model intended for combining with Gepard sea. NS->Nsea"""

    def __init__(self, **kwargs):
        # initial values of parameters and limits on their values
        self.parameters = {
              'Nsea' : 1.5,                                 
             'alS' : 1.13,                              
            'alpS' : 0.15,                              
              'MS' : 0.707,                               
              'rS' : 1.0,                               
              'bS' : 2.0,     'limit_bS' : (0.4, 5.0),
              'Nv' : 1.35,                              
             'alv' : 0.43,                              
            'alpv' : 0.85,                              
              'Mv' : 1.0,     'limit_Mv' : (0.4, 1.5),
              'rv' : 0.5,     'limit_rv' : (0., 8.),
              'bv' : 2.2,     'limit_bv' : (0.4, 5.),
               'C' : 7.0,      'limit_C' : (-10., 10.),
              'MC' : 1.3,     'limit_MC' : (0.4, 2.),
             'tNv' : 0.0,                             
             'tMv' : 2.7,    'limit_tMv' : (0.4, 2.),
             'trv' : 6.0,    'limit_trv' : (0., 8.),
             'tbv' : 3.0,    'limit_tbv' : (0.4, 5.)   }

        # order matters to fit.MinuitFitter, so it is defined by:
        self.parameter_names = ['Nsea', 'alS', 'alpS', 'MS', 'rS', 'bS',
                                'Nv', 'alv', 'alpv', 'Mv', 'rv', 'bv',
                                'C', 'MC',
                                'tNv', 'tMv', 'trv', 'tbv']

        # now do whatever else is necessary
        ComptonFormFactors.__init__(self, **kwargs)

    def subtraction(self, pt):
        return self.parameters['C']/(1.-pt.t/self.parameters['MC']**2)**2

    def ImH(self, pt, xi=0):
        """Imaginary part of CFF H."""
        p = self.parameters # just a shortcut
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
        val = ( (2.*4./9. + 1./9.) * p['Nv'] * p['rv'] * twox**(-p['alv']-p['alpv']*t) *
                 onex**p['bv'] / (1. - onex*t/(p['Mv']**2))  )
        sea = ( (2./9.) * p['Nsea'] * p['rS'] * twox**(-p['alS']-p['alpS']*t) *
                 onex**p['bS'] / (1. - onex*t/(p['MS']**2))**2 )
        return pi * (val + sea) / (1.+x)

    def ImHt(self, pt, xi=0):
        """Imaginary part of CFF Ht i.e. \tilde{H}."""
        p = self.parameters # just a shortcut
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
        val = ( (2.*4./9. + 1./9.) * p['tNv'] * p['trv'] * 
            # Regge trajectory params taken from H:
            twox**(-p['alv']-p['alpv']*t) *
                 onex**p['tbv'] / (1. - onex*t/(p['tMv']**2))  )
        return pi * val / (1.+x)

    def ImE(self, pt, xi=0):
        """Imaginary part of CFF E."""
        # Just changing function signature w.r.t. ComptonFormFactors
        # to make it compatible for dispersion integral
        return 0

    def ReEt(self, pt):
        """Instead of disp. rel. use pole formula."""
        return (2.2390424 * (1. - (1.7*(0.0196 - pt.t))/(1. 
            - pt.t/2.)**2))/((0.0196 - pt.t)*pt.xi)


class ComptonNeuralNets(Model):
    """Neural network CFFs"""

    def __init__(self, hidden_layers=[7], output_layer=['ImH', 'ReH'], endpointpower=None):
        """Model CFFs by neural networks.
        
        Neural network, created actually by Fitter instance, will have
        architecture:
        input layer: two neurons for xB and t
        hidden_layers: keyword argument hidden_layers is a list specifying
                       number of neurons in consequtive hidden layers (so
                       default is just one hidden layer with 7 neurons)
        output_layer:  keyword argument output_layer is a list specifying
                       names of CFFs will be given by neural nets. Rest are
                       zero.
        endpointpower: CFFs are defined as NN*(1-xB)**endpointpower to enforce
                       vanishing at xB=0 and to improve convergence
        
        """
        self.architecture = [2] + hidden_layers + [len(output_layer)]
        self.output_layer = output_layer
        self.nets = []
        self.parameters = {'nnet':0, 'outputvalue':None}
        self.parameter_names = ['nnet', 'outputvalue']
        self.endpointpower = endpointpower
        # now do whatever else is necessary
        #ComptonFormFactors.__init__(self)

    def __getattr__(self, name):
        """Return appropriate CFF function object."""
        # FIXME: I don't understand why I have to use this:
        if name in object.__getattribute__(self, 'output_layer'):
        # and this creates infinite recursion:
        #if name in self.output_layer:
            self.curname = name
            return self.CFF
        elif name in ComptonFormFactors.allCFFs:
            # if asked for CFF which is not in output_layer, return 0
            return self.zero
        elif name in ComptonFormFactors.allCFFeffs:
            # if asked for CFF which is not in output_layer, return 0
            return self.zero
        elif name == 'endpointpower':
            if self.__dict__.has_key('endpointpower'):
                return self.endpointpower
            else:
                return None
        else:
            raise AttributeError

    def zero(self, *args, **kwargs):
        return 0

    def CFF(self, pt, outputvalue=None, xi=0):
        # FIXME: This function is HEAVILY sub-optimal and non-pythonic!
        ind = self.output_layer.index(self.curname)
        if self.parameters['outputvalue'] != None:
            # this occurs during training: value is set by training
            # routine by calling with outputvalue set by training routine
            return self.parameters['outputvalue'][ind]
        if isinstance(xi, ndarray):
            # function was called with third argument that is xi nd array
            x = xi
        elif xi != 0:
            # function was called with third argument that is xi number
            x = array((xi,))
        else:
            # xi should be taken from pt object
            if isinstance(pt.xi, ndarray):
                x = pt.xi
            else:
                x = array((pt.xi,))
        xBs = 2.*x/(1.+x)
        res = []
        for xB in xBs:
            ar = []
            for net in self.nets:
                if self.endpointpower:
                    ar.append(net.activate([xB, pt.t])[ind]*(1-xB)**self.endpointpower)
                else:
                    ar.append(net.activate([xB, pt.t])[ind])
            all = array(ar).flatten()
            if self.parameters.has_key('nnet'):
                if self.parameters['nnet'] == 'ALL':
                    res.append(all)
                elif self.parameters['nnet'] == 'AVG':
                    res.append(all.mean())
                else: # we want particular net
                    try:
                        res.append(all[self.parameters['nnet']])
                    except IndexError:
                        raise IndexError, str(self)+' has only '+str(len(self.nets))+' nets!'
            # by default, we get mean value (FIXME:this should never occurr?)
            else:
                res.append(all.mean())
        res = array(res)
        if res.shape == (1,):
            # returns number
            return res[0]  
        else:
            # returns ndarray
            return res.transpose()


class ComptonGepard(ComptonFormFactors):
    """CFFs as implemented in gepard. 

    cutq2 - Q2 at which evolution is frozen (default = 0 GeV^2)
    
    """
    def __init__(self, cutq2=0.0, **kwargs):
        # initial values of parameters and limits on their values
        self.parameters = {
               'NS' : 0.15,
             'AL0S' : 1.0,
             'ALPS' : 0.15,
             'M02S' : 1.0,    'limit_M02S' : (0.3, 1.5),
           'DELM2S' : 0.0,
               'PS' : 2.0,
             'SECS' : 0.0,
             'THIS' : 0.0,
             'KAPS' : 0.0,
            'SKEWS' : 0.0,
               'NG' : 0.5,
             'AL0G' : 1.1,
             'ALPG' : 0.15,
             'M02G' : 0.7,    'limit_M02G' : (0.3, 1.5),
           'DELM2G' : 0.0,
               'PG' : 2.0,
             'SECG' : 0.0,
             'THIG' : 0.0,
             'KAPG' : 0.0,
            'SKEWG' : 0.0   }


        # gepard needs indices, not parameter names
        self.parameters_index = {
             11 : 'NS',
             12 : 'AL0S',
             13 : 'ALPS',
             14 : 'M02S',
             15 : 'DELM2S',
             16 : 'PS',
             17 : 'SECS',
             18 : 'KAPS',
             19 : 'SKEWS',
             21 : 'NG',  
             22 : 'AL0G',
             23 : 'ALPG',
             24 : 'M02G',
             25 : 'DELM2G',
             26 : 'PG',
             27 : 'SECG',
             28 : 'KAPG',
             29 : 'SKEWG',
             32 : 'THIS',
             42 : 'THIG' }

        # order matters to fit.MinuitFitter, so it is defined by:
        self.parameter_names = [ 
           'NS', 'AL0S', 'ALPS', 'M02S',
           'DELM2S', 'PS', 'SECS', 'THIS', 'KAPS', 'SKEWS',
           'NG', 'AL0G', 'ALPG', 'M02G',
           'DELM2G', 'PG', 'SECG', 'THIG', 'KAPG', 'SKEWG']

        # this was in Gepard's GEPARD.INI, which is not needed now
        # but look at it for documentation of what parameters below are
        g.parint.speed = 1
        g.parint.acc = 3
        g.parint.p = 0
        g.parint.nf = 4
        g.parint.czero = 1

        g.astrong.mu02 = 2.5
        g.astrong.asp = array([0.0606, 0.0518, 0.0488])

        g.parflt.q02 = 4.0
        g.parflt.rf2 = 1.0
        g.parflt.rr2 = 1.0

        g.mbcont.c = 0.35
        g.mbcont.phi = 1.57079632
        g.mbcont.cnd = -0.25
        g.mbcont.phind = 1.57

        g.parchr.scheme = array([c for c in 'CSBAR'])  # array(5)
        g.parchr.ansatz = array([c for c in 'FIT   ']) # array(6)

        # following two items usually came from driver file
        g.parchr.process = array([c for c in 'DVCS  '])  # array(6)
        g.parchr.fftype = array([c for c in 'SINGLET   ']) # array(10)

        g.init()
        # Cutting-off evolution  at Q2 = cutq2
        # Evaluate evolved C at this scale now.
        self.cutq2 = cutq2
        g.nqs.nqs = 1
        g.qs.qs[0] = self.cutq2
        g.kinematics.q2 = self.cutq2
        g.evolc(1, 1)
        g.evolc(2, 1)
        g.evolc(3, 1)
        self.qdict = {self.cutq2 : 1}

        g.newcall = 1
        # number of points on MB contour
        g.npts = 2**g.parint.acc * 12 / g.parint.speed
        self.g = g
        # now do whatever else is necessary
        ComptonFormFactors.__init__(self, **kwargs)

    def _evolve(self, pt):
        """Calculate evolution operator."""
        self.g.nqs.nqs += 1
        nqs = int(self.g.nqs.nqs)
        self.g.qs.qs[nqs-1] = pt.Q2
        self.qdict[pt.Q2] = nqs
        if pt.Q2 < self.cutq2:
            # just copy the evolved C from Q2=cutq2
            for k in range(self.g.npts):
                # both partial waves; for quarks and gluons:
                self.g.cgrid.cgrid[0,nqs-1,k,0] = self.g.cgrid.cgrid[0,0,k,0]
                self.g.cgrid.cgrid[1,nqs-1,k,0] = self.g.cgrid.cgrid[1,0,k,0]
                self.g.cgrid.cgrid[0,nqs-1,k,1] = self.g.cgrid.cgrid[0,0,k,1]
                self.g.cgrid.cgrid[1,nqs-1,k,1] = self.g.cgrid.cgrid[1,0,k,1]
        else:
            self.g.evolc(1, nqs)
            self.g.evolc(2, nqs)
            self.g.evolc(3, nqs)

    def _GepardCFFs(self, pt, xi=0):
        """Call gepard routine that calculates CFFs."""
        for i in self.parameters_index:
            g.par.par[i-1] = self.parameters[self.parameters_index[i]]

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

        g.kinematics.q2 = pt.Q2
        g.kinematics.xi = x
        g.kinematics.del2 = pt.t

        if not self.qdict.has_key(pt.Q2):
            self._evolve(pt)

        g.mt.nmts = 1
        g.mt.mtind = 0 
        g.mts.mts[0] = - pt.t

        g.cfff()
        g.newcall = 0

    def ImH(self, pt, xi=0):
        """Imaginary part of CFF H."""
        if self.g.newcall:
            self._GepardCFFs(pt, xi)
        return imag(g.cff.cff[g.parint.p])

    def ReH(self, pt):
        """Real part of CFF H."""
        if self.g.newcall:
            self._GepardCFFs(pt)
        return real(g.cff.cff[g.parint.p])

    def ImE(self, pt):
        """Imaginary part of CFF E."""
        if self.g.newcall:
            self._GepardCFFs(pt)
        return imag(g.cff.cffe[g.parint.p])

    def ReE(self, pt):
        """Real part of CFF E."""
        if self.g.newcall:
            self._GepardCFFs(pt)
        return real(g.cff.cffe[g.parint.p])


class ComptonHybrid(ComptonFormFactors):
    """This combines gepard for small xB and DR model for valence xB."""

    def __init__(self, instGepard, instDR, **kwargs):
        self.Gepard = instGepard  # instance of ComptonGepard
        self.DR = instDR  # instance of ComptonModelDR
        self.DR.parameters['Nsea'] = 0.  # sea comes from Gepard part
        self.parameters = hubDict(self.Gepard.parameters, self.DR.parameters)
        self.parameter_names = self.Gepard.parameter_names + self.DR.parameter_names

        self.g = g
        # now do whatever else is necessary
        ComptonFormFactors.__init__(self, **kwargs)


    def ImH(self, pt, xi=0):
        return  self.Gepard.ImH(pt, xi) + self.DR.ImH(pt, xi)

    def ReH(self, pt):
        return  self.Gepard.ReH(pt) + self.DR.ReH(pt)

    def ImE(self, pt, xi=0):
        return  self.Gepard.ImE(pt) + self.DR.ImE(pt, xi)

    def ReE(self, pt):
        return  self.Gepard.ReE(pt) + self.DR.ReE(pt)

    # tildes are not provided by Gepard

    def ImHt(self, pt, xi=0):
        return  self.DR.ImHt(pt, xi)

    def ReHt(self, pt):
        return  self.DR.ReHt(pt)

    def ImEt(self, pt):
        return  self.DR.ImEt(pt)

    def ReEt(self, pt):
        return  self.DR.ReEt(pt)

##  --- Complete models built from the above components ---

class ModelDR(ComptonModelDR, ElasticDipole):
    """Complete model as in arXiv:0904.0458.."""


class ModelNN(ComptonNeuralNets, ElasticDipole):
    """Complete model."""


class Hybrid(ComptonHybrid, ElasticDipole):
    """Complete hybrid model."""

class ModelDRPP(ComptonModelDRPP, ElasticDipole):
    """Complete model as in arXiv:0904.0458. + free pion pole normalization."""

