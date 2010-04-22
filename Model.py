""" Definitions of models. 

"Model" is a set of functions (typically CFFs or GPDs) which depend on
parameters (some of which can be provided by minimization routine in the
fitting procedure).  Theoretical "approach", when given an instance of a model
and parameter values can calculate observables.

"""
#from IPython.Debugger import Tracer; debug_here = Tracer()
import pickle, sys

from numpy import log, pi
from numpy import ndarray, array
from termcolor import colored

from quadrature import PVquadrature
from utils import AttrDict, flatten


class Model(object):
    """Base class for all models."""

    def __init__(self):
        # Intially all parameters are fixed and should be released by user
        exec('self.fixed = {' + ", ".join(map(lambda x: "'fix_%s': %s" % x, 
                    zip(self.parameter_names, len(self.parameter_names)*['True']))) + '}')
        self.parameters.update(self.fixed)

    def release_parameters(self, *args):
        """Release parameters for fitting.

        If allowed ranges have to be changed from default, user needs
        to modify parameters_dict dictionary directly.

        """
        for par in args:
            if par not in self.parameter_names:
                raise ValueError('Parameter "%s" is not defined in model %s' 
                        % (par, self))
            self.fixed['fix_'+par] = False
        self.parameters.update(self.fixed)

    def fix_parameters(self, *args):
        """Fix parameters so they are not fitting variables."""
        if args[0] == 'ALL':
            # fix 'em all
            for par in self.parameter_names:
                self.fixed['fix_'+par] = True
        else:
            for par in args:
                if par not in self.parameter_names:
                    raise ValueError('Parameter "%s" is not defined in model %s' 
                            % (par, self))
                self.fixed['fix_'+par] = True
        self.parameters.update(self.fixed)

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
            #if self.fixed['fix_'+name] == False:
            if self.parameters.has_key('limit_'+name):
                lo, hi = self.parameters['limit_'+name]
                if (abs((lo-value)*(hi-value)) < 0.001):
                    row = colored(row, 'red')
            if self.parameters['fix_'+name] == False:
                row = colored(row, 'green')
            for model in compare_with:
                value2 =  model.parameters[name]
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
    def __init__(self):
        self.allCFFs = ['ImH', 'ReH', 'ImE', 'ReE', 'ImHt', 'ReHt', 'ImEt', 'ReEt']
        # now do whatever else is necessary
        Model.__init__(self)



    def CFFvalues(self, pt):
        """Print values of CFFs. Pastable into Mathematica."""
        vals = map(lambda cff: str(getattr(self, cff)(pt)), self.allCFFs)
        s = "{" + 8*"%s -> %s, "
        s = s[:-2] + "}"
        return s % flatten(tuple(zip(self.allCFFs, vals)))


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
        res = PVquadrature(self.dispargV, 0, 1, (self.ImH, pt))
        pv = res + log(pt.xi**2 / (1.-pt.xi**2)) * self.ImH(pt)
        # P.V./pi - subtraction constant C/(1-t/MC^2)^2
        return pv/pi - self.subtraction(pt)

    def ReHt(self, pt):
        """ Real part of CFF Ht. 
        
        Given by dispersion integral over ImHt

        """

        res = PVquadrature(self.dispargA, 0, 1, (self.ImHt, pt))
        pv = res + log((1.+pt.xi)/(1.-pt.xi)) * self.ImHt(pt)
        return pv/pi   # this is P.V./pi 

    def ReE(self, pt):
        """Real part of CFF E.
        
        Given by dispersion integral over ImE + subtraction constant.
        
        """
        res = PVquadrature(self.dispargV, 0, 1, (self.ImE, pt))
        pv = res + log(pt.xi**2 / (1.-pt.xi**2)) * self.ImE(pt)
        # This is same subtraction constant
        # as for H, but with opposite sign
        return pv/pi + self.subtraction(pt)

    def ReEt(self, pt):
        """ Real part of CFF Et. 
        
        Given by dispersion integral over ImEt

        """
        res = PVquadrature(self.dispargA, 0, 1, (self.ImEt, pt))
        pv = res + log((1.+pt.xi)/(1.-pt.xi)) * self.ImEt(pt)
        return pv/pi   # this is P.V./pi 

    funcnames = ['ImH', 'ImE', 'ImHt', 'ImEt']
    # Initial definition of other CFFs. All just return zero.
    for name in funcnames:
        exec('def %s(self, pt): return 0.' % name)


class ComptonModelDRdict(ComptonDispersionRelations):
    """Model for CFFs as in arXiv:0904.0458."""

    def __init__(self):
        # initial values of parameters and limits on their values
        self.parameters = AttrDict({
              'NS' : 1.5,                                 
             'alS' : 1.13,                              
            'alpS' : 0.15,                              
              'MS' : 0.707,                               
              'rS' : 1.0,                               
              'bS' : 2.0,     'limit_bS' : (0.4, 5.0),
              'Nv' : 1.35,                              
             'alv' : 0.43,                              
            'alpv' : 0.85,                              
              'Mv' : 1.0,     'limit_Mv' : (0.9, 1.1),
              'rv' : 0.5,     'limit_rv' : (0., 8.),
              'bv' : 2.2,     'limit_bv' : (0.4, 5.),
               'C' : 7.0,      'limit_C' : (-10., 10.),
              'MC' : 1.3,     'limit_MC' : (0.4, 2.),
             'tNv' : 0.0,                             
             'tMv' : 2.7,    'limit_tMv' : (0.4, 2.),
             'trv' : 6.0,    'limit_trv' : (0., 8.),
             'tbv' : 3.0,    'limit_tbv' : (0.4, 5.)   })

        # order matters to fit.MinuitFitter, so it is defined by:
        self.parameter_names = ['NS', 'alS', 'alpS', 'MS', 'rS', 'bS',
                                'Nv', 'alv', 'alpv', 'Mv', 'rv', 'bv',
                                'C', 'MC',
                                'tNv', 'tMv', 'trv', 'tbv']

        # now do whatever else is necessary
        ComptonDispersionRelations.__init__(self)



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
        val = ( (2.*4./9. + 1./9.) * p.Nv * p.rv * twox**(-p.alv-p.alpv*t) *
                 onex**p.bv / (1. - onex*t/(p.Mv**2))  )
        sea = ( (2./9.) * p.NS * p.rS * twox**(-p.alS-p.alpS*t) *
                 onex**p.bS / (1. - onex*t/(p.MS**2))**2 )
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
        val = ( (2.*4./9. + 1./9.) * p.tNv * p.trv * 
            # Regge trajectory params taken from H:
            twox**(-p.alv-p.alpv*t) *
                 onex**p.tbv / (1. - onex*t/(p.tMv**2))  )
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


class ComptonModelDR(ComptonDispersionRelations):
    """Model for CFFs as in arXiv:0904.0458. -- no AttrDict!"""

    def __init__(self):
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
              'Mv' : 1.0,     'limit_Mv' : (0.9, 1.1),
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
        ComptonFormFactors.__init__(self)



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


class ComptonNNH(ComptonFormFactors):
    """Neural network CFF H -- both imaginary and real part given by nets."""

    def __init__(self, justImH=False):
        self.architecture = [2, 7, 2]
        self.justImH = justImH
        self.nets = []
        #sys.stderr.write('Neural nets loaded from nets.pkl')
        # nnet -  net index
        # outputvalue - given by neural net during training
        self.parameters = {'nnet':0, 'outputvalue':None}
        self.parameter_names = ['nnet', 'outputvalue']
        # now do whatever else is necessary
        ComptonFormFactors.__init__(self)

    def ImH(self, pt, outputvalue=None):
        if self.parameters['outputvalue'] != None:
            # this occurs during training
            return self.parameters['outputvalue'][0]
        ar = []
        for net in self.nets:
            ar.append(net.activate([pt.xB, pt.t])[0])
        all = array(ar).flatten()
        if self.parameters.has_key('nnet'):
            if self.parameters['nnet'] == 'ALL':
                return all
            elif self.parameters['nnet'] == 'AVG':
                return all.mean()
            else: # we want particular net
                try:
                    return all[self.parameters['nnet']]
                except IndexError:
                    raise IndexError, str(self)+' has only '+str(len(self.nets))+' nets!'
        # by default, we get mean value (FIXME:this should never occurr?)
        else:
            return all.mean()

    def ReH(self, pt, outputvalue=None):
        if self.justImH:
            return 0
        if self.parameters['outputvalue'] != None:
            # this occurs during training
            return self.parameters['outputvalue'][1]
        ar = []
        for net in self.nets:
            ar.append(net.activate([pt.xB, pt.t])[1])
        all = array(ar).flatten()
        if self.parameters.has_key('nnet'):
            if self.parameters['nnet'] == 'ALL':
                return all
            elif self.parameters['nnet'] == 'AVG':
                return all.mean()
            else: # we want particular net
                try:
                    return all[self.parameters['nnet']]
                except IndexError:
                    raise IndexError, str(self)+' has only '+str(len(self.nets))+' nets!'
        # by default, we get mean value
        else:
            return all.mean()

    funcnames = ['ImE', 'ReE', 'ImHt', 'ReHt', 'ImEt', 'ReEt']
    # Initial definition of other CFFs. All just return zero.
    for name in funcnames:
        exec('def %s(self, pt): return 0.' % name)


class ComptonNeuralNets(ComptonFormFactors):
    """Neural network CFFs"""

    def __init__(self, hidden_layers=[7], output_layer=['ImH', 'ReH']):
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
        
        """
        self.architecture = [2] + hidden_layers + [len(output_layer)]
        self.output_layer = output_layer
        self.nets = []
        self.parameters = {'nnet':0, 'outputvalue':None}
        self.parameter_names = ['nnet', 'outputvalue']
        # now do whatever else is necessary
        ComptonFormFactors.__init__(self)

    def __getattr__(self, name):
        """Return appropriate CFF function object."""
        # FIXME: I don't understand why I have to use this:
        if name in object.__getattribute__(self, 'output_layer'):
        # and this creates infinite recursion:
        #if name in self.output_layer:
            self.curname = name
            return self.CFF
        elif name in self.allCFFs:
            # if asked for CFF which is not in output_layer, return 0
            return self.zero
        else:
            raise AttributeError

    def zero(self, *args, **kwargs):
        return 0

    def CFF(self, pt, outputvalue=None):
        ind = self.output_layer.index(self.curname)
        if self.parameters['outputvalue'] != None:
            # this occurs during training: value is set by training
            # routine by calling with outputvalue set by training routine
            return self.parameters['outputvalue'][ind]
        ar = []
        for net in self.nets:
            #FIXME: network is unneccessarily run for each CFF. This could
            # be optimized by running it just once for each pt.
            ar.append(net.activate([pt.xB, pt.t])[ind])
        all = array(ar).flatten()
        if self.parameters.has_key('nnet'):
            if self.parameters['nnet'] == 'ALL':
                return all
            elif self.parameters['nnet'] == 'AVG':
                return all.mean()
            else: # we want particular net
                try:
                    return all[self.parameters['nnet']]
                except IndexError:
                    raise IndexError, str(self)+' has only '+str(len(self.nets))+' nets!'
        # by default, we get mean value (FIXME:this should never occurr?)
        else:
            return all.mean()

##  --- Complete models built from the above components ---

class ModelDR(ComptonModelDR, ElasticDipole):
    """Complete model as in arXiv:0904.0458.."""


class ModelNN(ComptonNeuralNets, ElasticDipole):
    """Complete model."""


class ModelNNH(ComptonNNH, ElasticDipole):
    """Complete model - devel version with just ImH and ReH."""

