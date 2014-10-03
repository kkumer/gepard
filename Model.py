""" Definitions of models. 

"Model" is a set of functions (typically CFFs or GPDs) which depend on
parameters (some of which can be provided by minimization routine in the
fitting procedure).  Theoretical "approach", when given an instance of a model
and parameter values can calculate observables.

"""
#from IPython.Debugger import Tracer; debug_here = Tracer()
import pickle, sys, logging

from numpy import log, pi, imag, real, sqrt, cos, sin
from numpy import ndarray, array
import scipy.stats
from scipy.special import j0, j1

from quadrature import PVquadrature, bquadrature, rthtquadrature
from utils import flatten, hubDict, stringcolor
from constants import tolerance2, GeVfm, Mp
import dispersion as DR

import pygepard as g1
import pygepard2 as g2
import pygepard3 as g3
import optModel

_lg = logging.getLogger('p.%s' % __name__)
#_lg.setLevel(logging.DEBUG)  #DEBUG, INFO, WARNING, ERROR, CRITICAL

class Model(object):
    """Base class for all models."""

    def __init__(self, **kwargs):
        """ 
        optimization -- use C/Fortran extensions or some such
        
        """
        self.optimization = kwargs.pop('optimization', False)
        # Intially all parameters are fixed and should be released by user
        exec('fixed = {' + ", ".join(map(lambda x: "'fix_%s': %s" % x, 
                    zip(self.parameter_names, len(self.parameter_names)*['True']))) + '}')
        # FIXME: duplication of stuff: parameters and ndparameters!
        self.parameters.update(fixed)
        # right-pad with zeros to the array of 20 elements needed by Fortran
        self.ndparameters = array([self.parameters[name] for name in self.parameter_names]+
                [0. for k in range(20-len(self.parameter_names))])
        #self.res = array([0. for k in range(18)])

    def release_parameters(self, *args):
        """Release parameters for fitting.

        If allowed ranges have to be changed from default, user needs
        to modify parameters_dict dictionary directly.
        FIXME:
        Note that relasing and fixing parameters for model instance is sensible
        only before creating Fitter instance! Afterwards, one has to fix and release
        params both for Fitter and for model.

        """
        for par in args:
            if par not in self.parameter_names:
                raise ValueError('Parameter "%s" is not defined in model %s' 
                        % (par, self))
            self.parameters['fix_'+par] = False

    def fix_parameters(self, *args):
        """Fix parameters so they are not fitting variables.

        FIXME:
        Note that relasing and fixing parameters for model instance is sensible
        only before creating Fitter instance! Afterwards, one has to fix and release
        params both for Fitter and for model.

        """
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

    def print_parameters(self, compare_with=[], exact=False, colors=True):
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
            if exact:
                row = '%5s -> %-g,' % (name, value)
            else:
                row = '%5s -> %-5.3g' % (name, value)
            if self.parameters.has_key('limit_'+name):
                lo, hi = self.parameters['limit_'+name]
                if (abs((lo-value)*(hi-value)) < 0.001):
                    row = stringcolor(row, 'red', colors)
            if self.parameters['fix_'+name] == False:
                row = stringcolor(row, 'green', colors)
            for model in compare_with:
                try:
                    value2 =  model.parameters[name]
                except KeyError:
                    # compared model doesn't have this parameter
                    value2 = 0
                app = '   %-5.3g' % value2
                #app = ('   '+parform) % value2
                # calculate relative diff, or absolute if value is zero
                diff = value - value2
                if value != 0:
                    diff = diff/abs(value)
                if diff > 0.05:
                    app = stringcolor(app, 'red', colors)
                elif diff < -0.05:
                    app = stringcolor(app, 'blue', colors)
                row += app
            row += '\n'
            s += row
        print s

    def print_parameters_fortran(self, output=sys.stdout):
        """Print model parameters in Fortran form
        
        """
        for ind in self.parameters_index:
            name = self.parameters_index[ind]
            value = self.parameters[name]
            output.write('      PAR(%i) = %-g\n' % (ind, value))

    def free_parameters(self):
        """Return just free (non-fixed) parameters."""

        # fitting parameters (not fixed)
        return [p for p in self.parameter_names if self.parameters['fix_'+p] == False]

    def print_parameters_errors(self, pvalues=False, ndof=0):
        """Print fitting parameters and their errors.

        Model must have covariance defined.


        """
        for p in self.free_parameters():
            val = self.parameters[p]
            err = sqrt(tolerance2)*sqrt(self.covariance[p,p])
            if pvalues:
                pval = 2*(1.-scipy.stats.t.cdf(abs(val/err), ndof))
                print '%5s = %8.3f +- %5.3f  (p = %.3g)' % (p, val, err, pval)
            else:
                print '%5s = %8.3f +- %5.3f' % (p, val, err)

    def print_covariance(self, colors=True, correlations=False):
        """Pretty-print covariance matrix

        """
        # fitting parameters (not fixed)
        pars = self.free_parameters()
        header = '      |' + len(pars)*'  %5s ' + '\n------+' + len(pars)*'--------'+'\n'
        sys.stdout.write(header % tuple(pars))
        for prow in pars:
            sys.stdout.write('%5s |' % prow)
            for pcol in pars:
                cov = self.covariance[prow, pcol]
                cor = cov/sqrt(self.covariance[prow, prow])/sqrt(
                        self.covariance[pcol, pcol])
                if correlations:
                    it = ' % 5.3f ' % cor
                else:
                    it = ' % 5.3f ' % cov
                if prow==pcol:
                    it = stringcolor(it, 'green', colors)
                if abs(cor) > 0.99 and prow!=pcol:
                    it = stringcolor(it, 'red', colors)
                elif abs(cor) < 0.01:
                    it = stringcolor(it, 'blue', colors)
                sys.stdout.write(it)
            sys.stdout.write('\n') 

    def print_correlation(self, colors=True):
        """Pretty-print correlation matrix

        """
        self.print_covariance(colors=colors, correlations=True)


    def _diff(self, f, p, pt, h=0.05):
        """Compute derivative of f w.r.t. model parameter p at point pt.
        
        Simple difference is used (f(p+h/2)-f(p-h/2))/h.
        f is string representing appropriate method of self.

        """
        fun = self.__getattribute__(f) 
        mem = self.parameters[p]
        self.parameters[p] = mem+h/2.
        up = fun(pt)
        self.parameters[p] = mem-h/2.
        down = fun(pt)
        self.parameters[p] = mem
        return (up-down)/h

    def uncert(self, f, pt):
        """Compute uncertainty of f from self's covariance dictionary.
        
        covariance dict could be provided by Fitter
        f is string representing appropriate method of self.
        Returns ndarray (mean-std, mean+std) for easy plotting.
        
        """
        fun = self.__getattribute__(f) 
        pars = [p for p in self.parameter_names if self.parameters['fix_'+p] == False]
        var = 0
        dfdp = {}
        for p in pars:
            dfdp[p] = self._diff(f, p, pt, h=sqrt(self.covariance[p,p]))
        for p1 in pars:
            for p2 in pars:
                var += dfdp[p1]*self.covariance[p1,p2]*dfdp[p2]
        mean = self.__getattribute__(f)(pt)
        std = sqrt(var)
        return array([mean-std, mean+std])
        
    def covariance_matrix(self, colors=True, correlations=False):
        """Return covariance matrix

        """
        # fitting parameters (not fixed)
        pars = [p for p in self.parameter_names if self.parameters['fix_'+p] == False]
        m = []
        for prow in pars:
            row = []
            for pcol in pars:
                row.append(self.covariance[prow, pcol])
            m.append(row)
        return array(m)
    

class ElasticFormFactors(Model):
    """Dirac and Pauli elastic form factors F_1 and F_2."""


class ElasticDipole(ElasticFormFactors):
    """Dipole approximation from DM's notebook."""

    def F1(self, t):
        """Dirac elastic proton form factor - dipole parametrization."""
        return (1.41 * (1.26 - t))/((0.71 - t)**2 * (3.53 - t))

    def F2(self, t):
        """Pauli elastic proton form factor - dipole parametrization."""
        return 3.2 / ((0.71 - t)**2 * (3.53 - t))


class ElasticKelly(ElasticFormFactors):
    """Kelly's approximation from DM's notebook."""

    def F1(self, t):
        """Dirac elastic proton form factor - Kelly's parametrization."""
        return ((1 + 0.06815437285120148*t)/(1 - 3.118062557942468*t + 
             1.0338391956016382*t**2 - 0.5031268669574522*t**3) - 
             (0.7931031653189349*(1 - 0.03407718642560074*t)*t)/
              (1 - 3.115222792407001*t + 1.520921000705686*t**2 - 
             0.14999913420898098*t**3))/(1 - 0.28397655354667284*t)

    def F2(self, t):
        """Pauli elastic proton form factor - Kelly's parametrization."""
        return (-((1 + 0.06815437285120148*t)/(1 - 3.118062557942468*t + 
             1.0338391956016382*t**2 - 0.5031268669574522*t**3)) + 
          (2.792847351*(1 - 0.03407718642560074*t))/(1 - 
              3.115222792407001*t + 1.520921000705686*t**2 - 
          0.14999913420898098*t**3))/ (1 - 0.28397655354667284*t)

    def F1n(self, t):
        """Dirac elastic neutron form factor - Kelly's parametrization."""
        return ((-0.4842637275288574*t)/((1 - 1.4084507042253522*t)**2*
          (1 - 0.9345440355820054*t)) + 
          (0.5417644379086957*(1 - 0.6598447281533554*t)*t)/
          (1 - 4.168632789020339*t + 1.9408278987597791*t**2 - 
           1.9100884849907935*t**3))/(1 - 0.2831951622975774*t)

    def F2n(self, t):
        """Pauli elastic neutron form factor - Kelly's parametrization."""
        return ((0.4842637275288574*t)/((1 - 1.4084507042253522*t)**2*
          (1 - 0.9345440355820054*t)) - (1.9130427*(1 - 0.6598447281533554*t))/
          (1 - 4.168632789020339*t + 1.9408278987597791*t**2 - 
          1.9100884849907935*t**3))/(1 - 0.2831951622975774*t)
 

class ComptonFormFactors(Model):
    """Twist-two, no-transversity set of 4 CFFs.

    They are set to be zero here. Actual models are built by subclassing this.

    """

    allCFFs = ['ImH', 'ReH', 'ImE', 'ReE', 'ImHt', 'ReHt', 'ImEt', 'ReEt']
    allCFFsb = ['ImH', 'ReH', 'ImE', 'ReE', 'ImHt', 'ReHt', 'ImEt', 'ReEb']
    allCFFeffs = ['ImHeff', 'ReHeff', 'ImEeff', 'ReEeff', 
                     'ImHteff', 'ReHteff', 'ImEteff', 'ReEteff']
    allGPDs = []


    def print_CFFs(self, pt, format=None):
        """Print values of CFFs at given kinematic point."""
        vals = map(lambda cff: getattr(self, cff)(pt), allCFFs)
        if format == 'mma':
            s = "{" + 8*"%s -> %f, "
            s = s[:-2] + "}"
        else:
            s = 8*"%4s = %5.2f\n"
        print s % flatten(tuple(zip(allCFFs, vals)))


    # Initial definition of all CFFs. All just return zero.
    for name in allCFFs:
        exec('def %s(self, pt): return 0.' % name)

    for name in allCFFeffs:
        exec('def %s(self, pt): return 0.' % name)

    # Define E-bar as xi*E-tilde
    def ReEb(self, pt):
        return (pt.xB/2.)*self.ReEt(pt)

    def is_within_model_kinematics(self, pt):
        return ( (1.5 <= pt.Q2 <= 5.) and 
                 (pt.tm < min(1., pt.Q2/4)) and
                 (1e-3 < pt.xB < 0.5)
               )
    
    
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
        #sys.stderr.write('u = %f, fun = %f' % (u, fun(pt,u)))
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
            pvpi = (res + log((1.+pt.xi)/(1.-pt.xi)) * self.ImEt(pt))/pi
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
    """Model for CFFs as in arXiv:0904.0458. + free pion pole"""

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
              'Mv' : 1.0,     'limit_Mv' : (0.4, 4.),
              'rv' : 0.5,     'limit_rv' : (0., 8.),
              'bv' : 2.2,     'limit_bv' : (0.4, 5.),
               'C' : 7.0,      'limit_C' : (-10., 10.),
              'MC' : 1.3,     'limit_MC' : (0.4, 4.),
             'tNv' : 0.0,                             
             'tMv' : 2.7,    'limit_tMv' : (0.4, 4.),
             'trv' : 6.0,    'limit_trv' : (0., 8.),
             'tbv' : 3.0,    'limit_tbv' : (0.4, 5.),
             'rpi' : 1.0,    'limit_rpi' : (-8, 8.),
             'Mpi' : 1.0,    'limit_Mpi' : (0.4, 4.)   }

        # order matters to fit.MinuitFitter, so it is defined by:
        self.parameter_names = ['NS', 'alS', 'alpS', 'MS', 'rS', 'bS',
                                'Nv', 'alv', 'alpv', 'Mv', 'rv', 'bv',
                                'C', 'MC',
                                'tNv', 'tMv', 'trv', 'tbv', 'rpi', 'Mpi']

        # now do whatever else is necessary
        ComptonFormFactors.__init__(self, **kwargs)

    def ReEt(self, pt):
        """Instead of disp. rel. use pole formula"""
        return self.parameters['rpi'] * 2.16444 / (0.0196 - pt.t) / (1. 
            - pt.t/self.parameters['Mpi']**2)**2 / pt.xi


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


class ComptonModelDRPPsea(ComptonModelDRsea):
    """As DRPP but with NS->Nsea. For combining with Gepard sea"""

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
              'Mv' : 1.0,     'limit_Mv' : (0.4, 4.),
              'rv' : 0.5,     'limit_rv' : (0., 8.),
              'bv' : 2.2,     'limit_bv' : (0.4, 5.),
               'C' : 7.0,      'limit_C' : (-10., 10.),
              'MC' : 1.3,     'limit_MC' : (0.4, 4.),
             'tNv' : 0.0,                             
             'tMv' : 2.7,    'limit_tMv' : (0.4, 4.),
             'trv' : 6.0,    'limit_trv' : (0., 8.),
             'tbv' : 3.0,    'limit_tbv' : (0.4, 5.),
             'rpi' : 1.0,    'limit_rpi' : (-8, 8.),
             'Mpi' : 1.0,    'limit_Mpi' : (0.4, 4.)   }

        # order matters to fit.MinuitFitter, so it is defined by:
        self.parameter_names = ['Nsea', 'alS', 'alpS', 'MS', 'rS', 'bS',
                                'Nv', 'alv', 'alpv', 'Mv', 'rv', 'bv',
                                'C', 'MC',
                                'tNv', 'tMv', 'trv', 'tbv', 'rpi', 'Mpi']

        # now do whatever else is necessary
        ComptonFormFactors.__init__(self, **kwargs)

    def ReEt(self, pt):
        """Instead of disp. rel. use pole formula"""
        return self.parameters['rpi'] * 2.16444 / (0.0196 - pt.t) / (1. 
            - pt.t/self.parameters['Mpi']**2)**2 / pt.xi


class ComptonNeuralNets(Model):
    """Neural network CFFs"""

    # FIXME: this variable should be purged out of the code most likely
    allCFFs = ['ImH', 'ReH', 'ImE', 'ReE', 'ImHt', 'ReHt', 'ImEt', 'ReEt']
    allGPDs = []

    def __init__(self, hidden_layers=[7], output_layer=['ImH', 'ReH'], endpointpower=None, useDR=None):
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
                useDR:  use dispersion relations for some Re(CFF)s.
                         E.g.  useDR = ['ReH', 'ReE', 'ReEt', 'ReHt']
        
        """
        self.architecture = [2] + hidden_layers + [len(output_layer)]
        self.output_layer = output_layer
        self.nets = []
        self.netsC = []  # for subtraction constants if useDR
        self.parameters = {'nnet':0, 'outputvalue':None, 'outputvalueC':None}
        self.parameter_names = ['nnet', 'outputvalue', 'outputvalueC']
        self.endpointpower = endpointpower
        self.useDR = useDR
        # now do whatever else is necessary
        #ComptonFormFactors.__init__(self)

    def __getattr__(self, name):
        """Return appropriate CFF function object."""
        #_lg.debug('NN model called with attr = %s\n' % name)
        # FIXME: I don't understand why I have to use this:
        if name in object.__getattribute__(self, 'output_layer'):
        # and this creates infinite recursion:
        #if name in self.output_layer:
            self.curname = name
            return self.CFF
                #and self.useDR \
        #elif hasattr(self, 'useDR') \
        elif self.__dict__.has_key('useDR') \
                and object.__getattribute__(self, 'useDR') \
                and name in object.__getattribute__(self, 'useDR'):
            self.curname = name
            return self.CFF
        elif name in ComptonFormFactors.allCFFs:
            # if asked for CFF which is not in output_layer, return 0
            self.curname = name
            return self.zero
        elif name in ComptonFormFactors.allCFFeffs:
            # if asked for CFF which is not in output_layer, return 0
            self.curname = name
            return self.zero
        elif name in ['endpointpower', 'optimization', 'useDR']:
            if self.__dict__.has_key(name):
                return self.__dict__[name]
            else:
                return None
        else:
            #syst.stderr.write('Possibly caught exception: AttErr: $s\n' % name)
            raise AttributeError, name

    #def zero(self, *args, **kwargs):
        #return 0

    def subtraction(self, pt):
        return 2.25  # temporary
        if self.parameters['outputvalueC'] != None:
            # this occurs during training: value is set by training
            # routine by calling with outputvalueC set by training routine
            return self.parameters['outputvalueC']
        ar = []
        for netC in self.netsC:
            ar.append(netC.activate([pt.t])[0])
        all = array(ar).flatten()
        if self.parameters.has_key('nnet'):
            if self.parameters['nnet'] == 'ALL':
                return all
            elif self.parameters['nnet'] == 'AVG':
                return all.mean()
            else: # we want particular netC
                try:
                    return all[self.parameters['nnet']]
                except IndexError:
                    raise IndexError, str(self)+' has only '+str(len(self.netsC))+' nets!'
        # by default, we get mean value (FIXME:this should never occurr?)
        else:
            return all.mean()

    def CFF(self, pt, xi=0):
        # FIXME: This function is HEAVILY sub-optimal and non-pythonic!
        #_lg.debug('NN model CFF called as = %s\n' % self.curname)
        if hasattr(self, 'useDR') and self.useDR and self.curname in self.useDR:
            #_lg.debug('Doing DR for CFF: %s\n' % self.curname)
            if self.curname == 'ReH':
                #debug_here()
                return DR.intV(self.ImH, pt) - self.subtraction(pt)
            elif self.curname == 'ReE':
                return DR.intV(self.ImE, pt) + self.subtraction(pt)
            elif self.curname in ['ReHt', 'ReEt']:
                return DR.intA(self.__getattr__('Im'+self.curname[2:]), pt)
            else:
                raise ValueError, 'Only Re(CFF) can be calculated via disp. rel.'
        ind = self.output_layer.index(self.curname)
        if isinstance(xi, ndarray) and len(self.nets)>0:
            # function was called with third argument that is xi nd array
            # and we already have some nets so disp.int. can be calculated
            x = xi
        elif self.parameters['outputvalue'] != None:
            # this occurs during training: value is set by training
            # routine by calling with outputvalue set by training routine
            return self.parameters['outputvalue'][ind]
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

    def zero(self, pt, xi=0):
        # FIXME: This function is HEAVILY sub-optimal and non-pythonic!
        # FIXME: It is also essentially a  copy of CFF!!! (This was
        #        copied in a hurry just to get right array shape.)
        #_lg.debug('NN model CFF called as = %s\n' % self.curname)
        if isinstance(xi, ndarray) and len(self.nets)>0:
            # function was called with third argument that is xi nd array
            # and we already have some nets so disp.int. can be calculated
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
                ar.append(0.)
            all = array(ar).flatten()
            if self.parameters.has_key('nnet'):
                if self.parameters['nnet'] == 'ALL':
                    res.append(all)
                elif self.parameters['nnet'] == 'AVG':
                    res.append(all.mean())
                else: # we want particular net
                    res.append(0)
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

    # To have different Gepard models available we have to
    # use separate modules - otherwise things clash
    gepardPool = [g1, g2, g3]  #  modules to choose from
    #gepardPool = [g1]  #  modules to choose from

    def __init__(self, cutq2=0.0, ansatz='FIT', p=0, scheme='MSBAR', speed=1, q02=4.0, **kwargs):
        _lg.debug('Creating %s.\n' % str(self))
        # initial values of parameters and limits on their values
        self.cutq2 = cutq2
        self.ansatz = ansatz
        self.p = p
        self.scheme = scheme
        self.speed = speed
        self.q02 = q02
        self.kwargs = kwargs
        self.parameters = {
        #  GPD H
               'NS' : 0.15,
             'AL0S' : 1.0,
             'ALPS' : 0.15,
             'M02S' : 1.0,    'limit_M02S' : (0.1, 1.5),
           'DELM2S' : 0.0,
               'PS' : 2.0,
             'SECS' : 0.0,
             'THIS' : 0.0,
             'KAPS' : 0.0,
            'SKEWS' : 0.0,
               'NG' : 0.5,
             'AL0G' : 1.1,
             'ALPG' : 0.15,
             'M02G' : 0.7,    'limit_M02G' : (0.1, 1.5),
           'DELM2G' : 0.0,
               'PG' : 2.0,
             'SECG' : 0.0,
             'THIG' : 0.0,
             'KAPG' : 0.0,
            'SKEWG' : 0.0,
             'DELB' : 0.0,
               'ND' : 1.0,    
             'AL0D' : 0.5,
             'ALPD' : 1.0,
             'M02D' : 1.0,    'limit_M02D' : (0.1, 1.5),
        # GPD E
             'EAL0S' : 1.0,
             'EALPS' : 0.15,
             'EM02S' : 1.0,    'limit_EM02S' : (0.1, 1.5),
           'EDELM2S' : 0.0,
               'EPS' : 2.0,
             'ESECS' : 0.0,
             'ETHIS' : 0.0,
             'EKAPS' : 0.0,
            'ESKEWS' : 0.0,
             'EAL0G' : 1.1,
             'EALPG' : 0.15,
             'EM02G' : 0.7,    'limit_EM02G' : (0.1, 1.5),
           'EDELM2G' : 0.0,
               'EPG' : 2.0,
             'ESECG' : 0.0,
             'ETHIG' : 0.0,
             'EKAPG' : 0.0,
            'ESKEWG' : 0.0   }


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
             37 : 'THIS',
             41 : 'ND',  
             42 : 'AL0D',
             43 : 'ALPD',
             44 : 'M02D',
             47 : 'THIG',
             48 : 'DELB',
            112 : 'EAL0S',
            113 : 'EALPS',
            114 : 'EM02S',
            115 : 'EDELM2S',
            116 : 'EPS',
            117 : 'ESECS',
            119 : 'ESKEWS',
            122 : 'EAL0G',
            123 : 'EALPG',
            124 : 'EM02G',
            125 : 'EDELM2G',
            126 : 'EPG',
            127 : 'ESECG',
            129 : 'ESKEWG',
            137 : 'ETHIS',
            147 : 'ETHIG' }

        # order matters to fit.MinuitFitter, so it is defined by:
        self.parameter_names = [ 
           'NS', 'AL0S', 'ALPS', 'M02S',
           'DELM2S', 'PS', 'SECS', 'THIS', 'KAPS', 'SKEWS',
           'NG', 'AL0G', 'ALPG', 'M02G',
           'DELM2G', 'PG', 'SECG', 'THIG', 'KAPG', 'SKEWG',
           'DELB',
           'ND', 'AL0D', 'ALPD', 'M02D',
           'EAL0S', 'EALPS', 'EM02S',
           'EDELM2S', 'EPS', 'ESECS', 'ETHIS', 'ESKEWS',
           'EAL0G', 'EALPG', 'EM02G',
           'EDELM2G', 'EPG', 'ESECG', 'ETHIG', 'ESKEWG']

        self.allGPDs = ['gpdHtrajQ', 'gpdHtrajG', 'gpdEtrajQ', 'gpdEtrajG',
                        'gpdHzeroQ', 'gpdHzeroG', 'gpdEzeroQ', 'gpdEzeroG',
                        'gpdHskewQ', 'gpdHskewG',
                        'gpdHQb', 'gpdHQbpol', 'gpdHGb']

        if ansatz == 'NSFIT' or ansatz == 'FITBP':
            self.parameters_index.update({
                 31 : 'NU',
                 32 : 'AL0U',
                 33 : 'ALPU',
                 34 : 'M02U',
                 35 : 'DELM2U',
                 36 : 'PU',
                 41 : 'ND',  
                 42 : 'AL0D',
                 43 : 'ALPD',
                 44 : 'M02D',
                 45 : 'DELM2D',
                 46 : 'PD'})
        elif ansatz not in ['FIT', 'FIT14', 'FITEXP', 'EPH', 'EPHEXP', 'EFL', 
                'EFLEXP', 'HOUCHE', 'NSPHOU', 'NSMHOU','TEST']:
            raise ValueError, "Invalid ansatz: %s\n" % ansatz
        
        if ansatz == 'FITEXP':
            self.parameters['ALPS'] = 0.0   # like in smallx.nb
            self.parameters['ALPG'] = 0.0
            self.parameters['M02G'] = 1./4.63 # from J/Psi production
            self.parameters['limit_M02S'] = (0.0, 1.5)
            self.parameters['limit_M02G'] = (0.0, 1.5)

        self._gepardinit(cutq2, ansatz, p, scheme, speed, q02, **kwargs)   # gepard init
        # now do whatever else is necessary
        ComptonFormFactors.__init__(self, **kwargs)

    def _gepardinit(self, cutq2=0.0, ansatz='FIT', p=0, scheme='MSBAR', speed=1, q02=4.0, **kwargs):
        """Initialize gepard part of model."""
        emptyPoolMessage = '''
        Pool of gepard modules is empty. 
        No new gepard models can be created. Restart everything!\n'''
        if kwargs.pop('newgepard', False):
            # Consume one gepard module from the gepardPool
            try:
                self.g = ComptonGepard.gepardPool.pop()
                _lg.debug('Consumed %s from the gepardPool. Leaving %i.\n' 
                    % (self.g.__file__, len(ComptonGepard.gepardPool)))
            except IndexError:
                sys.stderr.write(emptyPoolMessage)
                return -1
        else:
            # Just use the last gepard module, but leave it in the gepardPool
            try:
                self.g = ComptonGepard.gepardPool[-1]
                _lg.debug('Using %s from gepardPool, and leaving it there. Their number is %i.\n' 
                    % (self.g.__file__, len(ComptonGepard.gepardPool)))
            except IndexError:
                sys.stderr.write(emptyPoolMessage)
                return -1
        # this was in Gepard's GEPARD.INI, which is not needed now
        # but look at it for documentation of what parameters below are
        self.g.parint.speed = speed
        self.g.parint.acc = 3
        self.g.parint.p = p
        self.g.parint.nf = 4
        self.g.parint.czero = 1

        self.g.astrong.mu02 = 2.5
        self.g.astrong.asp = array([0.0606, 0.0518, 0.0488])

        self.g.parflt.q02 = q02
        self.g.parflt.rf2 = 1.0
        self.g.parflt.rr2 = 1.0
        self.g.parflt.rdaf2 = 1.0
        self.g.parflt.rgpdf2 = 1.0

        self.g.mbcont.c = 0.35
        self.g.mbcont.phi = 1.57079632
        self.g.mbcont.cnd = -0.25
        self.g.mbcont.phind = 1.57

        self.g.parchr.ansatz = array([c for c in ansatz + (6-len(ansatz))*' ']) # array(6)
        self.g.parchr.scheme = array([c for c in scheme + (5-len(scheme))*' ']) # array(5)

        self.g.init()
        # Cutting-off evolution  at Q2 = cutq2
        # Evaluate evolved C at this scale now.
        # FIXME: check that this works 
        for pid in range(-4,4):
            self.g.nqs.nqs[pid] = 1
            self.g.qs.qs[pid, 0] = self.cutq2
            self.g.kinematics.q2 = self.cutq2
            self.g.evolc(1, 1)
            self.g.evolc(2, 1)
            self.g.evolc(3, 1)

        self.g.newcall = 1
        # number of points on MB contour
        self.g.npts = 2**self.g.parint.acc * 12 / self.g.parint.speed
        # now do whatever else is necessary

    def __getstate__(self):
        # We have to remove unpicklable gepard module object
        _lg.debug('Shelving [ComptonGepard] %s.\n' % str(self))
        del self.g
        if self.kwargs.has_key('newgepard'):
            _lg.warning("Model with newgepard atribute saved. It will consume one GepardPool module when restored.\n")
        return self.__dict__

    def __setstate__(self, dict):
        _lg.debug('Unshelving %s.' % str(self))
        self.__dict__ = dict
        # We now have to reconstruct gepard module object
        self._gepardinit(cutqq2=self.cutq2, ansatz=self.ansatz,
                p=self.p, scheme=self.scheme, speed=self.speed, q02=self.q02, **self.kwargs)
        # FIXME: Next is for compatibility with old models saved in database.
        #        Should upgrade database and remoe this:
        #self._gepardinit(cutqq2=self.cutq2, ansatz=self.ansatz, speed=self.speed, 
                #q02=self.q02, **self.kwargs)
        #if hasattr(self, 'p'):
            #self.g.parint.p = self.p
        #if hasattr(self, 'scheme'):
            #self.g.parchr.scheme = self.scheme


    def return_gepard(self):
        _lg.debug('Returning %s to GepardPool.' % (self.g.__file__,))
        ComptonGepard.gepardPool.append(self.g)


    def _GepardFFs(self, pt, FF='cfff'):
        """Call gepard routine that calculates CFFs or TFFs."""
        for i in self.parameters_index:
            self.g.par.par[i-1] = self.parameters[self.parameters_index[i]]

        self.g.kinematics.q2 = pt.Q2
        self.g.kinematics.xi = pt.xi
        self.g.kinematics.del2 = pt.t

        self.g.mt.nmts = 1
        self.g.mt.mtind = 0 
        self.g.mts.mts[0] = - pt.t

        getattr(self.g, FF)()
        self.g.newcall = 0

    def DISF2(self, pt):
        """Call gepard routine that calculates DIS F2."""
        for i in self.parameters_index:
            self.g.par.par[i-1] = self.parameters[self.parameters_index[i]]
        self.g.parint.pid = 0
        self.g.kinematics.q2 = pt.Q2
        # Note a hack in Fortran gepard where for F2 calculation
        # XI should actually be xB, and not xB/2 !
        self.g.kinematics.xi = 2*pt.xi/(1.+pt.xi)
        self.g.kinematics.del2 = 0
        self.g.f2f()
        return self.g.f2.f2[self.g.parint.p]


    def gpdHtrajQ(self, pt):
        """GPD H^q on xi=x trajectory.  """
        self.g.parint.pid = -3
        self.g.newcall = 1
        return self.ImH(pt)/pi

    def gpdHtrajG(self, pt):
        """GPD H^g on xi=x trajectory."""
        self.g.parint.pid = -4
        self.g.newcall = 1
        return pt.xi*self.ImH(pt)/pi

    def gpdEtrajQ(self, pt):
        """GPD E on xi=x trajectory."""
        self.g.parint.pid = -3
        self.g.newcall = 1
        return self.ImE(pt)/pi

    def gpdEtrajG(self, pt):
        """GPD H^g on xi=x trajectory."""
        self.g.parint.pid = -4
        self.g.newcall = 1
        return pt.xi*self.ImE(pt)/pi

    def gpdHzeroQ(self, pt, ts=None):
        """GPD H^sea on xi=0 trajectory.
        FIXME: After this, calling self.ImH is broken!!
        FIXME: LOT of code duplication here!!
        """
        memsub = {'SECS': self.parameters['SECS'],    'SECG': self.parameters['SECG'], 
                'ESECS' : self.parameters['ESECS'], 'ESECG' : self.parameters['ESECG'],
                  'THIS': self.parameters['THIS'],    'THIG': self.parameters['THIG'],
                'ETHIS' : self.parameters['ETHIS'], 'ETHIG' : self.parameters['ETHIG']}
        zerosub = {'SECS': 0, 'SECG': 0, 'ESECS' : 0, 'ESECG' : 0,
                   'THIS': 0, 'THIG': 0, 'ETHIS' : 0, 'ETHIG' : 0}
        self.parameters.update(zerosub)
        self.g.parint.pid = -1
        if isinstance(ts, ndarray):
            tmem = pt.t
            res = []
            for t in ts:
                pt.t = t
                self.g.newcall = 1
                res.append(self.ImH(pt))
            self.parameters.update(memsub)
            pt.t = tmem
            return array(res)/pi
        else:
            self.g.newcall = 1
            res = self.ImH(pt)
            self.parameters.update(memsub)
            return res/pi

    def gpdHzeroG(self, pt, ts=None):
        """GPD H^g on xi=0 trajectory."""
        memsub = {'SECS': self.parameters['SECS'],    'SECG': self.parameters['SECG'], 
                'ESECS' : self.parameters['ESECS'], 'ESECG' : self.parameters['ESECG'],
                  'THIS': self.parameters['THIS'],    'THIG': self.parameters['THIG'],
                'ETHIS' : self.parameters['ETHIS'], 'ETHIG' : self.parameters['ETHIG']}
        zerosub = {'SECS': 0, 'SECG': 0, 'ESECS' : 0, 'ESECG' : 0,
                   'THIS': 0, 'THIG': 0, 'ETHIS' : 0, 'ETHIG' : 0}
        self.parameters.update(zerosub)
        self.g.parint.pid = -2
        if isinstance(ts, ndarray):
            tmem = pt.t
            res = []
            for t in ts:
                pt.t = t
                self.g.newcall = 1
                res.append(self.ImH(pt))
            self.parameters.update(memsub)
            pt.t = tmem
            return pt.xi*array(res)/pi
        else:
            self.g.newcall = 1
            res = self.ImH(pt)
            self.parameters.update(memsub)
            return pt.xi*res/pi

    def gpdEzeroQ(self, pt, ts=None):
        """GPD E^sea on xi=0 trajectory."""
        memsub = {'SECS': self.parameters['SECS'],    'SECG': self.parameters['SECG'], 
                'ESECS' : self.parameters['ESECS'], 'ESECG' : self.parameters['ESECG'],
                  'THIS': self.parameters['THIS'],    'THIG': self.parameters['THIG'],
                'ETHIS' : self.parameters['ETHIS'], 'ETHIG' : self.parameters['ETHIG']}
        zerosub = {'SECS': 0, 'SECG': 0, 'ESECS' : 0, 'ESECG' : 0,
                   'THIS': 0, 'THIG': 0, 'ETHIS' : 0, 'ETHIG' : 0}
        self.parameters.update(zerosub)
        self.g.parint.pid = -1
        if isinstance(ts, ndarray):
            tmem = pt.t
            res = []
            for t in ts:
                pt.t = t
                self.g.newcall = 1
                res.append(self.ImE(pt))
            self.parameters.update(memsub)
            pt.t = tmem
            return array(res)/pi
        else:
            self.g.newcall = 1
            res = self.ImE(pt)
            self.parameters.update(memsub)
            return res/pi

    def gpdEzeroG(self, pt, ts=None):
        """GPD E^g on xi=0 trajectory."""
        memsub = {'SECS': self.parameters['SECS'],    'SECG': self.parameters['SECG'], 
                'ESECS' : self.parameters['ESECS'], 'ESECG' : self.parameters['ESECG'],
                  'THIS': self.parameters['THIS'],    'THIG': self.parameters['THIG'],
                'ETHIS' : self.parameters['ETHIS'], 'ETHIG' : self.parameters['ETHIG']}
        zerosub = {'SECS': 0, 'SECG': 0, 'ESECS' : 0, 'ESECG' : 0,
                   'THIS': 0, 'THIG': 0, 'ETHIS' : 0, 'ETHIG' : 0}
        self.parameters.update(zerosub)
        self.g.parint.pid = -2
        if isinstance(ts, ndarray):
            tmem = pt.t
            res = []
            for t in ts:
                pt.t = t
                self.g.newcall = 1
                res.append(self.ImE(pt))
            self.parameters.update(memsub)
            pt.t = tmem
            return pt.xi*array(res)/pi
        else:
            self.g.newcall = 1
            res = self.ImE(pt)
            self.parameters.update(memsub)
            return pt.xi*res/pi

    def gpdHskewQ(self, pt):
        """Skewneess of GPD H_Q"""
        return self.gpdHtrajQ(pt)/self.gpdHzeroQ(pt)

    def gpdHskewG(self, pt):
        """Skewneess of GPD H_G"""
        return self.gpdHtrajG(pt)/self.gpdHzeroG(pt)



    def gpdHQb(self, pt):
        """GPD H^sea in b space in fm^-2."""
        b = sqrt(pt.bx**2 + pt.by**2)
        return bquadrature(lambda d: d*j0(b*d/GeVfm)*self.gpdHzeroQ(pt, ts=-d**2), 
                0.0, 2.5) / (2*pi*GeVfm**2)

    def gpdHQbpol(self, pt):
        """polarized GPD H^sea in b space in fm^-2."""
        b = sqrt(pt.bx**2 + pt.by**2)
        aux = bquadrature(lambda d: d**2*j1(b*d/GeVfm)*self.gpdEzeroQ(pt, ts=-d**2), 
                0.0, 2.5) / (2*pi*GeVfm**2)
        return self.gpdHQb(pt) + pt.bx*aux/(2*b*Mp)

    def gpdHQbpolpol(self, pt, r, tht):
        """Same as gpdHQbpol but in polar coordinate system.
        FIXME: some ugly coding here.
        """
        if not isinstance(r, ndarray):
            assert isinstance(r, (int, long, float)) 
            ra = array([r])
        else:
            ra = r
        res = []
        for rr in ra:
            pt.bx = rr*cos(tht)
            pt.by = rr*sin(tht)
            b = sqrt(pt.bx**2 + pt.by**2)
            aux = bquadrature(lambda d: d**2*j1(b*d/GeVfm)*self.gpdEzeroQ(pt, ts=-d**2), 
                    0.0, 2.5) / (2*pi*GeVfm**2)
            res.append(self.gpdHQb(pt) + pt.bx*aux/(2*b*Mp))
        if not isinstance(r, ndarray):
            r = ra[0]
        else:
            r = ra
        return array(res)

    def rth(self, pt, tht):
        """Calculate <r(tht)> """
        norm = rthtquadrature(lambda r: r*self.gpdHQbpolpol(pt, r, tht), 0.0, 2.5)
        prob = rthtquadrature(lambda r: r*r*self.gpdHQbpolpol(pt, r, tht), 0.0, 2.5)
        return prob/norm

    def gpdHQbIntg(self, pt, ba):
        """Same as gpdHQb but convenient as integrand."""
        res = []
        for b in ba:
            aux = bquadrature(lambda d: d*j0(b*d/GeVfm)*self.gpdHzeroQ(pt, ts=-d**2), 
                0.0, 2.5) / (2*pi*GeVfm**2)
            res.append(aux)
        return array(res)

    def bsq(self, pt):
        """Calculate unpolarized <b^2> for model m."""
        norm = rthtquadrature(lambda b: b*self.gpdHQbIntg(pt, b), 0.0, 2.5)
        aux = rthtquadrature(lambda b: b**3*self.gpdHQbIntg(pt, b), 0.0, 2.5)
        return aux/norm
    
    def gpdHGb(self, pt):
        """GPD H^G in b space in fm^-2."""
        b = sqrt(pt.bx**2 + pt.by**2)
        return bquadrature(lambda d: d*j0(b*d/GeVfm)*self.gpdHzeroG(pt, ts=-d**2), 
                0.0, 2.5) / (2*pi*GeVfm**2)

    def gpdHGbpol(self, pt):
        """polarized GPD H^sea in b space in fm^-2."""
        b = sqrt(pt.bx**2 + pt.by**2)
        aux = bquadrature(lambda d: d**2*j1(b*d/GeVfm)*self.gpdEzeroG(pt, ts=-d**2), 
                0.0, 2.5) / (2*pi*GeVfm**2)
        return self.gpdHGb(pt) + pt.by*aux/(2*b*Mp)


    def beff(self, pt):
        """slope of  CFFH."""
        tmem = pt.t
        pt.t = -0.7
        self.g.newcall = 1
        up = self.ImH(pt)
        pt.t = -0.1
        self.g.newcall = 1
        down = self.ImH(pt)
        pt.t = tmem
        return - 2 * log(up/down)/0.7


    def ImH(self, pt):
        """Imaginary part of CFF H."""
        if self.g.newcall:
            self._GepardFFs(pt)
        return imag(self.g.cff.cff[self.g.parint.p])

    def ReH(self, pt):
        """Real part of CFF H."""
        if self.g.newcall:
            self._GepardFFs(pt)
        return real(self.g.cff.cff[self.g.parint.p])

    def ImHrho(self, pt):
        """Imaginary part of TFF H_rho."""
        if self.g.newcall:
            self._GepardFFs(pt, 'tfff')
        return imag(self.g.tff.tffh[self.g.parint.p])

    def ReHrho(self, pt):
        """Real part of TFF H_rho."""
        if self.g.newcall:
            self._GepardFFs(pt, 'tfff')
        return real(self.g.tff.tffh[self.g.parint.p])

    def ImE(self, pt):
        """Imaginary part of CFF E."""
        if self.g.newcall:
            self._GepardFFs(pt)
        return imag(self.g.cff.cffe[self.g.parint.p])

    def ReE(self, pt):
        """Real part of CFF E."""
        if self.g.newcall:
            self._GepardFFs(pt)
        return real(self.g.cff.cffe[self.g.parint.p])


class ComptonHybrid(ComptonFormFactors):
    """This combines gepard for small xB and DR model for valence xB."""

    def __init__(self, instGepard, instDR, **kwargs):
        self.Gepard = instGepard  # instance of ComptonGepard
        self.DR = instDR  # instance of ComptonModelDR
        self.DR.parameters['Nsea'] = 0.  # sea comes from Gepard part
        self.DR.ndparameters[0] = 0.  # sea comes from Gepard part
        #self.parameters = hubDict(self.Gepard.parameters, self.DR.parameters)
        #self.parameter_names = self.Gepard.parameter_names + self.DR.parameter_names
        self.parameters = hubDict(self.DR.parameters, self.Gepard.parameters)
        self.parameter_names = self.DR.parameter_names + self.Gepard.parameter_names

        self.g = self.Gepard.g
        # now do whatever else is necessary
        ComptonFormFactors.__init__(self, **kwargs)

    def __getstate__(self):
        # We have to remove unpicklable gepard module object
        _lg.debug('Shelving [ComptonHybrid] %s.\n' % str(self))
        del self.g
        return self.__dict__

    def __setstate__(self, dict):
        _lg.debug('Unshelving %s.' % str(self))
        self.__dict__ = dict
        self.g = self.Gepard.g

    def close(self):
        self.Gepard.close()

    def is_within_model_kinematics(self, pt):
        # relaxing xBmin and removing Q2max
        return ( (1.5 <= pt.Q2) and 
                 (pt.tm < min(1., pt.Q2/4)) and
                 (1e-5 < pt.xB < 0.5)
               )

    def ImH(self, pt, xi=0):
        return  self.Gepard.ImH(pt) + self.DR.ImH(pt, xi)

    def ReH(self, pt):
        return  self.Gepard.ReH(pt) + self.DR.ReH(pt)

    def ImE(self, pt, xi=0):
        return  self.Gepard.ImE(pt) + self.DR.ImE(pt, xi)

    def ReE(self, pt):
        return  self.Gepard.ReE(pt) + self.DR.ReE(pt)

    # FIXME: tildes are not provided by Gepard

    def ImHt(self, pt, xi=0):
        return  self.DR.ImHt(pt, xi)

    def ReHt(self, pt):
        return  self.DR.ReHt(pt)

    def ImEt(self, pt):
        return  self.DR.ImEt(pt)

    def ReEt(self, pt):
        return  self.DR.ReEt(pt)

    # FIXME: rho production TFFs and DIS F2 are not provided by DR

    def ImHrho(self, pt):
        return  self.Gepard.ImHrho(pt)

    def ReHrho(self, pt):
        return  self.Gepard.ReHrho(pt)

    def DISF2(self, pt):
        return  self.Gepard.DISF2(pt)


class ComptonLocal(ComptonFormFactors):
    """For local fitting of CFFs which are themselves fit paramters."""

    def __init__(self, **kwargs):
        self.parameters = {
            'pImH' :  0.0,                                 
            'pReH' :  0.0,                                 
            'pImE' :  0.0,                                 
            'pReE' :  0.0,                                 
           'pImHt' :  0.0,                                 
           'pReHt' :  0.0,                                 
           'pImEt' :  0.0,                                 
           'pReEt' :  0.0}
        self.parameter_names = ['pImH', 'pReH', 'pImHt', 'pReHt', 'pImE', 'pReE',
           'pImEt', 'pReEt']

        # now do whatever else is necessary
        ComptonFormFactors.__init__(self, **kwargs)

    # Each of 8 CFFs just returns value of pCFF parameter
    for name in ComptonFormFactors.allCFFs:
        exec('def %s(self, pt): return self.parameters["p%s"]' % (name, name))

##  --- Complete models built from the above components ---

class ModelDR(ComptonModelDR, ElasticDipole):
    """Complete model as in arXiv:0904.0458.."""

class Gepard(ComptonGepard, ElasticKelly):
    """Complete model as in arXiv:0904.0458.."""


class ModelDRKelly(ComptonModelDR, ElasticKelly):
    """Same, but with Kelly elastic form factors."""


class ModelNN(ComptonNeuralNets, ElasticDipole):
    """Complete model."""


class HybridDipole(ComptonHybrid, ElasticDipole):
    """Complete hybrid model."""


class ModelDRPP(ComptonModelDRPP, ElasticDipole):
    """Complete model as in arXiv:0904.0458. + free pion pole normalization."""


class HybridKelly(ComptonHybrid, ElasticKelly):
    """As Hybrid, but with Kelly elasticd FFs."""

class ModelLocal(ComptonLocal, ElasticKelly):
    """Model for local fitting of CFFs themselves."""

    def ImEt(self, pt):
        #To get sensible numbers for parameter.
        return (2.-pt.xB)/pt.xB * self.parameters['pImEt']

    def ReEt(self, pt):
        #To get sensible numbers for parameter.
        return (2.-pt.xB)/pt.xB * self.parameters['pReEt']
