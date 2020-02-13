""" Definitions of models. 

"Model" is a set of functions (typically CFFs or GPDs) which depend on
parameters (some of which can be provided by minimization routine in the
fitting procedure).  Theoretical "approach", when given an instance of a model
and parameter values can calculate observables.

"""


#from IPython.core.debugger import set_trace
import pickle, sys, logging

from numpy import log, pi, imag, real, sqrt, cos, sin, exp
from numpy import ndarray, array, sum
import scipy.stats
from scipy.special import j0, j1, gamma, beta

from quadrature import PVquadrature, bquadrature, rthtquadrature
from utils import flatten, hubDictNew, stringcolor
from constants import tolerance2, GeVfm, Mp
import dispersion as DR

import pygepard as g1
import pygepard2 as g2
import pygepard3 as g3
import pygepard4 as g4
import pygepard5 as g5
import pygepard6 as g6
import pygepard7 as g7
import optModel

_lg = logging.getLogger('p.%s' % __name__)
_lg.setLevel(logging.WARNING)  #DEBUG, INFO, WARNING, ERROR, CRITICAL
#_lg.setLevel(logging.DEBUG)  #DEBUG, INFO, WARNING, ERROR, CRITICAL

class Model(object):
    """Base class for all models."""


    def __init__(self, **kwargs):
        """ 
        optimization -- use C/Fortran extensions or some such
        
        """
        self.optimization = kwargs.pop('optimization', False)
        # Intially all parameters are fixed and should be released by user
        fixed = {'fix_{}'.format(p):True for p in self.parameter_names}
        # FIXME: duplication of stuff: parameters and ndparameters!
        self.parameters.update(fixed)
        # right-pad with zeros to the array of 50 elements 
        self.ndparameters = array([self.parameters[name] for name in self.parameter_names]+
                [0. for k in range(50-len(self.parameter_names))])
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

    def print_parameters(self, compare_with=[], exact=False, colors=True, delta=0.05):
        """Pretty-print parameters and their values.

        Variable parameters are printed green, while parameters with values
        at the limits of their range are printed red.
        If additional models are given in compare_with list, their parameter
        values are also printed and differences larger than delta are denoted
        by blue and red coloring.

        """
        s = ""
        for name in self.parameter_names:
            value = self.parameters[name]
            if exact:
                row = '%5s -> %-g,' % (name, value)
            else:
                row = '%5s -> %-5.3g' % (name, value)
            if 'limit_'+name in self.parameters:
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
                if exact:
                    app = '   %-g' % value2
                else:
                    app = '   %-5.3g' % value2
                #app = ('   '+parform) % value2
                # calculate relative diff, or absolute if value is zero
                diff = value - value2
                if value != 0:
                    diff = diff/abs(value)
                if diff > delta:
                    app = stringcolor(app, 'red', colors)
                elif diff < -delta:
                    app = stringcolor(app, 'blue', colors)
                row += app
            row += '\n'
            s += row
        print(s)

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
                print('%5s = %8.3f +- %5.3f  (p = %.3g)' % (p, val, err, pval))
            else:
                print('%5s = %8.3f +- %5.3f' % (p, val, err))

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

    def nF1(self, t):
        """Dirac elastic neutron form factor - Kelly's parametrization."""
        return ((-0.4842637275288574*t)/((1 - 1.4084507042253522*t)**2*
          (1 - 0.9345440355820054*t)) + 
          (0.5417644379086957*(1 - 0.6598447281533554*t)*t)/
          (1 - 4.168632789020339*t + 1.9408278987597791*t**2 - 
           1.9100884849907935*t**3))/(1 - 0.2831951622975774*t)

    def nF2(self, t):
        """Pauli elastic neutron form factor - Kelly's parametrization."""
        return ((0.4842637275288574*t)/((1 - 1.4084507042253522*t)**2*
          (1 - 0.9345440355820054*t)) - (1.9130427*(1 - 0.6598447281533554*t))/
          (1 - 4.168632789020339*t + 1.9408278987597791*t**2 - 
          1.9100884849907935*t**3))/(1 - 0.2831951622975774*t)
 

class ElasticZero(ElasticFormFactors):
    """Set F1=F2=0 to get just DVCS^2."""

    def F1(self, t):
        return 0.

    def F2(self, t):
        return 0.



class ComptonFormFactors(Model):
    """Twist-two, no-transversity set of 4 CFFs.

    They are set to be zero here. Actual models are built by subclassing this.

    """

    
    parameters = {}
    parameter_names = []
    allCFFs = ['ImH', 'ReH', 'ImE', 'ReE', 'ImHt', 'ReHt', 'ImEt', 'ReEt']
    allnCFFs = ['nImH', 'nReH', 'nImE', 'nReE', 'nImHt', 'nReHt', 'nImEt', 'nReEt'] # neutron
    allCFFsb = ['ImH', 'ReH', 'ImE', 'ReE', 'ImHt', 'ReHt', 'ImEt', 'ReEb']
    allCFFeffs = ['ImHeff', 'ReHeff', 'ImEeff', 'ReEeff', 
                     'ImHteff', 'ReHteff', 'ImEteff', 'ReEteff']
    allGPDs = []


    def print_CFFs(self, pt, format=None):
        """Print values of CFFs at given kinematic point."""
        vals = [getattr(self, cff)(pt) for cff in self.allCFFs]
        if format == 'mma':
            s = "{" + 8*"%s -> %f, "
            s = s[:-2] + "}"
        else:
            s = 8*"%4s = %5.2f\n"
        print(s % flatten(tuple(zip(self.allCFFs, vals))))


    # Initial definition of all CFFs. All just return zero.
    for name in allCFFs:
        exec('def %s(self, pt): return 0.' % name)

    for name in allnCFFs:
        exec('def %s(self, pt): return 0.' % name)

    for name in allCFFeffs:
        exec('def %s(self, pt): return 0.' % name)

    # Define E-bar as xi*E-tilde
    def ReEb(self, pt):
        return (pt.xi)*self.ReEt(pt)

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

    @staticmethod
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
             'tal' : 0.43,                             
             'talp' : 0.85,                             
             'tMv' : 2.7,    'limit_tMv' : (0.4, 2.),
             'trv' : 6.0,    'limit_trv' : (0., 8.),
             'tbv' : 3.0,    'limit_tbv' : (0.4, 5.)   }

        # order matters to fit.MinuitFitter, so it is defined by:
        self.parameter_names = ['NS', 'alS', 'alpS', 'MS', 'rS', 'bS',
                                'Nv', 'alv', 'alpv', 'Mv', 'rv', 'bv',
                                'C', 'MC',
                                'tNv', 'tal', 'talp',
                                'tMv', 'trv', 'tbv']

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
        try:
            regge = (-p['tal']-p['talp']*t)
        except KeyError:
            # Old models take Regge trajectory params from H:
            regge = (-p['alv']-p['alpv']*t)
        val = ( (2.*4./9. + 1./9.) * p['tNv'] * p['trv'] * 
            twox**regge *
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
             'tal' : 0.43,                             
             'talp' : 0.85,                             
             'tMv' : 2.7,    'limit_tMv' : (0.4, 4.),
             'trv' : 6.0,    'limit_trv' : (0., 8.),
             'tbv' : 3.0,    'limit_tbv' : (0.4, 5.),
             'rpi' : 1.0,    'limit_rpi' : (-8, 8.),
             'Mpi' : 1.0,    'limit_Mpi' : (0.4, 4.)   }

        # order matters to fit.MinuitFitter, so it is defined by:
        self.parameter_names = ['NS', 'alS', 'alpS', 'MS', 'rS', 'bS',
                                'Nv', 'alv', 'alpv', 'Mv', 'rv', 'bv',
                                'C', 'MC',
                                'tNv', 'tal', 'talp',
                                'tMv', 'trv', 'tbv', 'rpi', 'Mpi']

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
             'tal' : 0.43,                             
             'talp' : 0.85,                             
             'tMv' : 2.7,    'limit_tMv' : (0.4, 2.),
             'trv' : 6.0,    'limit_trv' : (0., 8.),
             'tbv' : 3.0,    'limit_tbv' : (0.4, 5.)   }

        # order matters to fit.MinuitFitter, so it is defined by:
        self.parameter_names = ['Nsea', 'alS', 'alpS', 'MS', 'rS', 'bS',
                                'Nv', 'alv', 'alpv', 'Mv', 'rv', 'bv',
                                'C', 'MC',
                                'tNv', 'tal', 'talp',
                                'tMv', 'trv', 'tbv']

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

    def nImH(self, pt, xi=0):
        """Imaginary part of neutron CFF H."""
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
        # Just isospin rotation from proton valence ImH:
        val = ( (1.*4./9. + 2./9.) * p['Nv'] * p['rv'] * twox**(-p['alv']-p['alpv']*t) *
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
        try:
            regge = (-p['tal']-p['talp']*t)
        except KeyError:
            # Old models take Regge trajectory params from H:
            regge = (-p['alv']-p['alpv']*t)
        val = ( (2.*4./9. + 1./9.) * p['tNv'] * p['trv'] * 
            twox**regge *
                 onex**p['tbv'] / (1. - onex*t/(p['tMv']**2))  )
        return pi * val / (1.+x)

    def nImHt(self, pt, xi=0):
        """Imaginary part of neutron CFF Ht i.e. \tilde{H}."""
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
        try:
            regge = (-p['tal']-p['talp']*t)
        except KeyError:
            # Old models take Regge trajectory params from H:
            regge = (-p['alv']-p['alpv']*t)
        # Just isospin rotation from proton valence ImH:
        val = ( (1.*4./9. + 2./9.) * p['tNv'] * p['trv'] * 
            twox**regge *
                 onex**p['tbv'] / (1. - onex*t/(p['tMv']**2))  )
        return pi * val / (1.+x)


    def ImE(self, pt, xi=0):
        """Imaginary part of CFF E."""
        # Just changing function signature w.r.t. ComptonFormFactors
        # to make it compatible for dispersion integral
        return 0

    nImE = ImE
    # def nImE(self, pt, xi=0):
        # return 0*self.ImH(pt,xi)

    def ReEt(self, pt):
        """Instead of disp. rel. use pole formula."""
        return (2.2390424 * (1. - (1.7*(0.0196 - pt.t))/(1. 
            - pt.t/2.)**2))/((0.0196 - pt.t)*pt.xi)

    nReEt = ReEt     # take same pion pole for neutron

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
             'tal' : 0.43,                             
             'talp' : 0.85,                             
             'tMv' : 2.7,    'limit_tMv' : (0.4, 4.),
             'trv' : 6.0,    'limit_trv' : (0., 8.),
             'tbv' : 3.0,    'limit_tbv' : (0.4, 5.),
             'rpi' : 1.0,    'limit_rpi' : (-8, 8.),
             'Mpi' : 1.0,    'limit_Mpi' : (0.4, 4.)   }

        # order matters to fit.MinuitFitter, so it is defined by:
        self.parameter_names = ['Nsea', 'alS', 'alpS', 'MS', 'rS', 'bS',
                                'Nv', 'alv', 'alpv', 'Mv', 'rv', 'bv',
                                'C', 'MC',
                                'tNv', 'tal', 'talp',
                                'tMv', 'trv', 'tbv', 'rpi', 'Mpi']

        # now do whatever else is necessary
        ComptonFormFactors.__init__(self, **kwargs)

    def ReEt(self, pt):
        """Instead of disp. rel. use pole formula"""
        return self.parameters['rpi'] * 2.16444 / (0.0196 - pt.t) / (1. 
            - pt.t/self.parameters['Mpi']**2)**2 / pt.xi

    nReEt = ReEt     # take same pion pole for neutron

class GK12(ComptonDispersionRelations):
    """ Goloskokov-Kroll PDF model 

        From [arXiv:1210.6975], Kroll:2012sm
        Real part of CFFs is obtained using dispersion relations.

    """


    # We keep the forward-case code just in case it's useful

    def G(self, x, Q2):
        """GK12 model for gluon PDF."""
        Q02 = 4.
        L = log(Q2/Q02)
        cg = array([2.23+0.362*L, 5.43-7.00*L, -34.0+22.5*L, 40.6-21.6*L])
        delg = (1.10+0.06*L-0.0027*L**2) - 1
        ng = 2
        xs = array([x**(j/2.) for j in range(4)])
        return x**(-delg)*(1.-x)**(2*ng+1)*sum(cg*xs)

    def s(self, x, Q2):
        """GK12 model for strange PDF."""
        Q02 = 4.
        L = log(Q2/Q02)
        cs = array([0.123+0.0003*L, -0.327-0.004*L, 0.692-0.068*L, -0.486+0.038*L])
        delg = (1.10+0.06*L-0.0027*L**2) - 1
        ng = 2
        xs = array([x**(j/2.) for j in range(4)])
        return x**(-delg)*(1.-x)**(2*ng+1)*sum(cs*xs)

    def uval(self, x, Q2):
        """GK12 model for u_val PDF."""
        Q02 = 4.
        L = log(Q2/Q02)
        cuval = array([1.52+0.248*L, 2.88-0.940*L, -0.095*L, 0 ])
        delval = 0.48 - 1
        nval = 1
        xs = array([x**(j/2.) for j in range(4)])
        return x**(-delval)*(1.-x)**(2*nval+1)*sum(cuval*xs)

    def dval(self, x, Q2):
        """GK12 model for d_val PDF."""
        Q02 = 4.
        L = log(Q2/Q02)
        cdval = array([0.76+0.248*L, 3.11-1.36*L, -3.99+1.15*L, 0])
        delval = 0.48 - 1
        nval = 1
        xs = array([x**(j/2.) for j in range(4)])
        return x**(-delval)*(1.-x)**(2*nval+1)*sum(cdval*xs)
        
    def udsea(self, x, Q2):
        """GK12 model for u or d sea."""
        kaps = 1.+0.68/(1.+0.52*log(Q2/4.))
        return kaps * self.s(x,Q2)


    def SIG(self, x, Q2):
        """GK12 model for singlet quark."""
        return self.uval(x, Q2) + self.dval(x,Q2) + 4*self.udsea(x,Q2) + 2*self.s(x,Q2)

    def DISF2(self, pt):
        """DIS F2 structure function in GK12 model."""
        x = pt.xB
        Q2 = pt.Q2
        return ( (4./9.) * (self.uval(x, Q2) + 2*self.udsea(x,Q2)) 
               + (1./9.) * (self.dval(x, Q2) + 2*self.udsea(x,Q2) 
                                             + 2*self.s(x,Q2))  )

    def _intval(self, x, eta, alt, j, zero=False):
        """Analytic integral over DD for valence sector."""
        x1 = (x+eta)/(1+eta)
        p = - alt + j/2
        if zero:
            x2 = 0
        else:
            x2 = (x-eta)/(1-eta)
        prefac = (3./2)*gamma(1+p)/gamma(4+p)
        if not zero and (eta/x < 1e-4):
            # use Taylor up to O(eta)
            b = (-2./3)*(p**3+6*p**2+11*p+6)*(x-1)**3*x**p
        else:
            # exact (but instable for small eta)
            b = ((eta**2-x)*(x1**(2+p)-x2**(2+p))
             +(2+p)*eta*(1-x)*(x1**(2+p)+x2**(2+p)))/eta**3
        return prefac * b

    def _intsea(self, x, eta, alt, j, zero=False):
        """Analytic integral over DD for sea sector."""
        x1 = (x+eta)/(1+eta)
        p = - alt + j/2
        if zero:
            x2 = 0
        else:
            x2 = (x-eta)/(1-eta)
        if not zero and (eta/x < 1e-2):
            # use Taylor up to O(eta^2)
            b = eta**2*(p-1)*p+x*(eta**2*(p+6)*(p*(x-2)+7*x)+14*x)
            return (1-x)**5 * b * x**(p-2) / 14.
        else:
            # exact (but instable for small eta)
            prefac = (15./2)*gamma(1+p)/gamma(6+p)/eta**5
            brm = 3*(eta**2-x)**2+eta**2*(p**2+6*p+8)*(1-x)**2
            brp = 3*eta*(eta**2-x)*(3+p)*(1-x)
            b = brm*(x1**(3+p)-x2**(3+p))+brp*(x1**(3+p)+x2**(3+p))  
            return prefac * b

    def _val(self, x, eta, alt, j):
        #assert eta >= 0
        if x >= eta:
            return self._intval(x, eta, alt, j)
        elif -eta < x < eta:
            return self._intval(x, eta, alt, j, zero=True) 
        else:
            return 0

    def _sea(self, x, eta, alt, j):
        DMKILL = 1 # Set to 0 to agree with DM
        assert eta >= 0
        if x >= eta:
            return self._intsea(x, eta, alt, j)
        elif -eta < x < eta:
            return ( self._intsea(self, x, eta, alt, j, zero=True) 
                      - DMKILL*self._intsea(-x, eta, alt, j, zero=True) )
        else:
            return -DMKILL*self._intsea(-x, eta, alt, j)  # x < -eta

    ######  GPDs  #########

    def Huval(self, x, eta, t, Q2):
        """GK12 model for H^u_val GPD."""
        Q02 = 4.
        L = log(Q2/Q02)
        cuval = array([1.52+0.248*L, 2.88-0.940*L, -0.095*L, 0 ])
        al0, alp = (0.48, 0.9)
        alt = al0 + alp*t
        xs = array([self._val(x, eta, alt, j) for j in range(4)])
        return sum(cuval*xs)  

    def Hdval(self, x, eta, t, Q2):
        """GK12 model for H^d_val GPD."""
        Q02 = 4.
        L = log(Q2/Q02)
        cdval = array([0.76+0.248*L, 3.11-1.36*L, -3.99+1.15*L, 0])
        al0, alp = (0.48, 0.9)
        alt = al0 + alp*t
        xs = array([self._val(x, eta, alt, j) for j in range(4)])
        return sum(cdval*xs)

    def Hs(self, x, eta, t, Q2):
        """GK12 model for strange sea GPD."""
        Q02 = 4.
        Mp2 = 0.9383**2
        L = log(Q2/Q02)
        cs = array([0.123+0.0003*L, -0.327-0.004*L, 0.692-0.068*L, -0.486+0.038*L])
        al0, alp = (1.10+0.06*L-0.0027*L**2, 0.15)
        alt = al0 + alp*t
        xs = array([self._sea(x, eta, alt, j) for j in range(4)])
        b = 2.58 + 0.25*log(Mp2/(Q2+Mp2))
        return exp(b*t) * sum(cs*xs)

    def Hudsea(self, x, eta, t, Q2):
        """GK12 model for u or d sea GPDs."""
        Q02 = 4.
        kaps = 1.+0.68/(1.+0.52*log(Q2/Q02))
        return kaps * self.Hs(x, eta, t, Q2)

    def Hu(self, x, eta, t, Q2): 
        return self.Huval(x, eta, t, Q2) + self.Hudsea(x, eta, t, Q2)
        
    def Hd(self, x, eta, t, Q2): 
        return self.Hdval(x, eta, t, Q2) + self.Hudsea(x, eta, t, Q2)

    def _invAq(self, al, cs):
        coefs = [1, (1-al)/(5-al), (2-al)*(1-al)/(6-al)/(5-al)]
        return beta(1-al,4)*sum(cs*coefs)

    def Htuval(self, x, eta, t, Q2):
        """GK12 model for Ht^u_val GPD."""
        Q02 = 4.
        L = log(Q2/Q02)
        ctuval = array([0.17+0.03*L, 1.34-0.02*L, 0.12-0.4*L])
        al0, alp = (0.48, 0.45)
        alt = al0 + alp*t
        # val() expects series in beta^(j/2), so 0,2,4 -> 0,1,2
        xs = array([self._val(x, eta, alt, j) for j in [0,2,4]])
        return 0.926*sum(ctuval*xs)/self._invAq(al0, ctuval)  

    def Htdval(self, x, eta, t, Q2):
        """GK12 model for Ht^u_val GPD."""
        Q02 = 4.
        L = log(Q2/Q02)
        ctuval = array([-0.32-0.04*L, -1.427-0.176*L, 0.692-0.068*L])
        al0, alp = (0.48, 0.45)
        alt = al0 + alp*t
        # val() expects series in beta^(j/2), so 0,2,4 -> 0,1,2
        xs = array([self._val(x, eta, alt, j) for j in [0,2,4]])
        return -0.341*sum(ctuval*xs)/self._invAq(al0,ctuval)  

    def Euval(self, x, eta, t, Q2):
        """GK12 model for E^u_val GPD."""
        ceuval = array([1, -1])  # (1-rho)
        al0, alp = (0.48, 0.9)
        alt = al0 + alp*t
        # val() expects series in beta^(j/2), so 0,2 -> 0,1
        xs = array([self._val(x, eta, alt, j) for j in [0,2]])
        return 1.67*sum(ceuval*xs)/beta(1-al0,5)  

    def Edval(self, x, eta, t, Q2, DMbeta=False):
        """GK12 model for E^d_val GPD."""
        # (1-rho)^2.6 up to rho^8
        if DMbeta:
            # Choice of DM in his notebook
            betd = 6
            cedval = [1, -3, +3, -1]
            bmax = 3
        else:
            # GK choice
            betd = 5.6
            cedval = array([1, -2.6, 2.08, -0.416, -0.0416, -0.011648, -0.0046592, -0.00226304, -0.001244672])
            bmax = 8
        al0, alp = (0.48, 0.9)
        alt = al0 + alp*t
        # val() expects series in beta^(j/2)
        xs = array([self._val(x, eta, alt, j) for j in                                                   range(0, 2*bmax+1, 2)])
        return -2.03*sum(cedval*xs)/beta(1-al0,1+betd)  

    def Esea(self, x, eta, t, Q2):
        """GK12 model for strange sea GPD."""
        Q02 = 4.
        Mp2 = 0.9383**2
        L = log(Q2/Q02)
        ces = array([1, -2, 1])
        al0, alp = (1.10+0.06*L-0.0027*L**2, 0.15)
        alt = al0 + alp*t
        xs = array([self._sea(x, eta, alt, j) for j in range(0,5,2)])
        b = 0.9*(2.58 + 0.25*log(Mp2/(Q2+Mp2)))
        # Note that opposite sign is also possible:
        return -0.155*exp(b*t) * sum(ces*xs)

    def Eu(self, x, eta, t, Q2): 
        return self.Euval(x, eta, t, Q2) + self.Esea(x, eta, t, Q2)

    def Etuval(self, x, eta, t, Q2):
        """GK12 model for Et^u GPD (non-pole part)."""
        ces = array([1, -2, 1])
        al0, alp = (0.48, 0.45)
        alt = al0 + alp*t
        xs = array([self._val(x, eta, alt, j) for j in range(0,5,2)])
        b = 0.9
        return 14.*exp(b*t) * sum(ces*xs)

    def Etdval(self, x, eta, t, Q2):
        """GK12 model for Et^d GPD (non-pole part)."""
        ces = array([1, -2, 1])
        al0, alp = (0.48, 0.45)
        alt = al0 + alp*t
        xs = array([self._val(x, eta, alt, j) for j in range(0,5,2)])
        b = 0.9
        return 4.*exp(b*t) * sum(ces*xs)

    def _gegen32(self, x):
        return 15.*x**2/2 - 3./2

    def _phi(self, u, a2=0.22):
        return 6*u*(1-u)*(1+a2*self._gegen32(2*u+1))

    def Etupole(self, x, eta, t, Q2):
        if -eta < x < eta:
            Mp2 = 0.938272013**2
            gA0 = 1.267
            LAM2 = 0.44**2
            mpi2 = 0.1396**2
            numer = Mp2 * gA0 * (LAM2 - mpi2)
            denom = (mpi2-t)*(LAM2-t)
            return numer*self._phi((x+eta)/2/eta)/denom/eta
        else:
            return 0

    def Etdpole(self, x, eta, t, Q2):
        return - self.Etupole(x, eta, t, Q2)

    ######  CFFs  #########

    def ImH(self, pt, xi=0):
        """GK12 model for Im(CFFH)."""
        # FIXME: The following solution is not elegant
        if isinstance(xi, ndarray):
            # function was called with third argument that is xi nd array
            xs = xi
        elif xi != 0:
            # function was called with third argument that is xi number
            xs = [xi]
        else:
            # xi should be taken from pt object
            xs = [pt.xi]
        t = pt.t
        Q2 = pt.Q2
        res = []
        for x in xs:
            res.append( pi*(  (4./9)*(self.Hu(x, x, t, Q2) - self.Hu(-x, x, t, Q2))
                    +(1./9)*(self.Hd(x, x, t, Q2) - self.Hd(-x, x, t, Q2))
                    +(1./9)*(self.Hs(x, x, t, Q2) - self.Hs(-x, x, t, Q2)) ) )
        if len(res) == 1:
            return res[0]
        else:
            return array(res)

    #def ImHalt(self, pt):
        #"""GK12 model for Im(CFFH) (alternative expression)."""
        #x = pt.xi
        #t = pt.t
        #Q2 = pt.Q2
        #assert x>0
        #return pi*(  (4./9)*(self.Huval(x, x, t, Q2) +2*self.Hudsea(x, x, t, Q2))
                #+(1./9)*(self.Hdval(x, x, t, Q2)+2*self.Hudsea(x, x, t, Q2))
                #+(1./9)*(2*self.Hs(x, x, t, Q2)) )

    def ImHt(self, pt, xi=0):
        """GK12 model for Im(CFFHt)."""
        # FIXME: The following solution is not elegant
        if isinstance(xi, ndarray):
            # function was called with third argument that is xi nd array
            xs = xi
        elif xi != 0:
            # function was called with third argument that is xi number
            xs = [xi]
        else:
            # xi should be taken from pt object
            xs = [pt.xi]
        t = pt.t
        Q2 = pt.Q2
        res = []
        for x in xs:
            assert x>0
            res.append( pi*(  (4./9)*self.Htuval(x, x, t, Q2)
                +(1./9)*self.Htdval(x, x, t, Q2) ) )
        if len(res) == 1:
            return res[0]
        else:
            return array(res)

    def ImE(self, pt, xi=0, DMbeta=False):
        """GK12 model for Im(CFFE)."""
        # FIXME: The following solution is not elegant
        if isinstance(xi, ndarray):
            # function was called with third argument that is xi nd array
            xs = xi
        elif xi != 0:
            # function was called with third argument that is xi number
            xs = [xi]
        else:
            # xi should be taken from pt object
            xs = [pt.xi]
        t = pt.t
        Q2 = pt.Q2
        res = []
        for x in xs:
            assert x>0
            res.append( pi*(  (4./9)*(self.Euval(x, x, t, Q2) 
                                        +2*self.Esea(x, x, t, Q2))
                +(1./9)*(self.Edval(x, x, t, Q2, DMbeta)+2*self.Esea(x, x, t, Q2))
                +(1./9)*(2*self.Esea(x, x, t, Q2)) )  )
        if len(res) == 1:
            return res[0]
        else:
            return array(res)

    def ImEt(self, pt, xi=0):
        """GK12 model for Im(CFFEt)."""
        # FIXME: The following solution is not elegant
        if isinstance(xi, ndarray):
            # function was called with third argument that is xi nd array
            xs = xi
        elif xi != 0:
            # function was called with third argument that is xi number
            xs = [xi]
        else:
            # xi should be taken from pt object
            xs = [pt.xi]
        t = pt.t
        Q2 = pt.Q2
        res = []
        for x in xs:
            assert x>0
            res.append( pi*(  (4./9)*(self.Etuval(x, x, t, Q2) )
                +(1./9)*(self.Etdval(x, x, t, Q2)) ) ) 
        if len(res) == 1:
            return res[0]
        else:
            return array(res)

    def ReEtpole(self, pt):
        Mp2 = 0.938272013**2
        gA0 = 1.267
        LAM2 = 0.44**2
        mpi2 = 0.1396**2
        numer = Mp2 * gA0 * (LAM2 - mpi2)
        denom = (mpi2-pt.t)*(LAM2-pt.t)
        green = numer/denom
        a2 = 0.22
        return 2*(1+a2)*green/pt.xi

    def ReEtnonpole(self, pt):
        return ComptonDispersionRelations.ReEt(self, pt)

    def ReEt(self, pt):
        return self.ReEtpole(pt) + self.ReEtnonpole(pt)


    def _GKaux(self, x, ds, Q2=4.):
        """auxilliary integrand that accepts ndarray"""
        res = []
        for d in ds:
            res.append(self.Hu(x, 0, -d**2, Q2))
        return res

    def _GKauxE(self, x, ds, Q2=4.):
        """auxilliary integrand that accepts ndarray"""
        res = []
        for d in ds:
            res.append(self.Eu(x, 0, -d**2, Q2))
        return res

    def gpdb(self, x, b, Q2=4., inf=2.5):
        """GPD in b space ??? FIXME."""
        return bquadrature(lambda d: d*j0(b*d/GeVfm)*self._GKaux(x, d, Q2), 
                0.0, inf) / (2*pi*GeVfm**2)

    def gpdbpol(self, pol, x, bx, by, Q2=4., inf=2.5):
        """polarized GPD in b space"""
        b = sqrt(bx**2 + by**2)
        aux = bquadrature(lambda d: d**2*j1(b*d/GeVfm)*self._GKauxE(x, d, Q2), 
                0.0, inf) / (2*pi*GeVfm**2)
        return self.gpdb(x, b, Q2, inf) + pol*bx*aux/(2*b*Mp)

class GK12D(GK12):
    """ Goloskokov-Kroll PDF model with added D-term

        From [arXiv:1210.6975], Kroll:2012sm
        Real part of CFFs is obtained using dispersion relations.

    """

    def subtraction(self, pt):
        """D-term"""
        denom = (1 - pt.t/(0.841*0.487**2))**0.841
        return -(10./9.)*(-1.9/denom)

class GK0(GK12):
    """ Goloskokov-Kroll PDF model with only CFF H """

    def ImE(self, pt, xi=0): 
        return 0.

    def ReE(self, pt, xi=0): 
        return 0.

    def ImHt(self, pt, xi=0): 
        return 0.

    def ReHt(self, pt, xi=0): 
        return 0.

    def ImEt(self, pt, xi=0): 
        return 0.

    def ReEt(self, pt, xi=0): 
        return 0.


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
        endpointpower: Im(CFF)s are defined as NN*(1-xB)**endpointpower to enforce
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
        elif name in ComptonFormFactors.allCFFs:
            # if asked for CFF which is not in output_layer, return 0
            self.curname = name
            return self.zero
        elif name in ComptonFormFactors.allCFFeffs:
            # if asked for CFF which is not in output_layer, return 0
            self.curname = name
            return self.zero
        elif name in ['endpointpower', 'optimization', 'useDR']:
            if name in self.__dict__:
                return self.__dict__[name]
            else:
                return None
        else:
            #syst.stderr.write('Possibly caught exception: AttErr: $s\n' % name)
            raise AttributeError(name)

    #def zero(self, *args, **kwargs):
        #return 0

    def subtraction(self, pt):
        """This subtraction constant is negative of the one which is standard
            in the literature!
        """
        #return 2.25  # for testing, return fixed value
        if self.parameters['outputvalueC'] is not None:
            # this occurs during training: value is set by training
            # routine by calling with outputvalueC set by training routine
            return self.parameters['outputvalueC']
        ar = []
        for netC in self.netsC:
            ar.append(netC.activate([pt.t])[0])
        all = array(ar).flatten()
        if 'nnet' in self.parameters:
            if self.parameters['nnet'] == 'ALL':
                return all
            elif self.parameters['nnet'] == 'AVG':
                return all.mean()
            else: # we want particular netC
                try:
                    return all[self.parameters['nnet']]
                except IndexError:
                    raise IndexError(str(self)+' has only '+str(len(self.netsC))+' nets!')
        # by default, we get mean value (FIXME:this should never occurr?)
        else:
            _lg.debug('FIXME: This line should never be reached')
            return all.mean()

    def _pipole(self, pt):
        """FIXME: code duplication from GK12 model"""
        Mp2 = 0.938272013**2
        gA0 = 1.267
        LAM2 = 0.44**2
        mpi2 = 0.1396**2
        numer = Mp2 * gA0 * (LAM2 - mpi2)
        denom = (mpi2-pt.t)*(LAM2-pt.t)
        green = numer/denom
        a2 = 0.22
        return 2*(1+a2)*green/pt.xi

    def CFF(self, pt, xi=0):
        # FIXME: This function is HEAVILY sub-optimal and non-pythonic and messy!
        _lg.debug('NN model CFF called as = %s\n' % self.curname)
        # Uncomment the following block to switch on the pion pole
        #if self.curname == 'ReEt':
        #    return self._pipole(pt)
        if hasattr(self, 'useDR') and self.useDR and self.curname in self.useDR:
            _lg.debug('Doing DR for CFF: %s\n' % self.curname)
            if self.curname == 'ReH':
                #set_trace()
                a = DR.intVNN(self.ImH, pt)
                b = self.subtraction(pt)
                if isinstance(b, ndarray):
                    b = b.reshape(b.size,1)
                return a - b
                #return DR.intVNN(self.ImH, pt) - self.subtraction(pt)
            elif self.curname == 'ReE':
                a = DR.intVNN(self.ImE, pt)
                b = self.subtraction(pt)
                if isinstance(b, ndarray):
                    b = b.reshape(b.size,1)
                return a + b
                #return DR.intVNN(self.ImE, pt) + self.subtraction(pt)
            elif self.curname in ['ReHt', 'ReEt']:
                return DR.intANN(self.__getattr__('Im'+self.curname[2:]), pt)
            else:
                raise ValueError('Only Re(CFF) can be calculated via disp. rel.')
        ind = self.output_layer.index(self.curname)
        if isinstance(xi, ndarray) and len(self.nets)>0:
            # function was called with third argument that is xi nd array
            # and we already have some nets so disp.int. can be calculated
            x = xi
        elif self.parameters['outputvalue'] is not None:
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
                if self.endpointpower and self.curname[:2] == 'Im':
                    ar.append(net.activate([xB, pt.t])[ind]*(1-xB)**self.endpointpower)
                else:
                    ar.append(net.activate([xB, pt.t])[ind])
            all = array(ar).flatten()
            if 'nnet' in self.parameters:
                if self.parameters['nnet'] == 'ALL':
                    res.append(all)
                elif self.parameters['nnet'] == 'AVG':
                    res.append(all.mean())
                else: # we want particular net
                    try:
                        res.append(all[self.parameters['nnet']])
                    except IndexError:
                        raise IndexError(str(self)+' has only '+str(len(self.nets))+' nets!')
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
            if 'nnet' in self.parameters:
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
    gepardPool = [g1, g2, g3, g4, g5, g6, g7]  #  modules to choose from
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
             #'EKAPS' : 0.0,
            'ESKEWS' : 0.0,
             'EAL0G' : 1.1,
             'EALPG' : 0.15,
             'EM02G' : 0.7,    'limit_EM02G' : (0.1, 1.5),
           'EDELM2G' : 0.0,
               'EPG' : 2.0,
             'ESECG' : 0.0,
             'ETHIG' : 0.0,
             #'EKAPG' : 0.0,
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
           'EDELM2S', 'EPS', 'ESECS', 'ETHIS',# 'EKAPS', 
           'ESKEWS',
           'EAL0G', 'EALPG', 'EM02G',
           'EDELM2G', 'EPG', 'ESECG', 'ETHIG',#'EKAPG', 
           'ESKEWG']

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
            raise ValueError("Invalid ansatz: %s\n" % ansatz)
        
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
        if 'newgepard' in self.kwargs:
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
        """Skewness of GPD H_Q"""
        return self.gpdHtrajQ(pt)/self.gpdHzeroQ(pt)

    def gpdHskewG(self, pt):
        """Skewness of GPD H_G"""
        return self.gpdHtrajG(pt)/self.gpdHzeroG(pt)

    def gpdHQb(self, pt, inf=2.5):
        """GPD H^sea in b space in fm^-2."""
        b = sqrt(pt.bx**2 + pt.by**2)
        return bquadrature(lambda d: d*j0(b*d/GeVfm)*self.gpdHzeroQ(pt, ts=-d**2), 
                0.0, inf) / (2*pi*GeVfm**2)

    def gpdHQbpol(self, pt, inf=2.5):
        """polarized GPD H^sea in b space in fm^-2."""
        b = sqrt(pt.bx**2 + pt.by**2)
        aux = bquadrature(lambda d: d**2*j1(b*d/GeVfm)*self.gpdEzeroQ(pt, ts=-d**2), 
                0.0, inf) / (2*pi*GeVfm**2)
        return self.gpdHQb(pt) + pt.bx*aux/(2*b*Mp)

    def gpdHQbpolpol(self, pt, r, tht, inf=2.5):
        """Same as gpdHQbpol but in polar coordinate system.
        FIXME: some ugly coding here.
        """
        if not isinstance(r, ndarray):
            assert isinstance(r, (int, float)) 
            ra = array([r])
        else:
            ra = r
        res = []
        for rr in ra:
            pt.bx = rr*cos(tht)
            pt.by = rr*sin(tht)
            b = sqrt(pt.bx**2 + pt.by**2)
            aux = bquadrature(lambda d: d**2*j1(b*d/GeVfm)*self.gpdEzeroQ(pt, ts=-d**2), 
                    0.0, inf) / (2*pi*GeVfm**2)
            res.append(self.gpdHQb(pt) + pt.bx*aux/(2*b*Mp))
        if not isinstance(r, ndarray):
            r = ra[0]
        else:
            r = ra
        return array(res)

    def rth(self, pt, tht, inf=2.5):
        """Calculate <r(tht)> """
        norm = rthtquadrature(lambda r: r*self.gpdHQbpolpol(pt, r, tht), 0.0, inf)
        prob = rthtquadrature(lambda r: r*r*self.gpdHQbpolpol(pt, r, tht), 0.0, 2.5)
        return prob/norm

    def gpdHQbIntg(self, pt, ba, inf=2.5):
        """Same as gpdHQb but convenient as integrand."""
        res = []
        for b in ba:
            aux = bquadrature(lambda d: d*j0(b*d/GeVfm)*self.gpdHzeroQ(pt, ts=-d**2), 
                0.0, inf) / (2*pi*GeVfm**2)
            res.append(aux)
        return array(res)

    def bsq(self, pt, inf=2.5):
        """Calculate unpolarized <b^2> for model m."""
        norm = rthtquadrature(lambda b: b*self.gpdHQbIntg(pt, b), 0.0, inf)
        aux = rthtquadrature(lambda b: b**3*self.gpdHQbIntg(pt, b), 0.0, 2.5)
        return aux/norm
    
    def gpdHGb(self, pt, inf=2.5):
        """GPD H^G in b space in fm^-2."""
        b = sqrt(pt.bx**2 + pt.by**2)
        return bquadrature(lambda d: d*j0(b*d/GeVfm)*self.gpdHzeroG(pt, ts=-d**2), 
                0.0, inf) / (2*pi*GeVfm**2)

    def gpdHGbpol(self, pt, inf=2.5):
        """polarized GPD H^sea in b space in fm^-2."""
        b = sqrt(pt.bx**2 + pt.by**2)
        aux = bquadrature(lambda d: d**2*j1(b*d/GeVfm)*self.gpdEzeroG(pt, ts=-d**2), 
                0.0, inf) / (2*pi*GeVfm**2)
        return self.gpdHGb(pt) + pt.by*aux/(2*b*Mp)


    def beff(self, pt):
        """slope of  CFFH"""
        t1 = -0.03
        t2 = -1.5
        tmem = pt.t
        pt.t = t2
        self.g.newcall = 1
        up = self.ImH(pt)
        pt.t = t1
        self.g.newcall = 1
        down = self.ImH(pt)
        pt.t = tmem
        return log(up/down)/(t2-t1)


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
        self.parameters = hubDictNew(self.DR.parameters, self.Gepard.parameters)
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
                 (1e-5 < pt.xB < 0.65)
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

    # TODO: maybe neutron should be done using some particle attribute of CFFs

    def nImH(self, pt, xi=0):
        return  self.Gepard.ImH(pt) + self.DR.nImH(pt, xi)

    def nReH(self, pt):
        return  self.Gepard.ReH(pt) + self.DR.nReH(pt)

    def nImE(self, pt, xi=0):
        return  self.Gepard.ImE(pt) + self.DR.nImE(pt, xi)

    def nReE(self, pt):
        return  self.Gepard.ReE(pt) + self.DR.nReE(pt)

    def nImHt(self, pt, xi=0):
        return  self.DR.nImHt(pt, xi)

    def nReHt(self, pt):
        return  self.DR.nReHt(pt)

    def nImEt(self, pt):
        return  self.DR.nImEt(pt)

    def nReEt(self, pt):
        return  self.DR.nReEt(pt)


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


class BMP(ComptonFormFactors):
    """This is poor-man's higher twist model, according to Braun et al.

    In general, instance of any model has to return CFFs as defined by BMJ,
    so here rotation from BMP is defined

    """

    allCFFs = ['ImH', 'ReH', 'ImE', 'ReE', 'ImHt', 'ReHt', 'ImEt', 'ReEt']
    allCFFeffs = ['ImHeff', 'ReHeff', 'ImEeff', 'ReEeff',
                     'ImHteff', 'ReHteff', 'ImEteff', 'ReEteff']

    ###  BMJ CFFs:

    #  Fp = F_{+,+} = F_{-,-}
    #  F0 = F_{0,+} = F_{0,-}
    #  Fm = F_{-,+} = F_{+,-}  -- FIXME: not implemented yet


    ###  BMP CFFs:

    #  BFp, BF0, BFm  in complete analogy to above. They are
    #                 parameters of this model!

    #  Fp = BHp + (chi/2)*(BFp+BFm) - chi0*BF0
    #  F0 = -(1+chi)*BH0 + chi0*(BFp+BHm)

    # Prefixing parameters with 'p' below (like pIMBH0) is
    # not really necessary, but I do it for consistent naming
    # of "local" fitting models

    def __init__(self, **kwargs):
        self.parameters = {
            'pImHBp' :  0.0,   'pImHB0' :  0.0,
            'pReHBp' :  0.0,   'pReHB0' :  0.0,
            'pImEBp' :  0.0,   'pImEB0' :  0.0,
            'pReEBp' :  0.0,   'pReEB0' :  0.0,
           'pImHtBp' :  0.0,  'pImHtB0' :  0.0,
           'pReHtBp' :  0.0,  'pReHtB0' :  0.0,
           'pImEtBp' :  0.0,  'pImEtB0' :  0.0,
           'pReEtBp' :  0.0,  'pReEtB0' :  0.0}

        self.parameter_names = [
             'pImHBp', 'pReHBp', 'pImHtBp', 'pReHtBp',
             'pImEBp', 'pReEBp', 'pImEtBp', 'pReEtBp',
             'pImHB0', 'pReHB0', 'pImHtB0', 'pReHtB0',
             'pImEB0', 'pReEB0', 'pImEtB0', 'pReEtB0']

        # now do whatever else is necessary
        ComptonFormFactors.__init__(self, **kwargs)


    # BMJ leading twist CFFs:
    #     F = Fp i.e. F_{+,+}

    for name in allCFFs:
        exec('def {0}(self, pt): return self.BMJ_CFFp("{0}", pt)'.format(name))

    # BMJ eff. tw-3 CFFs:
    #    Feff = sqrt(1+eps2)*Q*(2-xB+xB*t/Q2)/sqrt(2)/Ktilde * F0

    for name in allCFFs:
        exec('def {0}eff(self, pt): return sqrt(1.+pt.eps2)*sqrt(pt.Q2)*(2.-pt.xB+pt.xB*pt.t/pt.Q2)/sqrt(2.)/pt.tK * self.BMJ_CFF0("{0}", pt)'.format(name))


    def BMJ_CFFp(self, CFF, pt):
        """ BMJ Fp i.e. F_{+,+} in terms of BMP CFFs """

        BCFFp = self.parameters['p{}Bp'.format(CFF)]
        BCFF0 = self.parameters['p{}B0'.format(CFF)]
        BCFFm = 0   # NOT IMPLEMENTED
        return BCFFp + (pt.chi)/2.*(BCFFp+BCFFm) - pt.chi0*BCFF0

    def BMJ_CFF0(self, CFF, pt):
        """ BMJ F0 in terms of BMP CFFs """

        BCFFp = self.parameters['p{}Bp'.format(CFF)]
        BCFF0 = self.parameters['p{}B0'.format(CFF)]
        BCFFm = 0   # NOT IMPLEMENTED
        return -(1.+pt.chi)*BCFF0 + pt.chi0*(BCFFp+BCFFm)

    def BMJ_CFFm(self, CFF, pt):
        """ BMJ Fm i.e. F_{-,+} in terms of BMP CFFs """

        # NOT IMPLEMENTED
        return 0.


##  --- Complete models built from the above components ---

class PureBetheHeitler(ComptonFormFactors, ElasticKelly):
    """Pure Bethe-Heitler (all CFFs=0) model."""

    def is_within_model_kinematics(self, pt):
        return ( (1.5 <= pt.Q2) and
                 (pt.tm < min(1., pt.Q2/4)) and
                 (1e-5 < pt.xB < 0.65)
               )

class PureBetheHeitlerNeutron(ComptonFormFactors, ElasticKelly):
    """Pure Bethe-Heitler (all CFFs=0) model for neutrons."""

    F1 = ElasticKelly.nF1
    F2 = ElasticKelly.nF2

class ModelDR(ComptonModelDR, ElasticDipole):
    """Complete model as in arXiv:0904.0458.."""

class Gepard(ComptonGepard, ElasticKelly):
    """Complete model as in arXiv:0904.0458.."""

class GK(GK12, ElasticKelly):
    """Goloskokov-Kroll model as described in arXiv:1210.6975."""

class GKD(GK12D, ElasticKelly):
    """Goloskokov-Kroll model as described in arXiv:1210.6975 + D-term."""

class GKonlyH(GK0, ElasticKelly):
    """Goloskokov-Kroll model with only CFF H != 0"""

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


class HybridKellyNeutron(ComptonHybrid, ElasticKelly):
    """As HybridKelly, but with n --> p isospin flip
    in order to be able to use unchanged proton DVCS observables code for
    calculation of neutron DVCS observables."""

    F1 = ElasticKelly.nF1
    F2 = ElasticKelly.nF2

    ImH = ComptonHybrid.nImH
    ImHt = ComptonHybrid.nImHt
    ImE = ComptonHybrid.nImE
    ImEt = ComptonHybrid.nImEt

    ReH = ComptonHybrid.nReH
    ReHt = ComptonHybrid.nReHt
    ReE = ComptonHybrid.nReE
    ReEt = ComptonHybrid.nReEt


class HybridZero(ComptonHybrid, ElasticZero):
    """As Hybrid, but with zero elasticd FFs, so giving only DVCS^2 contrib."""


class ModelLocal(ComptonLocal, ElasticKelly):
    """Model for local fitting of CFFs themselves."""

    def ImEt(self, pt):
        #To get sensible numbers for parameter.
        return (2.-pt.xB)/pt.xB * self.parameters['pImEt']

    def ReEt(self, pt):
        #To get sensible numbers for parameter.
        return (2.-pt.xB)/pt.xB * self.parameters['pReEt']

class ModelBMP(BMP, ElasticKelly):
    """Model for local fitting of CFFs in Braun et al. framework"""
