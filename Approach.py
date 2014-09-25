#from IPython.Debugger import Tracer; debug_here = Tracer()

import copy

from numpy import sin, cos, pi, sqrt, array, linspace, ndarray
from scipy.special import gammainc
from scipy.stats import scoreatpercentile

import utils, quadrature, Data
from constants import *

class Approach(object):
    """Class of approaches to calculation of observables.

    This is a base class and subclasses should implement
    various observables as methods. Instances of subclasses
    are instantiated with specific functions for modelling
    CFFs and FFs, so they contain complete algorithm for
    calculation of specific observable.

    Boolean variable 'optimization' (default is False) controls
    whether faster (but analytically more complicated) formulas
    are used. This is for development purposes only. When
    fitting optimization=True is better choice.

    Implemented subclases:  BMK, hotfixedBMK
    TODO: VGG, Guichon?
    """


    def __init__(self, m, optimization=False):
        """Initialize with m as instance of Model class."""
        self.model = m
        self.m = self.model  # shortcut
        self.optimization = optimization
        self.name = 'N/A'  # to be used as identifier in theory database
        self.description = 'N/A'  # something human-understandable

    def __repr__(self):
        return "<Theory: %s + %s at %s>" % utils.flatten(
                (str(self.__class__).split()[1][1:-2],
                  tuple(str(self.model).split()[::3])))

    def copy(self):
        """Return deep copy of itself."""
        return copy.deepcopy(self)

    def save(self, db):
        """Save theory to database."""
        db[self.name] = self

    def chisq(self, points, sigmas=False, **kwargs):
        """Return tuple (chi-square, d.o.f., probability). If the approach and model
           provide uncertainties, they are ignored - only experimental uncertainties
           are taken into account."""
        nfreepars=utils.npars(self.model)
        dof = len(points) - nfreepars
        allsigmas = [(self.predict(pt, observable=pt.yaxis, **kwargs) - pt.val) / pt.err for
                    pt in points]
        chi = sum(s*s for s in allsigmas)  # equal to m.fval if minuit fit is done
        fitprob = (1.-gammainc(dof/2., chi/2.)) # probability of this chi-sq
        if sigmas:
            return allsigmas
        else:
            return (chi, dof, fitprob)

    def scan(self, parname, points, npoints=5):
        """Scan chi-square dependence on single parameter."""
        mem = self.m.parameters[parname]
        for val in linspace(0.5*mem, 1.5 * mem, npoints):
            self.m.parameters[parname] = val
            self.m.ndparameters[self.m.parameter_names.index(parname)] = val
            chi, dof, fitprob = self.chisq(points)
            print '%s  ->  %s' % (val, chi)
        self.m.parameters[parname] = mem  # restore original value
        self.m.ndparameters[self.m.parameter_names.index(parname)] = mem


    def print_chisq(self, points, sigmas=False, **kwargs):
        """Pretty-print the chi-square."""
        if sigmas:
            print self.chisq(points, sigmas=True, **kwargs)
        print 'P(chi-square, d.o.f) = P(%1.2f, %2d) = %5.4f' % self.chisq(points, **kwargs)

        
    def predict(self, pt, error=False, CL=False, **kwargs):
        """Give prediction for DataPoint pt.

        Keyword arguments:
        error - if available, produce tuple (mean, error)
        CL - error is not std.dev., but 68% C.L. (mean, errplus, errminus)
        observable - string. Default is pt.yaxis. It is acceptable also
                     to pass CFF as observable, e.g., observable = 'ImH'
        parameters - dictionary which will temporarily update model's one
        orig_conventions - give prediction using original conventions of
                           the given DataPoint (e.g. for plotting)

        """
        m = self.model
        if kwargs.has_key('observable'):
            obs = kwargs['observable']
        else:
            obs = pt.yaxis

        if kwargs.has_key('parameters'):
            old = m.parameters.copy()
            m.parameters.update(kwargs['parameters'])
        #elif isinstance(m, Model.ComptonNeuralNets):
        #    # It is not training (which always uses 'parameters'), and
        #    # we are not asked for particular net (call would again come
        #    # with 'parameters'), so we want mean of all nets
        #    m.parameters['nnet'] = 'ALL'
        #    result = getattr(self, obs)(pt)
        #    if error:
        #        return (result.mean(), result.std())
        #    else:
        #        return result.mean()

        if obs in m.allCFFs or obs in m.allGPDs:
            # we need a model attribute
            fun = getattr(m, obs)
        else:
            # we need a "real" observable
            fun = getattr(self, obs)

        if error:
            try:
                pars = [p for p in m.parameter_names if m.parameters['fix_'+p] == False]
                var = 0
                dfdp = {}
                for p in pars:
                    # calculating dfdp = derivative of observable w.r.t. parameter:
                    h=sqrt(m.covariance[p,p])
                    mem = m.parameters[p]
                    m.parameters[p] = mem+h/2.
                    if self.model.__dict__.has_key('g'): self.m.g.newcall = 1
                    up = fun(pt)
                    m.parameters[p] = mem-h/2.
                    if self.model.__dict__.has_key('g'): self.m.g.newcall = 1
                    down = fun(pt)
                    m.parameters[p] = mem
                    dfdp[p] = (up-down)/h
                for p1 in pars:
                    for p2 in pars:
                        var += dfdp[p1]*m.covariance[p1,p2]*dfdp[p2]
                var = tolerance2 * var
                result = (fun(pt), sqrt(var))
            except KeyError:
                # we have neural net
                allnets = fun(pt)
                if CL:
                    # 68% confidence level
                    m = allnets.mean()
                    try:
                        result = (m, 
                              scoreatpercentile(allnets, 84)[0] - m,
                              m - scoreatpercentile(allnets, 16)[0])
                    except IndexError:
                        result = (m, 
                              scoreatpercentile(allnets, 84) - m,
                              m - scoreatpercentile(allnets, 16))
                else:
                    # one sigma
                    result = (allnets.mean(), allnets.std())
        else:
            if self.model.__dict__.has_key('g'): self.m.g.newcall = 1
            result = fun(pt)
            if isinstance(result, ndarray):
                # we have neural net
                result = result.mean()

        if kwargs.has_key('parameters'):
            # restore old values
            self.model.parameters.update(old)

        if kwargs.pop('orig_conventions', False):
            # express result in conventions of original datapoint
            try:
                result = (self.orig_conventions(pt, result[0]),) + result[1:]
            except IndexError:
                result = self.orig_conventions(pt, result)
        return result


class BMK(Approach):
    """Implementation of formulas from hep-ph/0112108  (BMK)"""

    ### Kinematics ###
    # (implemented as static and class methods)

    def tmin(Q2, xB, eps2):
        """BMK Eq. (31)"""
        return -Q2 * ( 2. * (1.-xB)*(1. - sqrt(1.+eps2)) + eps2 ) / (
                4. * xB * (1.-xB) + eps2 )
    tmin = staticmethod(tmin)

    def xBmin(s, Q2):
        """Constrained by xB=Q2/(s-Mp^2)/yMax, with 1-yMax+yMax^2*eps2/4=0."""
        yMax = 1 + Mp2*Q2/(s-Mp2)**2
        return Q2 / (s-Mp2) / yMax
    xBmin = staticmethod(xBmin)

    def K2(Q2, xB, t, y, eps2):
        """BMK Eq. (30)"""
        tm = BMK.tmin(Q2, xB, eps2)
        brace = sqrt(1.+eps2) + (4. * xB * (1.-xB) + eps2 ) / (
                4. * (1.-xB) ) * (t - tm) / Q2
        return -(t/Q2) * (1.-xB) * (1.-y-y*y*eps2/4.) * (
                1. - tm / t ) * brace
    K2 = staticmethod(K2)

    def J(Q2, xB, t, y, eps2):
        """BMK below Eq. (32)"""
        return (1.-y-y*eps2/2.) * (1. + t/Q2) - (1.-xB)*(2.-y)*t/Q2
    J = staticmethod(J)

    @staticmethod
    def is_within_phase_space(pt):
        """Is pt kinematics within allowed phase space?"""
        return (pt.xB > BMK.xBmin(pt.s, pt.Q2) and pt.t < BMK.tmin(pt.Q2, pt.xB, pt.eps2))

    def r(Q2, xB, t, y, eps2):
        """DM's fitting notes, below Eq. (13)"""
        K = sqrt(BMK.K2(Q2, xB, t, y, eps2))
        brace = (2.-y)**2 * K / (1.-y) + (1./K)*(t/Q2)*(1.-y)*(2.-xB)
        return - (2.-y) / (2.-2.*y+y**2) * brace
    r = staticmethod(r)

    def P1P2(pt):
        """ Product of Bethe-Heitler propagators, Eq(32) """
        P1 = - ( BMK.J(pt.Q2, pt.xB, pt.t, pt.y, pt.eps2) + 2. * 
                sqrt(BMK.K2(pt.Q2, pt.xB, pt.t, pt.y, pt.eps2)) * cos(pt.phi) ) / (
                        pt.y * (1. + pt.eps2) )
        P2 = 1. + pt.t / pt.Q2  - P1
        return P1 * P2
    P1P2 = staticmethod(P1P2)

    def anintP1P2(pt):
        """ Analitical integral of \int P1 P2 d\phi """
        xB, Q2, t, y, eps2, K2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2, pt.K2
        brace = ( (1 - y - (1+eps2/2.) * y**2 * eps2/2.) * (1. + t/Q2)**2 +
                  2.*K2 - (1.-xB)*(2.-y)**2 * (1. + xB*t/Q2) * t/Q2 )
        return -2. * pi * brace / (1+eps2)**2 / y**2
    anintP1P2 = staticmethod(anintP1P2)

    def PreFacSigma(self, pt):
        """ Prefactor of 4-fold xs. Take prefactor in Eq(22) times 2\pi because of proton Phi integration
        and times y/Q2 because of dy -> dQ2. Convert to nanobarns."""
        return alpha**3 * pt.xB * pt.y**2 / (8. * pi * pt.Q2**2 * sqrt(1.+pt.eps2)) * GeV2nb

    def PreFacBH(self, pt):
        """ Prefactor from Eq. (25), without e^6 """
        return 1./(pt.xB**2 * pt.y**2 * (1.+pt.eps2)**2 * pt.t * pt.P1P2)

    def PreFacDVCS(self, pt):
        """ Prefactor from Eq. (26), without e^6 """
        return 1./(pt.y**2 * pt.Q2 )

    def PreFacINT(self, pt):
        """ Prefactor from Eq. (27), without e^6 """
        return 1./(pt.xB * pt.y**3 * pt.t * pt.P1P2)

    def w(self, pt):
        """ Weight factor removing BH propagators from INT and BH amplitudes. 
        It is normalized to \int_0^2pi w  2pi as in BMK. """
        return 2.*pi*pt.P1P2 / pt.intP1P2

    def to_conventions(pt):
        """Transform stuff into BMK conventions."""
        # C1. azimutal angle phi should be in radians.
        if pt.has_key('phi') and pt.units['phi'][:3]=='deg':
            pt.phi = pt.phi * pi / 180.
            pt.newunits['phi'] = 'rad'
        # C2. phi_{Trento} -> (pi - phi_{BKM})
        if pt.has_key('frame') and pt.frame == 'Trento':
            if pt.has_key('phi'):
                pt.phi = pi - pt.phi
            elif pt.has_key('FTn'):
                if pt.FTn == 1 or pt.FTn == 3 or pt.FTn == -2:
                    pt.val = - pt.val
        # C3. varphi_{Trento} -> (varphi_{BKM} + pi)
            if pt.has_key('varphi'):
                pt.varphi = pt.varphi - pi
            elif pt.has_key('varFTn'):
                if pt.varFTn == 1 or pt.varFTn == -1:
                    pt.val = - pt.val
                else:
                    raise ValueError('varFTn = %d not allowed. Only +/-1!' % pt.varFTn)
            pt.newframe = 'BMK'
    to_conventions = staticmethod(to_conventions)

    def from_conventions(pt):
        """Transform stuff from Approach's conventions into original data's."""
        # C2. phi_{BKM} -> (pi - phi_{Trento})
        if pt.has_key('frame') and pt.frame == 'Trento':
            if pt.has_key('phi'):
                pt.phi = pi - pt.phi
            elif pt.has_key('FTn'):
                if pt.FTn == 1 or pt.FTn == 3:
                    pt.val = - pt.val
        # C3. varphi_{Trento} -> (varphi_{BKM} + pi)
            if pt.has_key('varphi'):
                pt.varphi = pt.varphi + pi
            elif pt.has_key('varFTn'):
                if pt.varFTn == 1 or pt.varFTn == -1:
                    pt.val = - pt.val
            pt.newframe = 'Trento'
        # C1. azimutal angle phi back to degrees
        if pt.has_key('phi') and pt.units['phi'][:3]=='deg':
            pt.phi = pt.phi / pi * 180.
        return pt
    from_conventions = staticmethod(from_conventions)

    def orig_conventions(pt, val):
        """Like from_conventions, but for the prediction val."""
        # This doesn't touches pt
        # C2. phi_{BKM} -> (pi - phi_{Trento})
        if pt.has_key('frame') and pt.frame == 'Trento' and pt.has_key('FTn'):
            if pt.FTn == 1 or pt.FTn == 3 or pt.FTn == -2:
                val = - val
        if pt.has_key('frame') and pt.frame == 'Trento' and pt.has_key('varFTn'):
            if pt.varFTn == 1 or pt.varFTn == -1:
                val = - val
        return val
    orig_conventions = staticmethod(orig_conventions)

    def prepare(pt):
        """Pre-calculate GPD-independent kinamatical constants and functions."""
        pt.y = (pt.W**2 + pt.Q2 - Mp2) / (pt.s - Mp2)
        pt.eps = 2. * pt.xB * Mp / sqrt(pt.Q2)
        pt.eps2 = pt.eps**2
        if pt.has_key('t'):
            pt.J = BMK.J(pt.Q2, pt.xB, pt.t, pt.y, pt.eps2)
            pt.K2 = BMK.K2(pt.Q2, pt.xB, pt.t, pt.y, pt.eps2)
            pt.K = sqrt(pt.K2)
            pt.tK2 = pt.K2*pt.Q2/(1-pt.y-pt.eps2*pt.y**2/4.)
            pt.tK = sqrt(pt.tK2)
            pt.r = BMK.r(pt.Q2, pt.xB, pt.t, pt.y, pt.eps2)
            # First option is numerical, second is analytical and faster
            #pt.intP1P2 = quadrature.Hquadrature(lambda phi: P1P2(pt, phi), 0, 2.0*pi)
            pt.intP1P2 = BMK.anintP1P2(pt)
        if pt.has_key('phi'):
            pt.P1P2 = BMK.P1P2(pt)
    prepare = staticmethod(prepare)


    ################################################
    #                                              #
    ###  Terms of squared ep->epgamma amplitude  ###
    #                                              #
    ################################################


    #### Bethe-Heitler amplitude squared Fourier coefficients

    ######  Unpolarized target 

    def cBH0unpSX(self, pt):
        """ BKM Eq. (35) - small-x approximation """
        return 16. * pt.K2 * (pt.Q2/pt.t) * ( 
                self.m.F1(pt.t)**2 - (pt.t/(4.0*Mp2)) * self.m.F2(pt.t)**2 
                  ) + 8. * (2. - pt.y)**2 * ( 
                self.m.F1(pt.t)**2 - (pt.t/(4.0*Mp2)) * self.m.F2(pt.t)**2 )

    def cBH0unp(self, pt):
        """ BKM Eq. (35) """
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        FE2 = self.m.F1(t)**2 - t * self.m.F2(t)**2 / (4.0 * Mp2) 
        FM2 = (self.m.F1(t) + self.m.F2(t))**2
        # braces are expressions in {..} in Eq. (35)
        brace1 = (2.+3.*eps2) * (Q2/t) * FE2 + 2.* xB**2 * FM2
        brace2 = (  (2.+eps2)*( (4.*xB**2*Mp2/t)*(1.+t/Q2)**2 +
                       4.*(1.-xB)*(1.+xB*t/Q2) ) * FE2 +
                    4.*xB**2*( xB + (1.-xB+eps2/2.)*(1.-t/Q2)**2 -
                       xB*(1.-2.*xB)*t**2/Q2**2 )  * FM2  )
        brace3 = 2.*eps2*(1.-t/(4.*Mp2)) * FE2 - xB**2*(1.-t/Q2)**2 * FM2
        return ( 8. * pt.K2 * brace1 + 
                (2.-y)**2 * brace2 + 
                 8. * (1.+eps2) * (1.-y-eps2*y**2/4.) * brace3  )

    def cBH1unp(self, pt):
        """ BKM Eq. (36) """
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        FE2 = self.m.F1(t)**2 - t * self.m.F2(t)**2 / (4.0 * Mp2) 
        FM2 = (self.m.F1(t) + self.m.F2(t))**2
        brace = ( (4.*xB**2*Mp2/t - 2.*xB - eps2) * FE2 + 
                   2.*xB**2*(1.-(1.-2.*xB)*t/Q2) * FM2 )
        return 8. * pt.K * (2.-y) * brace

    def cBH2unp(self, pt):
        """ BKM Eq. (37) """
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        FE2 = self.m.F1(t)**2 - t * self.m.F2(t)**2 / (4.0 * Mp2) 
        FM2 = (self.m.F1(t) + self.m.F2(t))**2
        brace = 4.*Mp2/t * FE2 + 2. * FM2
        return 8. * xB**2 * pt.K2 * brace

    def TBH2unp(self, pt):
        """ unp Bethe-Heitler amplitude squared. BKM Eq. (25)  """
        return  self.PreFacBH(pt) * ( self.cBH0unp(pt) + 
                   self.cBH1unp(pt)*cos(pt.phi) + self.cBH2unp(pt)*cos(2.*pt.phi) )

    ###### Transversely polarized target 

    def cBH0TP(self, pt):
        """ BKM Eq. (40) """
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        F1, F2 = self.m.F1(t), self.m.F2(t)
        sqrt1yeps = sqrt(1-y-eps2*y**2/4.)
        brace = ( xB**3*Mp2/Q2*(1-t/Q2)*(F1+F2) +
                (1-(1-xB)*t/Q2)*(xB**2*Mp2/t*(1-t/Q2)*F1+xB/2.*F2) )
        return (-8*pt.in1polarization*cos(pt.varphi)*(2-y)*y*sqrt(Q2)/Mp*
                sqrt(1+eps2)*pt.K/sqrt1yeps*(F1+F2)*brace)

    def cBH1TP(self, pt):
        """ BKM Eq. (41) """
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        F1, F2 = self.m.F1(t), self.m.F2(t)
        sqrt1yeps = sqrt(1-y-eps2*y**2/4.)
        brace = ( 2*pt.K**2*Q2/t/sqrt1yeps**2 * (xB*(1-t/Q2)*F1 +
            t/4./Mp2*F2) + (1+eps2)*xB*(1-t/Q2)*(F1+t/4./Mp2*F2) )
        return (-16*pt.in1polarization*cos(pt.varphi)*xB*y*
                sqrt1yeps*Mp/sqrt(Q2)*sqrt(1+eps2)*(F1+F2)*brace)

    def sBH1TP(self, pt):
        """ BKM Eq. (42) """
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        F1, F2 = self.m.F1(t), self.m.F2(t)
        sqrt1yeps = sqrt(1-y-eps2*y**2/4.)
        return (16*pt.in1polarization*sin(pt.varphi)*xB**2*y*
                sqrt1yeps*Mp/sqrt(Q2)*sqrt((1+eps2)**3)*
                (1-t/Q2)*(F1+F2)*(F1+t/4./Mp2*F2))

    def TBH2TP(self, pt):
        """ TP Bethe-Heitler amplitude squared. BKM Eq. (25)  """
        return  self.PreFacBH(pt) * ( self.cBH0TP(pt) + 
                   self.cBH1TP(pt)*cos(pt.phi) + self.sBH1TP(pt)*sin(pt.phi) )


    #### DVCS

    ##### {\cal C} coefficients

    def CCALDVCSunp(self, pt):
        """ BKM Eq. (66) """

        xB2 = pt.xB**2
        ReH = self.m.ReH(pt)
        ImH = self.m.ImH(pt)
        ReE = self.m.ReE(pt)
        ImE = self.m.ImE(pt)
        ReHt = self.m.ReHt(pt)
        ImHt = self.m.ImHt(pt)
        ReEt = self.m.ReEt(pt)
        ImEt = self.m.ImEt(pt)
        parenHH = ReH**2 + ImH**2 + ReHt**2 + ImHt**2
        parenEH = 2.*( ReE*ReH + ImE*ImH + ReEt*ReHt + ImEt*ImHt ) 
        parenEE =  ReE**2 + ImE**2 
        parenEtEt = ReEt**2 + ImEt**2
        brace = 4. * (1.-pt.xB) * parenHH - xB2 * parenEH - (xB2 
                + (2.-pt.xB)**2 * pt.t/(4.*Mp2)) * parenEE - xB2 * pt.t/(4.*Mp2) * parenEtEt
        return brace / (2.-pt.xB)**2


    def CCALDVCSTP(self, pt): 
        """ BKM Eq. (68) returns tuple (CTP+, Im(CTP-))"""
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        H = (self.m.ReH(pt)) + 1j * ( self.m.ImH(pt))
        EE = (self.m.ReE(pt)) + 1j * ( self.m.ImE(pt))
        tH = (self.m.ReHt(pt)) + 1j * ( self.m.ImHt(pt))
        tE = (self.m.ReEt(pt)) + 1j * ( self.m.ImEt(pt))
        HCC = (self.m.ReH(pt)) + 1j * ( -self.m.ImH(pt))
        EECC = (self.m.ReE(pt)) + 1j * ( -self.m.ImE(pt))
        tHCC = (self.m.ReHt(pt)) + 1j * ( -self.m.ImHt(pt))
        tECC = (self.m.ReEt(pt)) + 1j * ( -self.m.ImEt(pt))
        resp = ( 2*xB*(H*tECC+tE*HCC) - 2*(2-xB)*(tH*EECC+tHCC*EE) +
                xB**2*(EE*tECC+tE*EECC) ) / (2.-xB)**2
        resm = ((2-xB)*(H*EECC-EE*HCC) - xB*(tH*tECC-tE*tHCC))*2/(2.-xB)**2
        return (resp.real, resm.imag)  # resp.real = resp
           
    #### DVCS amplitude squared Fourier coefficients

    ######  Unpolarized target 

    def CDVCSunpPP(self, pt):
        """BMK Eq. (43)"""
        return 2. * (2. - 2.*pt.y + pt.y**2)


    def cDVCS0unp(self, pt):
        """ BKM Eq. (43) """
        return self.CDVCSunpPP(pt) * self.CCALDVCSunp(pt)

    def TDVCS2unp(self, pt):
        """ unp DVCS amplitude squared. BKM Eq. (26) - FIXME: only twist two now """
        return  self.PreFacDVCS(pt) * self.cDVCS0unp(pt)

    ###### Transversely polarized target 

    def cDVCS0TP(self, pt):
        """ BKM Eq. (49) """
        y = pt.y
        ReCCp, ImCCm = self.CCALDVCSTP(pt)
        bracket = ( -pt.in1polarization*y*(2-y)*cos(pt.varphi)*ReCCp +
                (2-2*y+y**2)*sin(pt.varphi)*ImCCm )
        return -sqrt(pt.Q2*pt.K2)/Mp/sqrt(1-y) * bracket

    def TDVCS2TP(self, pt):
        """ TP DVCS amplitude squared. BKM Eq. (26) - FIXME: only twist two now """
        return  self.PreFacDVCS(pt) * self.cDVCS0TP(pt)

    #### Interference

    ######  Unpolarized target 

    def ReCCALINTunp(self, pt):
        """ Real part of BKM Eq. (69) """

        return self.m.F1(pt.t)*self.m.ReH(pt) + pt.xB/(2.-pt.xB)*(self.m.F1(pt.t)+
                self.m.F2(pt.t))*self.m.ReHt(pt) - pt.t/(4.*Mp2)*self.m.F2(pt.t)*self.m.ReE(pt)

    def ImCCALINTunp(self, pt):
        """ Imag part of BKM Eq. (69) """

        return self.m.F1(pt.t)*self.m.ImH(pt) + pt.xB/(2.-pt.xB)*(self.m.F1(pt.t)+
                self.m.F2(pt.t))*self.m.ImHt(pt) - pt.t/(4.*Mp2)*self.m.F2(pt.t)*self.m.ImE(pt)

    def ReDELCCALINTunp(self, pt):
        """ Real part of BKM Eq. (72) """

        fx = pt.xB / (2. - pt.xB)
        return - fx * (self.m.F1(pt.t)+self.m.F2(pt.t)) * ( fx *(self.m.ReH(pt) 
            + self.m.ReE(pt)) + self.m.ReHt(pt) )

    def ImDELCCALINTunp(self, pt):
        """ Imag part of BKM Eq. (72) """

        fx = pt.xB / (2. - pt.xB)
        return - fx * (self.m.F1(pt.t)+self.m.F2(pt.t)) * ( fx *(self.m.ImH(pt) 
            + self.m.ImE(pt)) + self.m.ImHt(pt) )

    def ReCCALINTunpEFF(self, pt):
        return 0

    def ImCCALINTunpEFF(self, pt):
        return 0

    def cINT0unpSX(self, pt):
        """ BKM Eq. (53) - small-x approximation!! """
        return -8. * (2. - pt.y) * (2. - 2. * pt.y + pt.y**2) *  (
                -pt.t/pt.Q2) * self.ReCCALINTunp(pt)

    def cINT0unp(self, pt):
        """ BKM Eq. (53) """
        return -8. * (2. - pt.y) * (
                   (2.-pt.y)**2 * pt.K2 * self.ReCCALINTunp(pt) / (1.-pt.y) + 
                   (pt.t/pt.Q2) * (1.-pt.y) * (2.-pt.xB) * 
                      ( self.ReCCALINTunp(pt) + self.ReDELCCALINTunp(pt) ) )

    def cINT1unp(self, pt):
        """ BKM Eq. (54) """
        return -8. * pt.K * (2. - 2. * pt.y + pt.y**2) * self.ReCCALINTunp(pt)

    def sINT1unp(self, pt):
        """ BKM Eq. (54) """
        return  pt.in1polarization * 8. * pt.K * pt.y * (2.-pt.y) * self.ImCCALINTunp(pt)

    def cINT2unp(self, pt):
        """ BKM Eq. (55) """
        return  -16. * pt.K2 * (2. - pt.y) / (2.-pt.xB) * self.ReCCALINTunpEFF(pt)

    def sINT2unp(self, pt):
        """ BKM Eq. (55) """
        return  pt.in1polarization * 16. * pt.K2 * pt.y / (2.-pt.xB) * self.ImCCALINTunpEFF(pt)

    def TINTunp(self, pt):
        """ BH-DVCS interference. BKM Eq. (27) - FIXME: only twist two """
        return  - pt.in1charge * self.PreFacINT(pt) * ( self.cINT0unp(pt)  
                + self.cINT1unp(pt) * cos(pt.phi)
                #+ self.cINT2unp(pt) * cos(2.*pt.phi) 
                + self.sINT1unp(pt) * sin(pt.phi)
                #+ self.sINT2unp(pt) * sin(2.*pt.phi)
                )

    ###### Transversely polarized target 

    def ReCCALINTTPp(self, pt):
        """ Real part of BKM Eq. (71) """

        xB, t, F1, F2 = pt.xB, pt.t, self.m.F1(pt.t), self.m.F2(pt.t)
        H, E, Ht, Et = self.m.ReH(pt), self.m.ReE(pt), self.m.ReHt(pt), self.m.ReEt(pt)
        brace1 = xB**2/(2-xB)*(H+xB/2.*E) + xB*t/4./Mp2*E
        brace2 = 4*(1-xB)/(2-xB)*F2*Ht - (xB*F1+xB**2/(2-xB)*F2)*Et
        return (F1+F2)*brace1 - xB**2/(2-xB)*F1*(Ht+xB/2.*Et) + t/4./Mp2*brace2

    def ImCCALINTTPp(self, pt):
        """ Imag part of BKM Eq. (71) FIXME: code duplication """

        xB, t, F1, F2 = pt.xB, pt.t, self.m.F1(pt.t), self.m.F2(pt.t)
        H, E, Ht, Et = self.m.ImH(pt), self.m.ImE(pt), self.m.ImHt(pt), self.m.ImEt(pt)
        brace1 = xB**2/(2-xB)*(H+xB/2.*E) + xB*t/4./Mp2*E
        brace2 = 4*(1-xB)/(2-xB)*F2*Ht - (xB*F1+xB**2/(2-xB)*F2)*Et
        return (F1+F2)*brace1 - xB**2/(2-xB)*F1*(Ht+xB/2.*Et) + t/4./Mp2*brace2
           
    def ReCCALINTTPm(self, pt):
        """ Real part of BKM Eq. (71) """

        xB, t, F1, F2 = pt.xB, pt.t, self.m.F1(pt.t), self.m.F2(pt.t)
        H, E, Ht, Et = self.m.ReH(pt), self.m.ReE(pt), self.m.ReHt(pt), self.m.ReEt(pt)
        xBaux = xB**2/(2-xB)
        paren1 = xB**2*F1 - (1-xB)*t/Mp2*F2
        brace = t/4./Mp2*((2-xB)*F1 + xBaux*F2) + xBaux*F1
        return paren1*H/(2-xB) + brace*E - xBaux*(F1+F2)*(Ht+t/4./Mp2*Et)
           
    def ImCCALINTTPm(self, pt):
        """ Imag part of BKM Eq. (71)  FIXME: code duplication"""

        xB, t, F1, F2 = pt.xB, pt.t, self.m.F1(pt.t), self.m.F2(pt.t)
        H, E, Ht, Et = self.m.ImH(pt), self.m.ImE(pt), self.m.ImHt(pt), self.m.ImEt(pt)
        xBaux = xB**2/(2-xB)
        paren1 = xB**2*F1 - (1-xB)*t/Mp2*F2
        brace = t/4./Mp2*((2-xB)*F1 + xBaux*F2) + xBaux*F1
        return paren1*H/(2-xB) + brace*E - xBaux*(F1+F2)*(Ht+t/4./Mp2*Et)

    def ReDELCCALINTTPp(self, pt):
        """ Real part of BKM Eq. (74) """

        xB, t, F1, F2 = pt.xB, pt.t, self.m.F1(pt.t), self.m.F2(pt.t)
        Ht, Et = self.m.ReHt(pt), self.m.ReEt(pt)
        return -t/Mp2*(F2*Ht - xB/(2-xB)*(F1+xB*F2/2)*Et)

    def ImDELCCALINTTPp(self, pt):
        """ Imag part of BKM Eq. (74) FIXME: code duplication """

        xB, t, F1, F2 = pt.xB, pt.t, self.m.F1(pt.t), self.m.F2(pt.t)
        Ht, Et = self.m.ImHt(pt), self.m.ImEt(pt)
        return -t/Mp2*(F2*Ht - xB/(2-xB)*(F1+xB*F2/2)*Et)

    def ReDELCCALINTTPm(self, pt):
        """ Real part of BKM Eq. (75) """

        xB, t, F1, F2 = pt.xB, pt.t, self.m.F1(pt.t), self.m.F2(pt.t)
        H, E = self.m.ReH(pt), self.m.ReE(pt)
        return t/Mp2*(F2*H - F1*E)

    def ImDELCCALINTTPm(self, pt):
        """ Imag part of BKM Eq. (75) FIXME: code duplication """

        xB, t, F1, F2 = pt.xB, pt.t, self.m.F1(pt.t), self.m.F2(pt.t)
        H, E = self.m.ImH(pt), self.m.ImE(pt)
        return t/Mp2*(F2*H - F1*E)

    def cINT0TP(self, pt):
        """ BKM Eq. (61) """
        y = pt.y
        brace1 = ((2-y)**2/(1-y)+2)*self.ReCCALINTTPp(pt) + self.ReDELCCALINTTPp(pt)
        brace2 = (2-y)**2/(1-y)*self.ImCCALINTTPm(pt) + self.ImDELCCALINTTPm(pt)
        return 8*Mp*sqrt((1-y)*pt.K2/pt.Q2) * (
            -pt.in1polarization*y*cos(pt.varphi)*brace1
            +(2-y)*sin(pt.varphi)*brace2    )

    def cINT1TP(self, pt):
        """ BKM Eq. (62) """
        y = pt.y
        brace1 = -pt.in1polarization*y*(2-y)*self.ReCCALINTTPp(pt)
        brace2 = (2-2*y+y**2)*self.ImCCALINTTPm(pt)
        return 8*Mp*sqrt((1-y)/pt.Q2) * (
            cos(pt.varphi)*brace1 + sin(pt.varphi)*brace2    )

    def sINT1TP(self, pt):
        """ BKM Eq. (62) """
        y = pt.y
        brace1 = (2-2*y+y**2)*self.ImCCALINTTPp(pt)
        brace2 = pt.in1polarization*y*(2-y)*self.ReCCALINTTPm(pt)
        return 8*Mp*sqrt((1-y)/pt.Q2) * (
            cos(pt.varphi)*brace1 + sin(pt.varphi)*brace2    )

    def TINTTP(self, pt):
        """ BH-DVCS interference. BKM Eq. (27) - FIXME: only twist two """
        return  - pt.in1charge * self.PreFacINT(pt) * ( self.cINT0TP(pt)  
                + self.cINT1TP(pt) * cos(pt.phi)
                + self.sINT1TP(pt) * sin(pt.phi)
                )

## Placeholder for original BMK longitudinally polarized target formulas

    def TBH2LP(self, pt):
        raise ValueError('XLP not implemented for BMK model! Use BM10')

    def TDVCS2LP(self, pt):
        raise ValueError('XTP not implemented for BMK model! Use BM10')

    def TINTLP(self, pt):
        raise ValueError('XTP not implemented for BMK model! Use BM10')

#####   Observables   ##

## Cross-sections

    def XS(self, pt, **kwargs):
        """Differential 5-fold e p --> e p gamma cross section. 
        
        """
        # Overriding pt kinematics with those from kwargs
        if kwargs.has_key('vars'):
            ptvars = Data.DummyPoint(init=kwargs['vars'])
            kin = utils.fill_kinematics(ptvars, old=pt)
            BMK.prepare(kin)
        else:
            # just copy everything from pt
            ptempty = Data.DummyPoint()
            kin = utils.fill_kinematics(ptempty, old=pt)
            BMK.prepare(kin)
            ## Nothing seems to be gained by the following approach:
            #kin = dict((i, getattr(pt, i)) for i in 
            #        ['xB', 'Q2', 'W', 's', 't', 'mt', 'phi', 'in1charge',
            #            'in1polarization', 'in2particle'])

        # Copy non-kinematical info
        for atr in ['in1charge', 'in1polarization', 'in2polarization']:
            if pt.has_key(atr):
                setattr(kin, atr, getattr(pt, atr))

        # For efficient calculation of XS with unpolarized beam
        if kwargs.has_key('zeropolarized') and kwargs['zeropolarized']:
            kin.in1polarization = 0

        # Flipping spins and/or charges for asymmetries
        if kwargs.has_key('flip') and kwargs['flip']:
            if isinstance(kwargs['flip'], list):
                for item in kwargs['flip']:
                    setattr(kin, item, - getattr(pt, item))
            else:
                setattr(kin, kwargs['flip'], - getattr(pt, kwargs['flip']))

        # Weighting the integrand by BH propagators
        if kwargs.has_key('weighted') and kwargs['weighted']:
            wgh = self.w(kin)
        else:
            wgh = 1

        # Gepard may need resetting
        if self.model.__dict__.has_key('g'): 
	    self.m.g.newcall = 1
	    self.m.g.parint.pid = 1

        # Finally, we build up the cross-section
        # 1. unpolarized target part
        aux = self.TBH2unp(kin) + self.TINTunp(kin) + self.TDVCS2unp(kin)
        if hasattr(pt, 'in2polarizationvector'):
            # 2. longitudinally polarized target part
            if pt.in2polarizationvector == 'L':
                aux += kin.in2polarization*(
                        self.TBH2LP(kin) + self.TINTLP(kin) + self.TDVCS2LP(kin))
            elif pt.in2polarizationvector == 'T':
            # 3. transversally polarized target part
            # We directly take cos(varphi) or sin(varphi) terms depending if
            # varFTn is specified
                if hasattr(pt, 'varFTn'):
                    # FIXME: should deal properly with Trento/BKM differences
                    # Also, this gives result in Trento convention now,
                    # and is completely O.K. for Trento asymmetries
                    #kin.varphi = (3*pt.varFTn+1)*pi/4.  # pi for cos, -pi/2 for sin
                    # And this is for BMK
                    kin.varphi = (1-pt.varFTn)*pi/4.  # 0 for cos, pi/2 for sin
                aux += kin.in2polarization*(
                        self.TBH2TP(kin) + self.TINTTP(kin) + self.TDVCS2TP(kin))
            else:
                raise ValueError('in2polarizationvector must be either L or T!')
        return wgh * self.PreFacSigma(kin) * aux


    def Xunp(self, pt, **kwargs):
        """ Calculate 4-fold differential cross section for unpolarized target. 
        
        """
        # set target polarization to zero, but first write it down
        mem = pt.__dict__.pop('in2polarization', None)
        pt.in2polarization = 0
        res = self.XS(pt, **kwargs)
        # restore old value
        if mem: pt.in2polarization = mem
        return res
           
    def XLP(self, pt, **kwargs):
        """ Differential cross section - part for transversely polarized target."""
        pt.in2polarizationvector = 'L'
        pt.in2polarization = 1
        pol = kwargs.copy()
        pol.update({'flip':'in2polarization'})
        o =  self.XS(pt, **kwargs)
        f =  self.XS(pt, **pol)
        return (o-f)/2.

    def XTP(self, pt, **kwargs):
        """Differential cross section - part for transversely polarized target."""
        pt.in2polarizationvector = 'T'
        pt.in2polarization = 1
        pol = kwargs.copy()
        pol.update({'flip':'in2polarization'})
        o =  self.XS(pt, **kwargs)
        f =  self.XS(pt, **pol)
        return (o-f)/2.

    def F2(self, pt):
        """DIS F2 form factor."""

        self.m.g.parint.pid = 0
        self.m.g.newcall = 1
        res = self.m.DISF2(pt)
        return res

    def _XDVCStApprox(self, pt):
        """Partial DVCS cross section w.r.t. Mandelstam t.
         Approx. formula used in NPB10 paper."""

        W2 = pt.W * pt.W
        self.m.g.parint.pid = 1
        self.m.g.newcall = 1
        # Simplified formula used also in Fortran gepard code
        res = 260.5633976788416 * W2 * ( 
                (self.m.ImH(pt)**2 + self.m.ReH(pt)**2)
                - pt.t/(4.*Mp2)*(self.m.ReE(pt)**2 + self.m.ImE(pt)**2)) / (
            (W2 + pt.Q2) * (2.0 * W2 + pt.Q2)**2 )
        return res

    def _XrhotApprox(self, pt):
        """Partial DVrhoP cross section w.r.t. Mandelstam t.
         Approx. formula for small xB."""

        self.m.g.parint.pid = 2
        self.m.g.newcall = 1
        res = 112175.5 * pt.xB**2 * ( 
                self.m.ImHrho(pt)**2 + self.m.ReHrho(pt)**2) / pt.Q2**2
        return res


    def _XDVCStEx(self, pt):
        """Partial DVCS cross section w.r.t. Mandelstam t."""

        eps2 = 4. * pt.xB**2 * Mp2 / pt.Q2
        self.m.g.parint.pid = 1
        self.m.g.newcall = 1
        ImH, ReH, ImE, ReE = self.m.ImH(pt), self.m.ReH(pt), self.m.ImE(pt), self.m.ReE(pt)
        res = 65.14079453579676 * ( pt.xB**2 / pt.Q2**2 / (1-pt.xB) / (2-pt.xB)**2 /
                sqrt(1 + eps2) * (
                    4 * (1 - pt.xB) * (ImH**2 + ReH**2)
                - pt.xB**2 * (ReE**2+ImE**2 + 2*ReE*ReH + 2*ImE*ImH)
                - (2-pt.xB)**2 *pt.t/4/Mp2*(ImE**2+ReE**2) ) )
        return res

    #_XDVCSt = _XDVCStApprox
    _XDVCSt = _XDVCStEx
    _Xrhot = _XrhotApprox

    def _Xt4int(self, t, pt):
        """Same as _XDVCSt/_Xrhot but with additional variable t 
        to facilitate integration over it.
        
        """
        aux = []
        for t_single in t:
            pt.t = t_single
            if hasattr(pt, 'process') and pt.process == 'gammastarp2rho0p':
            	res = self._Xrhot(pt)
	    else:
            	res = self._XDVCSt(pt)
            del pt.t
            #if debug == 2: print "t = %s  =>  dsig/dt = %s" % (t_single, res)
            aux.append(res)
        return array(aux)


    def X(self, pt):
        """Total DVCS or DVMP cross section. """
        if pt.has_key('t') or pt.has_key('tm'):
            # partial XS w.r.t momentum transfer t
            if hasattr(pt, 'process') and pt.process == 'gammastarp2rho0p':
		return self._Xrhot(pt)
	    else:
		return self._XDVCSt(pt)

        else:
            # total XS
            if pt.has_key('tmmax'):
                tmmax = pt.tmmax
            else:
                tmmax = 1.  # default -t cuttoff in GeV^2
            res = quadrature.tquadrature(lambda t: self._Xt4int(t, pt), -tmmax, 0)
            return res


## General assymetries
## TODO: Lot of code duplication here - this should be united in one clever function

    def _phiharmonic(self, fun, pt, **kwargs):
        """Return fun evaluated for phi=pt.phi, or harmonic of fun
        corresponding to pt.FTn.
        
        """
        if pt.has_key('phi') or (kwargs.has_key('vars')
                and kwargs['vars'].has_key('phi')):
            return fun(pt, **kwargs)
        elif pt.has_key('FTn'):
            if pt.FTn < 0:
                res = quadrature.Hquadrature(lambda phi: 
                        fun(pt, vars={'phi':phi}, **kwargs) * sin(-pt.FTn*phi), 0, 2*pi)
            elif pt.FTn > 0:
                res = quadrature.Hquadrature(lambda phi: 
                        fun(pt, vars={'phi':phi}, **kwargs) * cos(pt.FTn*phi), 0, 2*pi)
            elif pt.FTn == 0:
                res = quadrature.Hquadrature(lambda phi: 
                        fun(pt, vars={'phi':phi}, **kwargs), 0, 2*pi)/2.
            else:
                raise ValueError('FTn = % is weird!' % str(pt.FTn))
            return  res / pi
        else:
            raise ValueError('[%s] has neither azimuthal angle phi nor harmonic FTn defined!' % pt)

    def _TSA(self, pt, **kwargs):
        """Target spin asymmetry (transversal or longitudinal)."""
        pol = kwargs.copy()
        pol.update({'flip':'in2polarization'})
        o =  self.XS(pt, **kwargs)
        p =  self.XS(pt, **pol)
        return (o-p)/(o+p)

    def _TTSAlong(self, pt, **kwargs):
        """Calculate target spin asymmetry (TSA).
        
        According to 1004.0177 Eq. (1.6)

        """
        bpol = kwargs.copy()
        bpol.update({'flip':'in1polarization'})
        tpol = kwargs.copy()
        tpol.update({'flip':'in2polarization'})
        both = kwargs.copy()
        both.update({'flip':['in1polarization', 'in2polarization']})
        o =  self.XS(pt, **kwargs)
        p =  self.XS(pt, **bpol)
        t =  self.XS(pt, **tpol)
        b =  self.XS(pt, **both)
        return ((p+o) - (b+t)) / ((p+o) + (b+t))

    def TSA(self, pt, **kwargs):
        """Target spin asymmetry (transversal or longitudinal) or its harmonics."""
        return self._phiharmonic(self._TSA, pt, **kwargs)

    def _BTSA(self, pt, **kwargs):
        """Calculate beam-target spin asymmetry (BTSA).
        
        According to 1004.0177 Eq. (1.8)

        """
        bpol = kwargs.copy()
        bpol.update({'flip':'in1polarization'})
        tpol = kwargs.copy()
        tpol.update({'flip':'in2polarization'})
        both = kwargs.copy()
        both.update({'flip':['in1polarization', 'in2polarization']})
        o =  self.XS(pt, **kwargs)
        p =  self.XS(pt, **bpol)
        t =  self.XS(pt, **tpol)
        b =  self.XS(pt, **both)
        return ((o+b) - (p+t)) / ((o+b) + (p+t))

    def BTSA(self, pt, **kwargs):
        """Calculate beam-target spin asymmetry or its harmonics."""
        return self._phiharmonic(self._BTSA, pt, **kwargs)

    def _CBTSA(self, pt, chargepar=-1, **kwargs):
        """Calculate charge-beam spin-target spin asymmetry (CBTSA).
        
        According to 1106.2990 Eq. (18)(chargepar=1) and (19)(chargepar=-1)

        """
        T_args = kwargs.copy()
        T_args.update({'flip':'in2polarization'})
        B_args = kwargs.copy()
        B_args.update({'flip':'in1polarization'})
        BT_args = kwargs.copy()
        BT_args.update({'flip':['in1polarization', 'in2polarization']})
        C_args = kwargs.copy()
        C_args.update({'flip':'in1charge'})
        CT_args = kwargs.copy()
        CT_args.update({'flip':['in1charge', 'in2polarization']})
        CB_args = kwargs.copy()
        CB_args.update({'flip':['in1charge', 'in1polarization']})
        CBT_args = kwargs.copy()
        CBT_args.update({'flip':['in1charge' ,'in1polarization', 'in2polarization']})
        o =  self.XS(pt, **kwargs)
        t =  self.XS(pt, **T_args)
        b =  self.XS(pt, **B_args)
        bt =  self.XS(pt, **BT_args)
        c =  self.XS(pt, **C_args)
        ct =  self.XS(pt, **CT_args)
        cb =  self.XS(pt, **CB_args)
        cbt =  self.XS(pt, **CBT_args)
        return ((o-t-b+bt) + chargepar*(c-ct-cb+cbt)) / (o+t+b+bt+c+ct+cb+cbt)

    def ALTI(self, pt, **kwargs):
        """Calculate {A_LT,I} as defined by HERMES 1106.2990 Eq. (19) or its harmonics."""
        return self._phiharmonic(self._CBTSA, pt, **kwargs)

    def ALTBHDVCS(self, pt, **kwargs):
        """Calculate {A_LT,BHDVCS} as defined by HERMES 1106.2990 Eq. (18) or its harmonics."""
        return self._phiharmonic(self._CBTSA, pt, chargepar=+1, **kwargs)

    def BSD(self, pt, **kwargs):
        """Calculate 4-fold helicity-dependent cross section measured by HALL A """

        R = kwargs.copy()
        R.update({'flip':'in1polarization'})
        return ( self.Xunp(pt, **kwargs) 
                - self.Xunp(pt, **R) ) / 2.

    def BSS(self, pt, **kwargs):
        """4-fold helicity-independent cross section measured by HALL A """
        R = kwargs.copy()
        R.update({'flip':'in1polarization'})
        return ( self.Xunp(pt, **kwargs) 
                + self.Xunp(pt, **R) ) / 2.

    def _ALUI(self, pt, **kwargs):
        """Calculate BSA as defined by HERMES 0909.3587 Eq. (2.2) """
        pol = kwargs.copy()
        pol.update({'flip':'in1polarization'})
        chg = kwargs.copy()
        chg.update({'flip':'in1charge'})
        both = kwargs.copy()
        both.update({'flip':['in1polarization', 'in1charge']})
        o =  self.XS(pt, **kwargs)
        p =  self.XS(pt, **pol)
        c =  self.XS(pt, **chg)
        b =  self.XS(pt, **both)
        return ((o-p) - (c-b)) / ((o+p) + (c+b))

    def ALUI(self, pt, **kwargs):
        """Calculate BSA as defined by HERMES 0909.3587 Eq. (2.2) or its harmonics."""
        return self._phiharmonic(self._ALUI, pt, **kwargs)

    def _ALUDVCS(self, pt, **kwargs):
        """Calculate BSA as defined by HERMES 0909.3587 Eq. (2.3) """

        pol = kwargs.copy()
        pol.update({'flip':'in1polarization'})
        chg = kwargs.copy()
        chg.update({'flip':'in1charge'})
        both = kwargs.copy()
        both.update({'flip':['in1polarization', 'in1charge']})
        o =  self.XS(pt, **kwargs)
        p =  self.XS(pt, **pol)
        c =  self.XS(pt, **chg)
        b =  self.XS(pt, **both)
        return ((o-p) + (c-b)) / ((o+p) + (c+b))

    def ALUDVCS(self, pt, **kwargs):
        """Calculate BSA as defined by HERMES 0909.3587 Eq. (2.3) or its harmonics."""
        return self._phiharmonic(self._ALUDVCS, pt, **kwargs)

    def _AUTI(self, pt, **kwargs):
        """Calculate TTSA as defined by HERMES 0802.2499 Eq. (15) """

        pol = kwargs.copy()
        pol.update({'flip':'in2polarization'})
        chg = kwargs.copy()
        chg.update({'flip':'in1charge'})
        both = kwargs.copy()
        both.update({'flip':['in2polarization', 'in1charge']})
        o =  self.XS(pt, **kwargs)
        p =  self.XS(pt, **pol)
        c =  self.XS(pt, **chg)
        b =  self.XS(pt, **both)
        return ((o-p) - (c-b)) / ((o+p) + (c+b))

    def AUTI(self, pt, **kwargs):
        """Calculate TTSA as defined by HERMES 0802.2499 Eq. (15) or its phi-harmonics."""
        return self._phiharmonic(self._AUTI, pt, **kwargs)

    def _AUTDVCS(self, pt, **kwargs):
        """Calculate TTSA as defined by HERMES 0802.2499 Eq. (14) """

        pol = kwargs.copy()
        pol.update({'flip':'in2polarization'})
        chg = kwargs.copy()
        chg.update({'flip':'in1charge'})
        both = kwargs.copy()
        both.update({'flip':['in2polarization', 'in1charge']})
        o =  self.XS(pt, **kwargs)
        p =  self.XS(pt, **pol)
        c =  self.XS(pt, **chg)
        b =  self.XS(pt, **both)
        return ((o-p) + (c-b)) / ((o+p) + (c+b))

    def AUTDVCS(self, pt, **kwargs):
        """Calculate TTSA as defined by HERMES 0802.2499 Eq. (15) or its phi-harmonics."""
        return self._phiharmonic(self._AUTDVCS, pt, **kwargs)

    def _BSA(self, pt, **kwargs):
        """Calculate beam spin asymmetry (BSA)."""
        return self.BSD(pt, **kwargs) / self.BSS(pt, **kwargs)

    def BSAold(self, pt, **kwargs):
        """Calculate beam spin asymmetry (BSA) or its harmonics."""
        if pt.has_key('phi'):
            return self._BSA(pt, **kwargs)
        elif pt.has_key('FTn') and pt.FTn == -1:
            # FIXME: faster shortcut (approximate!)
            if kwargs.has_key('vars'):
                kwargs['vars'].update({'phi':pi/2.})
            else:
                kwargs['vars'] = {'phi':pi/2.}
            return  self._BSA(pt, **kwargs) 
        else:
            raise ValueError('[%s] has neither azimuthal angle phi\
 nor harmonic FTn = -1 defined!' % pt)
        ### Exact but slower:
            #res = quadrature.Hquadrature(lambda phi: 
            #        self._BSA(pt, vars={'phi':phi}) * sin(phi), 0, 2*pi)
            #return  res / pi

    def BSAexact(self, pt, **kwargs):
        """Calculate beam spin asymmetry (BSA) or its harmonics."""
        if pt.has_key('phi'):
            return self._BSA(pt, **kwargs)
        elif pt.has_key('FTn') and pt.FTn == -1:
            res = quadrature.Hquadrature(lambda phi: 
                    self._BSA(pt, vars={'phi':phi}) * sin(phi), 0, 2*pi)
        else:
            raise ValueError('[%s] has neither azimuthal angle phi\
 nor harmonic FTn == -1 defined!' % pt)
            return  res / pi

    def BSAnew(self, pt, **kwargs):
        """Calculate beam spin asymmetry (BSA) or its harmonics."""
        res = self._phiharmonic(self._BSA, pt, **kwargs)
        return  res

    BSA = BSAold

    def _BCA(self, pt, **kwargs):
        """Calculate beam charge asymmetry (BCA)."""

        kwargs.update({'zeropolarized':True})
        R = kwargs.copy()
        R.update({'flip':'in1charge'})
        return (
           self.Xunp(pt, **kwargs) 
             - self.Xunp(pt, **R) )/(
           self.Xunp(pt, **kwargs )
             + self.Xunp(pt, **R) )
        # optimized formula (remove parts which cancel anyway)
        # return  self.TINTunp(pt, phi, 0, 1) / ( 
        #               self.TBH2unp(pt, phi) + self.TDVCS2unp(pt, phi) )

    def BCA(self, pt, **kwargs):
        """Calculate beam charge asymmetry (BCA) or its harmonics."""
        res = self._phiharmonic(self._BCA, pt, **kwargs)
        # FIXME: the following has to be dealt with during 
        # conventions translation, and not here?
        #if pt.has_key('FTn') and (pt.FTn == 1 or pt.FTn == 3):
        #    res = -res
        return  res

    def BCSD(self, pt, **kwargs):
        """4-fold beam charge-spin cross section difference measured by COMPASS """
        R = kwargs.copy()
        R.update({'flip':['in1polarization', 'in1charge']})
        return (self.Xunp(pt, **kwargs) 
                - self.Xunp(pt, **R))/2.

    def BCSS(self, pt, **kwargs):
        """4-fold beam charge-spin cross section sum measured by COMPASS. """
        R = kwargs.copy()
        R.update({'flip':['in1polarization', 'in1charge']})
        return (self.Xunp(pt, **kwargs) 
                + self.Xunp(pt, **R))/2.

    def BCSA(self, pt, **kwargs):
        """Beam charge-spin asymmetry as measured by COMPASS. """
        return  self.BCSD(pt, **kwargs) / self.BCSS(pt, **kwargs)

    def BCA0minusr1(self, pt):
        return self.BCAcos0(pt) - pt.r * self.BCAcos1(pt)

    def BSDw2C(self, pt):
        """Im(C^I) as defined by HALL A """
        return self.ImCCALINTunp(pt)

    def BSSw2C(self, pt):
        """Re(C^I) or Re(C^I + Del C^I) as defined by HALL A.

        FIXME: Although it is attributed to FTn=0, Re(C^I + Del C^I)
        is only a part of zeroth harmonic.

        """
        if pt.FTn == 0:
            return self.ReCCALINTunp(pt) + self.ReDELCCALINTunp(pt)
        elif pt.FTn == 1:
            return self.ReCCALINTunp(pt)

    def BSDw(self, pt, **kwargs):
        """Calculate weighted beam spin difference (BSD) or its harmonics."""
        kwargs['weighted'] = True
        return self._phiharmonic(self.BSD, pt, **kwargs)

    def BSSw(self, pt, **kwargs):
        """Calculate weighted beam spin sum (BSS) or its harmonics."""
        kwargs['weighted'] = True
        return self._phiharmonic(self.BSS, pt, **kwargs)

    def XwA(self, pt):
        """Ratio of first two cos harmonics of w-weighted cross section. In BMK, not Trento??"""
        b0 = quadrature.Hquadrature(lambda phi: self.BSS(pt, vars={'phi':phi}, weighted=True), 
                0, 2.0*pi) / (2.0*pi)
        b1 = quadrature.Hquadrature(lambda phi: self.BSS(pt, vars={'phi':phi}, weighted=True) * cos(phi), 
                0, 2.0*pi) / pi
        return b1/b0



class hotfixedBMK(BMK):
    """Some BMK formulas are hotfixed here according to BM"""


    def CDVCSunpPP(self, pt):
        """BMK Eq. (43) + BM hotfix"""

        return 2.*( (2. - 2.*pt.y + pt.y**2 + 
            (pt.eps2*pt.y**2)/2.) / (1. + pt.eps2) )

    def CINTunpPP0(self, pt):
        """ Obtained by saving from DM's notebook """

        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return (
           (-((8.*(2. - y)*(1. + sqrt(1. + eps2)))/(2.*(1. + eps2)**2.5)))*
         ((pt.K2*(2. - y)**2)/(1. - y - (y**2*eps2)/4.) + (t/Q2)*(2. - xB)*
           sqrt(1. + eps2)*(1. - y - (y**2*eps2)/4.)*
           (1. + (eps2 + (2.*xB*t*(2. - xB + eps2/(2.*xB) + 0.5*(-1. + sqrt(
                    1. + eps2))))/Q2)/((2. - xB)*(1. + sqrt(1. + eps2)))))
         )

    def CINTunpPP0q(self, pt):
        """ Obtained by saving from DM's notebook """
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return (
        (-((8.*(2. - y)*(1. + sqrt(1. + eps2)))/(2.*(1. + eps2)**2.5)))*(t/Q2)*
         (2. - xB)*sqrt(1. + eps2)*(1. - y - (y**2*eps2)/4.)*
         (1. + (eps2 + (2.*xB*t*(2. - xB + eps2/(2.*xB) + 
               0.5*(-1. + sqrt(1. + eps2))))/Q2)/((2. - xB)*(1. + sqrt(1. + eps2))))
         )

    def CINTunpPP1(self, pt):
        """ Obtained by saving from DM's notebook """
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return (
        (-((8.*pt.K*(2. - 2.*y + y**2 + (y**2*eps2)/2.))/(1. + eps2)**2.5))*
         (((1. - eps2 + sqrt(1. + eps2))/2.)*(1. - ((1. - 3.*xB)*t)/Q2 + 
            (xB*t*(1. + 3.*eps2 - sqrt(1. + eps2)))/
             (Q2*(1. - eps2 + sqrt(1. + eps2)))) + 
          ((2.*(1. - y - (y**2*eps2)/4.))/(2. - 2.*y + y**2 + (y**2*eps2)/2.))*
           (-((3.*eps2)/4.) + (xB*t*(1. + eps2/(4.*xB) + 
               ((1. - xB)*(-1. + sqrt(1. + eps2)))/(2.*xB)))/Q2))
           )

    def CINTunpPM3(self, pt): 
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        """ Obtained by saving from DM's notebook """
        return -((8.*Q2*pt.K**3)/(Mp2*(2. - xB)**2))

    def SINTunpPP1(self, pt):
        """ Obtained by saving from DM's notebook """
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return (
                # OLD INCORRECT: (1. + eps2)**2*((8.*pt.K*(2. - y)*y)/(1. + eps2))*
                # CORRECT: ((8.*pt.K*(2. - y)*y)/(1. + eps2))*
                ((8.*pt.K*(2. - y)*y)/(1. + eps2))*
         (1. + ((t - self.tmin(Q2, xB, eps2))*(1. - xB + 0.5*(-1. + sqrt(1. + eps2))))/
           (Q2*(1. + eps2)))
         )


    def cINT0unp(self, pt):
        """ hotfixed BKM Eq. (53) """

        return (self.CINTunpPP0(pt) * self.ReCCALINTunp(pt)
                + self.CINTunpPP0q(pt) * self.ReDELCCALINTunp(pt))

    def cINT1unp(self, pt):
        """ hotfixed BKM Eq. (54) """
        return self.CINTunpPP1(pt) * self.ReCCALINTunp(pt)

    def sINT1unp(self, pt):
        """ hotfixed BKM Eq. (54) """

        return pt.in1polarization * self.SINTunpPP1(pt) * self.ImCCALINTunp(pt)


class BM10ex(hotfixedBMK):
    """According to BM arXiv:1005.5209 [hep-ph] - exact formulas
    
    Only this Approach implements observables for longitudinally polarized target!
    """


    #### Bethe-Heitler - also longitudinally polarized target  (LP)

    def cBH0LP(self, pt):
        """ BKM Eq. (38) """
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        FE = self.m.F1(t) + t * self.m.F2(t) / (4.0 * Mp2) 
        FM = self.m.F1(t) + self.m.F2(t)
        # brackets are expressions in [..] in Eq. (39)
        bracket1 = (xB/2.)*(1.-t/Q2) - t/(4.*Mp2)
        bracket2 = ( 2 - xB - 2.*(1.-xB)**2 * t/Q2 + eps2*(1-t/Q2) -
                xB*(1-2*xB)*t**2/Q2**2 )
        bracket3 = ( (xB**2*Mp2/t) * (1+t/Q2)**2 +
                (1.-xB)*(1.+xB*t/Q2) )
        return ( 8.*pt.in1polarization*xB*(2.-y) *
                y*sqrt(1.+eps2)/(1.-t/(4.*Mp2)) * FM * (
                    0.5*bracket1*bracket2*FM +
                    (1.-(1.-xB)*t/Q2)*bracket3*FE ) )

    def cBH1LP(self, pt):
        """ BKM Eq. (39) """
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        FE = self.m.F1(t) + t * self.m.F2(t) / (4.0 * Mp2) 
        FM = self.m.F1(t) + self.m.F2(t)
        bracket1 = t/(2.*Mp2) - xB*(1.-t/Q2)
        bracket2 = ( 1.+xB-(3.-2.*xB)*(1.+xB*t/Q2) - 
                4.*xB**2*Mp2/t*(1.+t**2/Q2**2) )
        return ( -8.*pt.in1polarization*xB*y*pt.K * 
                sqrt(1.+eps2)/(1.-t/(4.*Mp2)) * FM * (
                    bracket1 * (1.-xB+xB*t/Q2) * FM + bracket2 * FE ) )

    def TBH2LP(self, pt):
        """ Bethe-Heitler amplitude squared for polarized target. BKM Eq. (25)  """
        return  self.PreFacBH(pt) * ( self.cBH0LP(pt) + self.cBH1LP(pt)*cos(pt.phi) )

    #### DVCS - also longitudinally polarized target  (LP)

    def CDVCSunpPPeff(self, pt):
        """BM10 (2.18)"""

        return 16.*pt.K2/(2.-pt.xB)**2/(1+pt.eps2)
                
    def cDVCS0unp(self, pt):
        """ BM10 (2.18)"""

        return (self.CDVCSunpPP(pt) * self.CCALDVCSunp(pt) +
                self.CDVCSunpPPeff(pt) * self.CCALDVCSunp(pt, leff=1, reff=1))

    def cDVCS1unp(self, pt):
        """ BM10 (2.19)"""

        PP = 8.*pt.K/(2.-pt.xB)/(1+pt.eps2)
        return PP*(2.-pt.y) * self.CCALDVCSunp(pt, im=0, leff=1)

    def sDVCS1unp(self, pt):
        """ BM10 (2.19)"""

        PP = 8.*pt.K/(2.-pt.xB)/(1+pt.eps2)
        return PP*(- pt.in1polarization * pt.y 
                * sqrt(1.+pt.eps2)) * self.CCALDVCSunp(pt, im=1, leff=1)

    def cDVCS0LP(self, pt):
        """ BM10 (2.20)"""

        return 2.*pt.in1polarization*pt.y*(
                2.-pt.y)/sqrt(1.+pt.eps2) * self.CCALDVCSLP(pt) 

    def cDVCS1LP(self, pt):
        """ BM10 (2.21)"""

        PP = - 8.*pt.K/(2.-pt.xB)/(1+pt.eps2)
        return PP*(- pt.in1polarization * pt.y 
                * sqrt(1.+pt.eps2)) * self.CCALDVCSLP(pt, im=0, leff=1)

    def sDVCS1LP(self, pt):
        """ BM10 (2.21)"""

        PP = - 8.*pt.K/(2.-pt.xB)/(1+pt.eps2)
        return PP*(2.-pt.y) * self.CCALDVCSLP(pt, im=1, leff=1)

    def CCALDVCSunp(self, pt, im=0, leff=0, reff=0): 
        """ BM10 (2.22), from DM's notebook """
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        if leff:
            H = (self.m.ReHeff(pt)) + 1j * ( self.m.ImHeff(pt))
            EE = (self.m.ReEeff(pt)) + 1j * ( self.m.ImEeff(pt))
            tH = (self.m.ReHteff(pt)) + 1j * ( self.m.ImHteff(pt))
            tE = (self.m.ReEteff(pt)) + 1j * ( self.m.ImEteff(pt))
        else:
            H = (self.m.ReH(pt)) + 1j * ( self.m.ImH(pt))
            EE = (self.m.ReE(pt)) + 1j * ( self.m.ImE(pt))
            tH = (self.m.ReHt(pt)) + 1j * ( self.m.ImHt(pt))
            tE = (self.m.ReEt(pt)) + 1j * ( self.m.ImEt(pt))
        if reff:
            HCC = (self.m.ReHeff(pt)) + 1j * ( -self.m.ImHeff(pt))
            EECC = (self.m.ReEeff(pt)) + 1j * ( -self.m.ImEeff(pt))
            tHCC = (self.m.ReHteff(pt)) + 1j * ( -self.m.ImHteff(pt))
            tECC = (self.m.ReEteff(pt)) + 1j * ( -self.m.ImEteff(pt))
        else:
            HCC = (self.m.ReH(pt)) + 1j * ( -self.m.ImH(pt))
            EECC = (self.m.ReE(pt)) + 1j * ( -self.m.ImE(pt))
            tHCC = (self.m.ReHt(pt)) + 1j * ( -self.m.ImHt(pt))
            tECC = (self.m.ReEt(pt)) + 1j * ( -self.m.ImEt(pt))
        res = (Q2*(Q2 + t*xB)*(4*H*HCC*(1 - xB) - 
           ((EECC*H + EE*HCC)*(Q2 + t)**2*xB**2)/(Q2*(Q2 + t*xB)) - 
           (Q2*t*tE*tECC*xB**2)/(4*Mp2*(Q2 + t*xB)) - 
           (Q2*(tECC*tH + tE*tHCC)*xB**2)/(Q2 + t*xB) + 
           4*tH*tHCC*(1 - xB + (eps2*(2*Q2 + t))/(4*(Q2 + t*xB))) + 
           EE*EECC*(-(((Q2 + t)**2*xB**2)/(Q2*(Q2 + t*xB))) - 
             (t*(Q2*(2 - xB) + t*xB)**2)/(4*Mp2*Q2*(Q2 + t*xB))))
           )/(Q2*(2 - xB) + t*xB)**2
        if im:
            return res.imag
        else:
            return res.real
 
    def CCALDVCSLP(self, pt, im=0, leff=0, reff=0): 
        """ BM10 (2.23), from DM's notebook """
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        if leff:
            H = (self.m.ReHeff(pt)) + 1j * ( self.m.ImHeff(pt))
            EE = (self.m.ReEeff(pt)) + 1j * ( self.m.ImEeff(pt))
            tH = (self.m.ReHteff(pt)) + 1j * ( self.m.ImHteff(pt))
            tE = (self.m.ReEteff(pt)) + 1j * ( self.m.ImEteff(pt))
        else:
            H = (self.m.ReH(pt)) + 1j * ( self.m.ImH(pt))
            EE = (self.m.ReE(pt)) + 1j * ( self.m.ImE(pt))
            tH = (self.m.ReHt(pt)) + 1j * ( self.m.ImHt(pt))
            tE = (self.m.ReEt(pt)) + 1j * ( self.m.ImEt(pt))
        if reff:
            HCC = (self.m.ReHeff(pt)) + 1j * ( -self.m.ImHeff(pt))
            EECC = (self.m.ReEeff(pt)) + 1j * ( -self.m.ImEeff(pt))
            tHCC = (self.m.ReHteff(pt)) + 1j * ( -self.m.ImHteff(pt))
            tECC = (self.m.ReEteff(pt)) + 1j * ( -self.m.ImEteff(pt))
        else:
            HCC = (self.m.ReH(pt)) + 1j * ( -self.m.ImH(pt))
            EECC = (self.m.ReE(pt)) + 1j * ( -self.m.ImE(pt))
            tHCC = (self.m.ReHt(pt)) + 1j * ( -self.m.ImHt(pt))
            tECC = (self.m.ReEt(pt)) + 1j * ( -self.m.ImEt(pt))
        res = (Q2*(Q2 + t*xB)*(4*H*HCC*(1 - xB) - 
           ((EECC*H + EE*HCC)*(Q2 + t)**2*xB**2)/(Q2*(Q2 + t*xB)) - 
           (Q2*t*tE*tECC*xB**2)/(4*Mp2*(Q2 + t*xB)) - 
           (Q2*(tECC*tH + tE*tHCC)*xB**2)/(Q2 + t*xB) + 
           4*tH*tHCC*(1 - xB + (eps2*(2*Q2 + t))/(4*(Q2 + t*xB))) + 
           EE*EECC*(-(((Q2 + t)**2*xB**2)/(Q2*(Q2 + t*xB))) - 
             (t*(Q2*(2 - xB) + t*xB)**2)/(4*Mp2*Q2*(Q2 + t*xB))))
            )/ (Q2*(2 - xB) + t*xB)**2
        if im:
            return res.imag
        else:
            return res.real

    def TDVCS2unp(self, pt):
        """ DVCS amplitude squared. BM10 Eq. (2.17)"""
        return  self.PreFacDVCS(pt) * (self.cDVCS0unp(pt)
                + self.cDVCS1unp(pt) * cos(pt.phi)
                + self.sDVCS1unp(pt) * sin(pt.phi) )

    def TDVCS2LP(self, pt):
        """ DVCS amplitude squared. BM10 Eq. (2.17)"""
        return  self.PreFacDVCS(pt) * (self.cDVCS0LP(pt)
                + self.cDVCS1LP(pt) * cos(pt.phi)
                + self.sDVCS1LP(pt) * sin(pt.phi) )

    #### INTERFERENCE - also longitudinally polarized target  (LP)

    def CCALINTunp(self, pt, im=0, eff=0):
        """ BM10 Eq. (2.28) ... set im=1 for imag part, eff=1 for F_eff """
    
        if eff:
            CFFH = (self.m.ReHeff(pt)) + 1j * ( self.m.ImHeff(pt))
            CFFE = (self.m.ReEeff(pt)) + 1j * ( self.m.ImEeff(pt))
            CFFHt = (self.m.ReHteff(pt)) + 1j * ( self.m.ImHteff(pt))
            CFFEt = (self.m.ReEteff(pt)) + 1j * ( self.m.ImEteff(pt))
        else:
            CFFH = (self.m.ReH(pt)) + 1j * ( self.m.ImH(pt))
            CFFE = (self.m.ReE(pt)) + 1j * ( self.m.ImE(pt))
            CFFHt = (self.m.ReHt(pt)) + 1j * ( self.m.ImHt(pt))
            CFFEt = (self.m.ReEt(pt)) + 1j * ( self.m.ImEt(pt))
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        res = ( 
    self.m.F1(t)*CFFH - (t/(4*Mp2))*self.m.F2(t)*CFFE + 
     (xB/((2 - xB) + (t*xB)/pt.Q2))*(self.m.F1(t) + self.m.F2(t))*CFFHt
    )
        if im:
            return res.imag
        else:
            return res.real
    
    def CCALINTunpV(self, pt, im=0, eff=0):
        """ BM10 Eq. (2.29) ... set im=1 for imag part, eff=1 for F_eff """
    
        if eff:
            CFFH = (self.m.ReHeff(pt)) + 1j * ( self.m.ImHeff(pt))
            CFFE = (self.m.ReEeff(pt)) + 1j * ( self.m.ImEeff(pt))
            CFFHt = (self.m.ReHteff(pt)) + 1j * ( self.m.ImHteff(pt))
            CFFEt = (self.m.ReEteff(pt)) + 1j * ( self.m.ImEteff(pt))
        else:
            CFFH = (self.m.ReH(pt)) + 1j * ( self.m.ImH(pt))
            CFFE = (self.m.ReE(pt)) + 1j * ( self.m.ImE(pt))
            CFFHt = (self.m.ReHt(pt)) + 1j * ( self.m.ImHt(pt))
            CFFEt = (self.m.ReEt(pt)) + 1j * ( self.m.ImEt(pt))
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        res = ( 
    (xB/((2 - xB) + (t*xB)/pt.Q2))*(self.m.F1(t) + self.m.F2(t))*
     (CFFH + CFFE)
    )
        if im:
            return res.imag
        else:
            return res.real
    

    def CCALINTunpA(self, pt, im=0, eff=0):
        """ ... set im=1 for imag part, eff=1 for F_eff """
    
        if eff:
            CFFH = (self.m.ReHeff(pt)) + 1j * ( self.m.ImHeff(pt))
            CFFE = (self.m.ReEeff(pt)) + 1j * ( self.m.ImEeff(pt))
            CFFHt = (self.m.ReHteff(pt)) + 1j * ( self.m.ImHteff(pt))
            CFFEt = (self.m.ReEteff(pt)) + 1j * ( self.m.ImEteff(pt))
        else:
            CFFH = (self.m.ReH(pt)) + 1j * ( self.m.ImH(pt))
            CFFE = (self.m.ReE(pt)) + 1j * ( self.m.ImE(pt))
            CFFHt = (self.m.ReHt(pt)) + 1j * ( self.m.ImHt(pt))
            CFFEt = (self.m.ReEt(pt)) + 1j * ( self.m.ImEt(pt))
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        res = ( 
    (xB/((2 - xB) + (t*xB)/pt.Q2))*(self.m.F1(t) + self.m.F2(t))*CFFHt
    )
        if im:
            return res.imag
        else:
            return res.real
    

    def CCALINTLP(self, pt, im=0, eff=0):
        """ ... set im=1 for imag part, eff=1 for F_eff """
    
        if eff:
            CFFH = (self.m.ReHeff(pt)) + 1j * ( self.m.ImHeff(pt))
            CFFE = (self.m.ReEeff(pt)) + 1j * ( self.m.ImEeff(pt))
            CFFHt = (self.m.ReHteff(pt)) + 1j * ( self.m.ImHteff(pt))
            CFFEt = (self.m.ReEteff(pt)) + 1j * ( self.m.ImEteff(pt))
        else:
            CFFH = (self.m.ReH(pt)) + 1j * ( self.m.ImH(pt))
            CFFE = (self.m.ReE(pt)) + 1j * ( self.m.ImE(pt))
            CFFHt = (self.m.ReHt(pt)) + 1j * ( self.m.ImHt(pt))
            CFFEt = (self.m.ReEt(pt)) + 1j * ( self.m.ImEt(pt))
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        res = ( 
    (xB/((2 - xB) + (t*xB)/pt.Q2))*(self.m.F1(t) + self.m.F2(t))*
      (CFFH + (xB/2)*(1 - t/pt.Q2)*
        CFFE) + (1 + (Mp2/pt.Q2)*
        ((2*xB**2)/((2 - xB) + (t*xB)/pt.Q2))*
        (3 + t/pt.Q2))*self.m.F1(t)*CFFHt - 
     (t/pt.Q2)*((xB*(1 - 2*xB))/((2 - xB) + 
        (t*xB)/pt.Q2))*self.m.F2(t)*CFFHt - 
     (xB/((2 - xB) + (t*xB)/pt.Q2))*
      ((xB/2)*(1 - t/pt.Q2)*self.m.F1(t) + (t/(4*Mp2))*self.m.F2(t))*
      CFFEt
    )
        if im:
            return res.imag
        else:
            return res.real
    

    def CCALINTLPV(self, pt, im=0, eff=0):
        """ ... set im=1 for imag part, eff=1 for F_eff """
    
        if eff:
            CFFH = (self.m.ReHeff(pt)) + 1j * ( self.m.ImHeff(pt))
            CFFE = (self.m.ReEeff(pt)) + 1j * ( self.m.ImEeff(pt))
            CFFHt = (self.m.ReHteff(pt)) + 1j * ( self.m.ImHteff(pt))
            CFFEt = (self.m.ReEteff(pt)) + 1j * ( self.m.ImEteff(pt))
        else:
            CFFH = (self.m.ReH(pt)) + 1j * ( self.m.ImH(pt))
            CFFE = (self.m.ReE(pt)) + 1j * ( self.m.ImE(pt))
            CFFHt = (self.m.ReHt(pt)) + 1j * ( self.m.ImHt(pt))
            CFFEt = (self.m.ReEt(pt)) + 1j * ( self.m.ImEt(pt))
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        res = ( 
    (xB/((2 - xB) + (t*xB)/pt.Q2))*(self.m.F1(t) + self.m.F2(t))*
     (CFFH + (xB/2)*(1 - t/pt.Q2)*
       CFFE)
    )
        if im:
            return res.imag
        else:
            return res.real
    

    def CCALINTLPA(self, pt, im=0, eff=0):
        """ ... set im=1 for imag part, eff=1 for F_eff """
    
        if eff:
            CFFH = (self.m.ReHeff(pt)) + 1j * ( self.m.ImHeff(pt))
            CFFE = (self.m.ReEeff(pt)) + 1j * ( self.m.ImEeff(pt))
            CFFHt = (self.m.ReHteff(pt)) + 1j * ( self.m.ImHteff(pt))
            CFFEt = (self.m.ReEteff(pt)) + 1j * ( self.m.ImEteff(pt))
        else:
            CFFH = (self.m.ReH(pt)) + 1j * ( self.m.ImH(pt))
            CFFE = (self.m.ReE(pt)) + 1j * ( self.m.ImE(pt))
            CFFHt = (self.m.ReHt(pt)) + 1j * ( self.m.ImHt(pt))
            CFFEt = (self.m.ReEt(pt)) + 1j * ( self.m.ImEt(pt))
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        res = ( 
    (xB/((2 - xB) + (t*xB)/pt.Q2))*(self.m.F1(t) + self.m.F2(t))*
     ((1 + (2*Mp2*xB)/pt.Q2)*CFFHt + 
      (xB/2)*CFFEt)
    )
        if im:
            return res.imag
        else:
            return res.real
    

    CINT = {} # name mapping dictionary - just for convenience
    SINT = {} # name mapping dictionary - just for convenience

    def CINTLP210(self, pt):
        """Same as CINT["LP", (-1, 1), (0)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( -((4*pt.in1polarization*y)/
      (1 + eps2)**(5/2.))*
     ((pt.tK2*(2 - y)**2*(-1 + sqrt(1 + eps2)))/pt.Q2 + 
      (1 - y - (y**2*eps2)/4)*((xB*t)/pt.Q2 - 
        (1/2.)*(1 - t/pt.Q2)*eps2)*
       (-1 + sqrt(1 + eps2) + (t*(1 - 2*xB + sqrt(1 + eps2)))/
         pt.Q2))
    )
    
    CINT['LP', (-1, 1), (0)] = CINTLP210

    def CINTLP211(self, pt):
        """Same as CINT["LP", (-1, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( -((4*pt.K*y*(2 - y)*pt.in1polarization)/
      (1 + eps2)**(5/2.))*(-1 + eps2 + sqrt(1 + eps2) - 
      (t*(-1 + eps2 + sqrt(1 + eps2) + 
         2*xB*(2 - sqrt(1 + eps2))))/pt.Q2)
    )
    
    CINT['LP', (-1, 1), (1)] = CINTLP211

    def CINTLP212(self, pt):
        """Same as CINT["LP", (-1, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
     -((2*pt.in1polarization*y*(1 - y - (y**2*eps2)/4))/
      (1 + eps2)**(5/2.))*(eps2*(1 + sqrt(1 + eps2)) - 
      (t/pt.Q2)*((t*(2*xB + eps2)*
          (-1 + 2*xB + sqrt(1 + eps2)))/pt.Q2 + 
        2*(eps2 + xB*(1 - eps2 + sqrt(1 + eps2)))))
    )
    
    CINT['LP', (-1, 1), (2)] = CINTLP212

    def CINTLP213(self, pt):
        """Same as CINT["LP", (-1, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 0
    )
    
    CINT['LP', (-1, 1), (3)] = CINTLP213

    def CINTLP010(self, pt):
        """Same as CINT["LP", (0, 1), (0)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
     ((8*sqrt(2)*pt.K*(1 - xB)*y*sqrt(1 - y - (y**2*eps2)/4)*pt.in1polarization)/
      (1 + eps2)**2)*(t/pt.Q2)
    )
    
    CINT['LP', (0, 1), (0)] = CINTLP010

    def CINTLP011(self, pt):
        """Same as CINT["LP", (0, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
     -((8*sqrt(2)*(2 - y)*y*sqrt(1 - y - (y**2*eps2)/4)*pt.in1polarization)/
      (1 + eps2)**2)*(pt.tK2/pt.Q2)
    )
    
    CINT['LP', (0, 1), (1)] = CINTLP011

    def CINTLP012(self, pt):
        """Same as CINT["LP", (0, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
     -((8*sqrt(2)*pt.K*y*sqrt(1 - y - (y**2*eps2)/4)*pt.in1polarization)/
      (1 + eps2)**2)*(1 + (xB*t)/pt.Q2)
    )
    
    CINT['LP', (0, 1), (2)] = CINTLP012

    def CINTLP013(self, pt):
        """Same as CINT["LP", (0, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 0
    )
    
    CINT['LP', (0, 1), (3)] = CINTLP013

    def CINTLP110(self, pt):
        """Same as CINT["LP", (1, 1), (0)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-((4*(1 + sqrt(1 + eps2))*y*pt.in1polarization
        )/(1 + eps2)**(5/2.)))*
     ((2 - y)**2*(pt.tK2/pt.Q2) + (1 - y - (y**2*eps2)/4)*
       ((xB*t)/pt.Q2 - (1/2.)*(1 - t/pt.Q2)*
         eps2)*(1 + ((sqrt(1 + eps2) - 1 + 2*xB)/
          (1 + sqrt(1 + eps2)))*(t/pt.Q2)))
    )
    
    CINT['LP', (1, 1), (0)] = CINTLP110

    def CINTLP111(self, pt):
        """Same as CINT["LP", (1, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-((8*pt.K*(2 - y)*y*pt.in1polarization)/
       (1 + eps2)**(5/2.)))*((1 + sqrt(1 + eps2) - 
       eps2)/2)*(1 - (1 - (2*xB*(2 + sqrt(1 + eps2)))/
         (1 - eps2 + sqrt(1 + eps2)))*(t/pt.Q2))
    )
    
    CINT['LP', (1, 1), (1)] = CINTLP111

    def CINTLP112(self, pt):
        """Same as CINT["LP", (1, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-((4*(1 - y - (y**2*eps2)/4)*y*pt.in1polarization
        )/(1 + eps2)**(5/2.)))*
     ((xB*t)/pt.Q2 - (1 - t/pt.Q2)*
       (eps2/2.))*(1 - sqrt(1 + eps2) - 
      (1 + sqrt(1 + eps2) - 2*xB)*(t/pt.Q2))
    )
    
    CINT['LP', (1, 1), (2)] = CINTLP112

    def CINTLP113(self, pt):
        """Same as CINT["LP", (1, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 0
    )
    
    CINT['LP', (1, 1), (3)] = CINTLP113

    def CINTLPA210(self, pt):
        """Same as CINT["LPA", (-1, 1), (0)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
    (4*t*xB*y*(2*(2 - y)**2*((t*(1 - xB)*(1 + (t*xB)/pt.Q2))/
          pt.Q2 + ((1 + t/pt.Q2)**2*eps2)/
          4) - (1 - (t*(1 - 2*xB))/pt.Q2)*
        (1 - y - (y**2*eps2)/4)*(1 - sqrt(1 + eps2) - 
         (t*(1 - 2*xB + sqrt(1 + eps2)))/pt.Q2))*
      pt.in1polarization)/(pt.Q2*
      (1 + eps2)**(5/2.))
    )
    
    CINT['LPA', (-1, 1), (0)] = CINTLPA210

    def CINTLPA211(self, pt):
        """Same as CINT["LPA", (-1, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-16*pt.K*t*xB*(2 - y)*y*
      (1 - (t*(1 - 2*xB))/pt.Q2)*pt.in1polarization)/
     (pt.Q2*(1 + eps2)**(5/2.))
    )
    
    CINT['LPA', (-1, 1), (1)] = CINTLPA211

    def CINTLPA212(self, pt):
        """Same as CINT["LPA", (-1, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
    (-4*t*xB*y*(1 - (t*(1 - 2*xB))/pt.Q2)*
      (1 - y - (y**2*eps2)/4)*(1 + sqrt(1 + eps2) + 
       (t*(-1 + 2*xB + sqrt(1 + eps2)))/pt.Q2)*
      pt.in1polarization)/(pt.Q2*
      (1 + eps2)**(5/2.))
    )
    
    CINT['LPA', (-1, 1), (2)] = CINTLPA212

    def CINTLPA213(self, pt):
        """Same as CINT["LPA", (-1, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 0
    )
    
    CINT['LPA', (-1, 1), (3)] = CINTLPA213

    def CINTLPA010(self, pt):
        """Same as CINT["LPA", (0, 1), (0)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-8*sqrt(2)*pt.K*t*xB*y*(1 + t/pt.Q2)*
      sqrt(1 - y - (y**2*eps2)/4)*pt.in1polarization)/
     (pt.Q2*(1 + eps2)**2)
    )
    
    CINT['LPA', (0, 1), (0)] = CINTLPA010

    def CINTLPA011(self, pt):
        """Same as CINT["LPA", (0, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 0
    )
    
    CINT['LPA', (0, 1), (1)] = CINTLPA011

    def CINTLPA012(self, pt):
        """Same as CINT["LPA", (0, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (8*sqrt(2)*pt.K*t*xB*y*(1 + t/pt.Q2)*
      sqrt(1 - y - (y**2*eps2)/4)*pt.in1polarization)/
     (pt.Q2*(1 + eps2)**2)
    )
    
    CINT['LPA', (0, 1), (2)] = CINTLPA012

    def CINTLPA013(self, pt):
        """Same as CINT["LPA", (0, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 0
    )
    
    CINT['LPA', (0, 1), (3)] = CINTLPA013

    def CINTLPA110(self, pt):
        """Same as CINT["LPA", (1, 1), (0)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (8*t*xB*y*((pt.tK2*(2 - y)**2)/pt.Q2 + 
       ((1 - (t*(1 - 2*xB))/pt.Q2)*
         (1 - y - (y**2*eps2)/4)*(1 + sqrt(1 + eps2))*
         (1 + (t*(-1 + 2*xB + sqrt(1 + eps2)))/(pt.Q2*
            (1 + sqrt(1 + eps2)))))/2)*pt.in1polarization)/
     (pt.Q2*(1 + eps2)**(5/2.))
    )
    
    CINT['LPA', (1, 1), (0)] = CINTLPA110

    def CINTLPA111(self, pt):
        """Same as CINT["LPA", (1, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (16*pt.K*t*xB*(2 - y)*y*
      (1 - (t*(1 - 2*xB))/pt.Q2)*pt.in1polarization)/
     (pt.Q2*(1 + eps2)**(5/2.))
    )
    
    CINT['LPA', (1, 1), (1)] = CINTLPA111

    def CINTLPA112(self, pt):
        """Same as CINT["LPA", (1, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (4*t*xB*y*(1 - (t*(1 - 2*xB))/pt.Q2)*
      (1 - y - (y**2*eps2)/4)*(1 - sqrt(1 + eps2) - 
       (t*(1 - 2*xB + sqrt(1 + eps2)))/pt.Q2)*pt.in1polarization
      )/(pt.Q2*(1 + eps2)**(5/2.))
    )
    
    CINT['LPA', (1, 1), (2)] = CINTLPA112

    def CINTLPA113(self, pt):
        """Same as CINT["LPA", (1, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 0
    )
    
    CINT['LPA', (1, 1), (3)] = CINTLPA113

    def CINTLPV210(self, pt):
        """Same as CINT["LPV", (-1, 1), (0)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
    (2*t*y*((2*pt.tK2*(2 - y)**2*(-1 + 2*xB + sqrt(1 + eps2)))/
        pt.Q2 + (4 - 2*xB + 3*eps2)*
        (1 - y - (y**2*eps2)/4)*
        (1 + (t*(4*(1 - xB)*xB + eps2))/(pt.Q2*
           (4 - 2*xB + 3*eps2)))*(-1 + sqrt(1 + eps2) + 
         (t*(1 - 2*xB + sqrt(1 + eps2)))/pt.Q2))*
      pt.in1polarization)/(pt.Q2*
      (1 + eps2)**(5/2.))
    )
    
    CINT['LPV', (-1, 1), (0)] = CINTLPV210

    def CINTLPV211(self, pt):
        """Same as CINT["LPV", (-1, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-4*pt.K*t*(2 - y)*y*(5 - 4*xB + 3*eps2 - 
       sqrt(1 + eps2) - (t*(1 - eps2 - 
          sqrt(1 + eps2) + 2*xB*(-4 + 4*xB + 
            sqrt(1 + eps2))))/pt.Q2)*pt.in1polarization
      )/(pt.Q2*(1 + eps2)**(5/2.))
    )
    
    CINT['LPV', (-1, 1), (1)] = CINTLPV211

    def CINTLPV212(self, pt):
        """Same as CINT["LPV", (-1, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-2*t*y*(1 - y - (y**2*eps2)/4)*
      (4 - 2*xB + 3*eps2 + (t*(4*xB - 4*xB**2 + eps2))/
        pt.Q2)*(1 + sqrt(1 + eps2) + 
       (t*(-1 + 2*xB + sqrt(1 + eps2)))/pt.Q2)*
      pt.in1polarization)/(pt.Q2*
      (1 + eps2)**(5/2.))
    )
    
    CINT['LPV', (-1, 1), (2)] = CINTLPV212

    def CINTLPV213(self, pt):
        """Same as CINT["LPV", (-1, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 0
    )
    
    CINT['LPV', (-1, 1), (3)] = CINTLPV213

    def CINTLPV010(self, pt):
        """Same as CINT["LPV", (0, 1), (0)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-8*sqrt(2)*pt.K*t*y*
      (-xB + (t*(1 - 2*xB))/pt.Q2)*
      sqrt(1 - y - (y**2*eps2)/4)*pt.in1polarization)/
     (pt.Q2*(1 + eps2)**2)
    )
    
    CINT['LPV', (0, 1), (0)] = CINTLPV010

    def CINTLPV011(self, pt):
        """Same as CINT["LPV", (0, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (8*sqrt(2)*t*pt.tK2*(2 - y)*y*
      sqrt(1 - y - (y**2*eps2)/4)*pt.in1polarization)/
     (pt.Q2**2*(1 + eps2)**2)
    )
    
    CINT['LPV', (0, 1), (1)] = CINTLPV011

    def CINTLPV012(self, pt):
        """Same as CINT["LPV", (0, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (8*sqrt(2)*pt.K*t*(1 - xB)*y*
      sqrt(1 - y - (y**2*eps2)/4)*pt.in1polarization)/
     (pt.Q2*(1 + eps2)**2)
    )
    
    CINT['LPV', (0, 1), (2)] = CINTLPV012

    def CINTLPV013(self, pt):
        """Same as CINT["LPV", (0, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 0
    )
    
    CINT['LPV', (0, 1), (3)] = CINTLPV013

    def CINTLPV110(self, pt):
        """Same as CINT["LPV", (1, 1), (0)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (4*t*y*(1 + sqrt(1 + eps2))*
      ((pt.tK2*(2 - y)**2*(1 - 2*xB + sqrt(1 + eps2)))/
        (pt.Q2*(1 + sqrt(1 + eps2))) + 
       (2 - xB + (3*eps2)/2)*(1 - y - (y**2*eps2)/4)*
        (1 + (t*(4*(1 - xB)*xB + eps2))/(pt.Q2*
           (4 - 2*xB + 3*eps2)))*
        (1 + (t*(-1 + 2*xB + sqrt(1 + eps2)))/(pt.Q2*
           (1 + sqrt(1 + eps2)))))*pt.in1polarization)/
     (pt.Q2*(1 + eps2)**(5/2.))
    )
    
    CINT['LPV', (1, 1), (0)] = CINTLPV110

    def CINTLPV111(self, pt):
        """Same as CINT["LPV", (1, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (8*pt.K*t*(2 - y)*y*(2*(1 - xB) + 
       sqrt(1 + eps2))*
      (1 - ((t - self.tmin(Q2, xB, eps2))*(1 + (1 - eps2)/sqrt(1 + eps2) - 
          2*xB*(1 + (4*(1 - xB))/sqrt(1 + eps2))))/
        (2*pt.Q2*(2*(1 - xB) + sqrt(1 + eps2))))*
      pt.in1polarization)/(pt.Q2*(1 + eps2)**2)
    )
    
    CINT['LPV', (1, 1), (1)] = CINTLPV111

    def CINTLPV112(self, pt):
        """Same as CINT["LPV", (1, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-2*t*y*(4 - 2*xB + 3*eps2)*
      (1 - y - (y**2*eps2)/4)*(1 + (t*(4*(1 - xB)*xB + eps2))/
        (pt.Q2*(4 - 2*xB + 3*eps2)))*
      (-1 + sqrt(1 + eps2) + (t*(1 - 2*xB + sqrt(1 + eps2)))/
        pt.Q2)*pt.in1polarization)/
     (pt.Q2*(1 + eps2)**(5/2.))
    )
    
    CINT['LPV', (1, 1), (2)] = CINTLPV112

    def CINTLPV113(self, pt):
        """Same as CINT["LPV", (1, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 0
    )
    
    CINT['LPV', (1, 1), (3)] = CINTLPV113

    def CINTunp210(self, pt):
        """Same as CINT["unp", (-1, 1), (0)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
    (8*(2 - y)*((pt.tK2*(2 - y)**2*(-1 + sqrt(1 + eps2)))/
        (2*pt.Q2*(1 + eps2)) + 
       (t*(t - self.tmin(Q2, xB, eps2))*xB*(1 - y - (y**2*eps2)/4)*
         (1 - xB + eps2/(2*xB) + (1 - sqrt(1 + eps2))/2))/
        (pt.Q2**2*sqrt(1 + eps2))))/
     (1 + eps2)**(3/2.)
    )
    
    CINT['unp', (-1, 1), (0)] = CINTunp210

    def CINTunp211(self, pt):
        """Same as CINT["unp", (-1, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
    (8*pt.K*((2*(1 - y - (y**2*eps2)/4)*
         ((t*(1 - (3*xB)/2 + (xB + eps2/2.)/(2*sqrt(1 + eps2))))/pt.Q2 + (1 + eps2/2. - 
            sqrt(1 + eps2))/(2*sqrt(1 + eps2))))/
        sqrt(1 + eps2) + ((2 - y)**2*(2 - sqrt(1 + eps2))*
         (-((t*xB)/pt.Q2) + ((1 - t/pt.Q2)*
            (-1 + eps2 + sqrt(1 + eps2)))/
           (2*(2 - sqrt(1 + eps2)))))/(1 + eps2)))/
     (1 + eps2)**(3/2.)
    )
    
    CINT['unp', (-1, 1), (1)] = CINTunp211

    def CINTunp212(self, pt):
        """Same as CINT["unp", (-1, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (4*(2 - y)*(1 - y - (y**2*eps2)/4)*
      (1 + sqrt(1 + eps2))*
      (eps2*(1 + (t*(xB + (t*(1 - xB))/pt.Q2 + 
            sqrt(1 + eps2)))/(pt.Q2*
           (1 + sqrt(1 + eps2)))) + 
       (t*(2 - 3*xB + (t*xB*(1 - 2*xB + (2*(1 - xB))/(1 + sqrt(
                1 + eps2))))/pt.Q2))/
        pt.Q2))/(1 + eps2)**(5/2.)
    )
    
    CINT['unp', (-1, 1), (2)] = CINTunp212

    def CINTunp213(self, pt):
        """Same as CINT["unp", (-1, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-8*pt.K*(1 - y - (y**2*eps2)/4)*
      (1 + eps2/2. + sqrt(1 + eps2))*
      (1 + (t*xB*(1 + eps2/(2*xB) + sqrt(1 + eps2)))/
        (pt.Q2*(1 + eps2/2. + sqrt(1 + eps2)))))/
     (1 + eps2)**(5/2.)
    )
    
    CINT['unp', (-1, 1), (3)] = CINTunp213

    def CINTunp010(self, pt):
        """Same as CINT["unp", (0, 1), (0)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (12*sqrt(2)*pt.K*(2 - y)*
      sqrt(1 - y - (y**2*eps2)/4)*(eps2 + 
       (t*(2 - 6*xB - eps2))/(3*pt.Q2)))/
     (1 + eps2)**(5/2.)
    )
    
    CINT['unp', (0, 1), (0)] = CINTunp010

    def CINTunp011(self, pt):
        """Same as CINT["unp", (0, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (8*sqrt(2)*sqrt(1 - y - (y**2*eps2)/4)*
      (((t - self.tmin(Q2, xB, eps2))*(2 - y)**2*(1 - xB + ((t - self.tmin(Q2, xB, eps2))*((1 - xB)*xB + 
             eps2/4))/(pt.Q2*sqrt(1 + eps2))))/
        pt.Q2 + ((1 - (t*(1 - 2*xB))/pt.Q2)*
         (1 - y - (y**2*eps2)/4)*(eps2 - 
          (2*t*xB*(1 + eps2/(2*xB)))/pt.Q2))/
        sqrt(1 + eps2)))/(1 + eps2)**2
    )
    
    CINT['unp', (0, 1), (1)] = CINTunp011

    def CINTunp012(self, pt):
        """Same as CINT["unp", (0, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-8*sqrt(2)*pt.K*(2 - y)*(1 + eps2/2.)*
      sqrt(1 - y - (y**2*eps2)/4)*
      (1 + (t*xB*(1 + eps2/(2*xB)))/(pt.Q2*
         (1 + eps2/2.))))/(1 + eps2)**(5/2.)
    )
    
    CINT['unp', (0, 1), (2)] = CINTunp012

    def CINTunp013(self, pt):
        """Same as CINT["unp", (0, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 0
    )
    
    CINT['unp', (0, 1), (3)] = CINTunp013

    def CINTunp110(self, pt):
        """Same as CINT["unp", (1, 1), (0)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-4*(2 - y)*(1 + sqrt(1 + eps2))*
      ((pt.tK2*(2 - y)**2)/(pt.Q2*sqrt(1 + eps2)) + 
       (t*(2 - xB)*(1 - y - (y**2*eps2)/4)*
         (1 + (eps2 + (2*t*xB*(2 - xB + eps2/(2*xB) + 
               (-1 + sqrt(1 + eps2))/2))/pt.Q2)/
           ((2 - xB)*(1 + sqrt(1 + eps2)))))/pt.Q2))/
     (1 + eps2)**2
    )
    
    CINT['unp', (1, 1), (0)] = CINTunp110

    def CINTunp111(self, pt):
        """Same as CINT["unp", (1, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-4*pt.K*(2 - 2*y + y**2 + (y**2*eps2)/2)*
       (1 - eps2 + sqrt(1 + eps2))*
       (1 - (t*(1 - 3*xB))/pt.Q2 + 
        (t*xB*(1 + 3*eps2 - sqrt(1 + eps2)))/
         (pt.Q2*(1 - eps2 + sqrt(1 + eps2)))))/
      (1 + eps2)**(5/2.) - (16*pt.K*(1 - y - (y**2*eps2)/4)*
       ((-3*eps2)/4 + (t*xB*(1 + eps2/(4*xB) + 
           ((1 - xB)*(-1 + sqrt(1 + eps2)))/(2*xB)))/
         pt.Q2))/(1 + eps2)**(5/2.)
    )
    
    CINT['unp', (1, 1), (1)] = CINTunp111

    def CINTunp112(self, pt):
        """Same as CINT["unp", (1, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (8*(2 - y)*(1 - y - (y**2*eps2)/4)*
      ((2*pt.tK2*eps2)/(pt.Q2*(1 + eps2 + 
          sqrt(1 + eps2))) + (t*(t - self.tmin(Q2, xB, eps2))*xB*
         (1 - xB + eps2/(2*xB) + (1 - sqrt(1 + eps2))/2))/
        pt.Q2**2))/(1 + eps2)**2
    )
    
    CINT['unp', (1, 1), (2)] = CINTunp112

    def CINTunp113(self, pt):
        """Same as CINT["unp", (1, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-8*pt.K*(1 - y - (y**2*eps2)/4)*
      (-1 + sqrt(1 + eps2))*((t*(1 - xB))/pt.Q2 + 
       ((1 + t/pt.Q2)*(-1 + sqrt(1 + eps2)))/2))/
     (1 + eps2)**(5/2.)
    )
    
    CINT['unp', (1, 1), (3)] = CINTunp113

    def CINTunpA210(self, pt):
        """Same as CINT["unpA", (-1, 1), (0)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
    (-4*t*(2 - y)*(((t - self.tmin(Q2, xB, eps2))*(1 - y - (y**2*eps2)/4)*
         (-2*xB**2 + eps2 + xB*(3 - sqrt(1 + eps2))))/
        pt.Q2 + (pt.tK2*(-4 + 2*xB*(-2 + y)**2 + 4*y + 
          y**2*(-1 + sqrt(1 + eps2) + eps2*
             sqrt(1 + eps2))))/(pt.Q2*
         sqrt(1 + eps2))))/(pt.Q2*(1 + eps2)**2)
    )
    
    CINT['unpA', (-1, 1), (0)] = CINTunpA210

    def CINTunpA211(self, pt):
        """Same as CINT["unpA", (-1, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
    (4*pt.K*t*((2 - 2*y + y**2 + (y**2*eps2)/2)*
        (5 - 4*xB + 3*eps2 - sqrt(1 + eps2) - 
         (t*(1 - eps2 - sqrt(1 + eps2) - 
            2*xB*(4 - 4*xB - sqrt(1 + eps2))))/pt.Q2) + 
       (1 - y - (y**2*eps2)/4)*(8 + 5*eps2 + 
         2*xB*(-3 + sqrt(1 + eps2)) - 
         (t*(2 - eps2 + 2*sqrt(1 + eps2) - 
            4*xB*(3 - 3*xB + sqrt(1 + eps2))))/pt.Q2)))/
     (pt.Q2*(1 + eps2)**(5/2.))
    )
    
    CINT['unpA', (-1, 1), (1)] = CINTunpA211

    def CINTunpA212(self, pt):
        """Same as CINT["unpA", (-1, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-16*t*(2 - y)*(1 - y - (y**2*eps2)/4)*
      (-((pt.tK2*(1 - 2*xB))/(pt.Q2*(1 + eps2))) + 
       ((1 - xB)*(2*xB**2 - eps2 - xB*(3 + sqrt(1 + eps2))))/
        (4*(1 - xB)*xB + eps2) - 
       ((t - self.tmin(Q2, xB, eps2))*(-2*xB**2 + eps2 + 
          xB*(3 + sqrt(1 + eps2))))/(4*pt.Q2*
         sqrt(1 + eps2))))/(pt.Q2*
      (1 + eps2)**(3/2.))
    )
    
    CINT['unpA', (-1, 1), (2)] = CINTunpA212

    def CINTunpA213(self, pt):
        """Same as CINT["unpA", (-1, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (16*pt.K*t*(1 - y - (y**2*eps2)/4)*
      (1 - xB + ((t - self.tmin(Q2, xB, eps2))*((1 - xB)*xB + eps2/4))/
        (pt.Q2*sqrt(1 + eps2))))/
     (pt.Q2*(1 + eps2)**2)
    )
    
    CINT['unpA', (-1, 1), (3)] = CINTunpA213

    def CINTunpA010(self, pt):
        """Same as CINT["unpA", (0, 1), (0)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (4*sqrt(2)*pt.K*t*(2 - y)*
      (8 - 6*xB + 5*eps2)*sqrt(1 - y - (y**2*eps2)/4)*
      (1 - (t*(2 - 12*(1 - xB)*xB - eps2))/(pt.Q2*
         (8 - 6*xB + 5*eps2))))/(pt.Q2*
      (1 + eps2)**(5/2.))
    )
    
    CINT['unpA', (0, 1), (0)] = CINTunpA010

    def CINTunpA011(self, pt):
        """Same as CINT["unpA", (0, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (8*sqrt(2)*t*sqrt(1 - y - (y**2*eps2)/4)*
      ((pt.tK2*(1 - 2*xB)*(2 - y)**2)/pt.Q2 + 
       (1 - (t*(1 - 2*xB))/pt.Q2)*
        (1 - y - (y**2*eps2)/4)*(4 - 2*xB + 3*eps2 + 
         (4*t*((1 - xB)*xB + eps2/4))/pt.Q2)))/
     (pt.Q2*(1 + eps2)**(5/2.))
    )
    
    CINT['unpA', (0, 1), (1)] = CINTunpA011

    def CINTunpA012(self, pt):
        """Same as CINT["unpA", (0, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (8*sqrt(2)*pt.K*t*(2 - y)*
      sqrt(1 - y - (y**2*eps2)/4)*(1 - xB + 
       (2*(t - self.tmin(Q2, xB, eps2))*((1 - xB)*xB + eps2/4))/(pt.Q2*
         sqrt(1 + eps2))))/(pt.Q2*(1 + eps2)**2)
    )
    
    CINT['unpA', (0, 1), (2)] = CINTunpA012

    def CINTunpA013(self, pt):
        """Same as CINT["unpA", (0, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 0
    )
    
    CINT['unpA', (0, 1), (3)] = CINTunpA013

    def CINTunpA110(self, pt):
        """Same as CINT["unpA", (1, 1), (0)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
    (8*t*(2 - y)*((pt.tK2*(2 - y)**2*(1 - 2*xB + sqrt(1 + eps2)))/
        (2*pt.Q2*sqrt(1 + eps2)) + 
       (1 - y - (y**2*eps2)/4)*((-2*pt.tK2)/pt.Q2 + 
         ((1 + sqrt(1 + eps2))*(1 - xB + sqrt(1 + eps2) + 
            (t*(-1 + sqrt(1 + eps2) + (xB*(3 - 2*xB + 
                  sqrt(1 + eps2)))/(1 + sqrt(1 + eps2))))/
             pt.Q2))/2)))/(pt.Q2*
      (1 + eps2)**2)
    )
    
    CINT['unpA', (1, 1), (0)] = CINTunpA110

    def CINTunpA111(self, pt):
        """Same as CINT["unpA", (1, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
    (-16*pt.K*t*((1 - y - (y**2*eps2)/4)*
        (1 - (t*(1 - 2*xB))/pt.Q2 + 
         ((t - self.tmin(Q2, xB, eps2))*(4*(1 - xB)*xB + eps2))/(4*pt.Q2*
           sqrt(1 + eps2))) - (2 - y)**2*(1 - xB/2 + 
         ((t - self.tmin(Q2, xB, eps2))*(4*(1 - xB)*xB + eps2))/(2*pt.Q2*
           sqrt(1 + eps2)) + ((1 - t/pt.Q2)*
           (1 - 2*xB + sqrt(1 + eps2)))/4)))/
     (pt.Q2*(1 + eps2)**2)
    )
    
    CINT['unpA', (1, 1), (1)] = CINTunpA111

    def CINTunpA112(self, pt):
        """Same as CINT["unpA", (1, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (4*t*(2 - y)*(1 - y - (y**2*eps2)/4)*
      ((4*pt.tK2*(1 - 2*xB))/(pt.Q2*sqrt(1 + eps2)) - 
       ((t - self.tmin(Q2, xB, eps2))*xB*(3 - 2*xB + eps2/xB - sqrt(1 + eps2)))/
        pt.Q2))/(pt.Q2*(1 + eps2)**2)
    )
    
    CINT['unpA', (1, 1), (2)] = CINTunpA112

    def CINTunpA113(self, pt):
        """Same as CINT["unpA", (1, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (16*pt.K*t*(t - self.tmin(Q2, xB, eps2))*((1 - xB)*xB + eps2/4)*
      (1 - y - (y**2*eps2)/4))/(pt.Q2**2*
      (1 + eps2)**(5/2.))
    )
    
    CINT['unpA', (1, 1), (3)] = CINTunpA113

    def CINTunpV210(self, pt):
        """Same as CINT["unpV", (-1, 1), (0)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
    (-4*t*xB*(2 - y)*((-2*pt.tK2*(2 - 2*y + y**2 + (y**2*eps2)/2))/
        pt.Q2 + (1 - (t*(1 - 2*xB))/pt.Q2)*
        (1 - y - (y**2*eps2)/4)*(1 + sqrt(1 + eps2))*
        ((-1 + sqrt(1 + eps2))/(1 + sqrt(1 + eps2)) + 
         (t*(1 - 2*xB + sqrt(1 + eps2)))/(pt.Q2*
           (1 + sqrt(1 + eps2))))))/(pt.Q2*
      (1 + eps2)**(5/2.))
    )
    
    CINT['unpV', (-1, 1), (0)] = CINTunpV210

    def CINTunpV211(self, pt):
        """Same as CINT["unpV", (-1, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
    (8*pt.K*t*xB*(2*(1 - (t*(1 - 2*xB))/pt.Q2)*
        (2 - 2*y + y**2 + (y**2*eps2)/2) + 
       (1 - y - (y**2*eps2)/4)*(3 - sqrt(1 + eps2) - 
         (t*(3*(1 - 2*xB) + sqrt(1 + eps2)))/pt.Q2)))/
     (pt.Q2*(1 + eps2)**(5/2.))
    )
    
    CINT['unpV', (-1, 1), (1)] = CINTunpV211

    def CINTunpV212(self, pt):
        """Same as CINT["unpV", (-1, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (4*t*xB*(2 - y)*(1 - y - (y**2*eps2)/4)*
      (1 + (4*pt.tK2)/pt.Q2 + sqrt(1 + eps2) + 
       (t*((t*(1 - 2*xB)*(1 - 2*xB - sqrt(1 + eps2)))/
           pt.Q2 - 2*(1 - xB*(2 + sqrt(1 + eps2)))))/
        pt.Q2))/(pt.Q2*(1 + eps2)**(5/2.))
    )
    
    CINT['unpV', (-1, 1), (2)] = CINTunpV212

    def CINTunpV213(self, pt):
        """Same as CINT["unpV", (-1, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (8*pt.K*t*xB*(1 - y - (y**2*eps2)/4)*
      (1 + sqrt(1 + eps2))*
      (1 - (t*(1 - 2*xB - sqrt(1 + eps2)))/(pt.Q2*
         (1 + sqrt(1 + eps2)))))/(pt.Q2*
      (1 + eps2)**(5/2.))
    )
    
    CINT['unpV', (-1, 1), (3)] = CINTunpV213

    def CINTunpV010(self, pt):
        """Same as CINT["unpV", (0, 1), (0)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (24*sqrt(2)*pt.K*t*xB*(2 - y)*
      (1 - (t*(1 - 2*xB))/pt.Q2)*
      sqrt(1 - y - (y**2*eps2)/4))/(pt.Q2*
      (1 + eps2)**(5/2.))
    )
    
    CINT['unpV', (0, 1), (0)] = CINTunpV010

    def CINTunpV011(self, pt):
        """Same as CINT["unpV", (0, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
    (16*sqrt(2)*t*xB*sqrt(1 - y - (y**2*eps2)/4)*
      ((pt.tK2*(2 - y)**2)/pt.Q2 + 
       (1 - (t*(1 - 2*xB)*(2 - (t*(1 - 2*xB))/pt.Q2))/
          pt.Q2)*(1 - y - (y**2*eps2)/4)))/
     (pt.Q2*(1 + eps2)**(5/2.))
    )
    
    CINT['unpV', (0, 1), (1)] = CINTunpV011

    def CINTunpV012(self, pt):
        """Same as CINT["unpV", (0, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (8*sqrt(2)*pt.K*t*xB*(2 - y)*
      (1 - (t*(1 - 2*xB))/pt.Q2)*
      sqrt(1 - y - (y**2*eps2)/4))/(pt.Q2*
      (1 + eps2)**(5/2.))
    )
    
    CINT['unpV', (0, 1), (2)] = CINTunpV012

    def CINTunpV013(self, pt):
        """Same as CINT["unpV", (0, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 0
    )
    
    CINT['unpV', (0, 1), (3)] = CINTunpV013

    def CINTunpV110(self, pt):
        """Same as CINT["unpV", (1, 1), (0)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
    (8*t*xB*(2 - y)*((pt.tK2*(2 - y)**2)/(pt.Q2*
         sqrt(1 + eps2)) + ((1 + t/pt.Q2)*
         (1 - y - (y**2*eps2)/4)*(1 + sqrt(1 + eps2))*
         (1 + (t*(-1 + 2*xB + sqrt(1 + eps2)))/(pt.Q2*
            (1 + sqrt(1 + eps2)))))/2))/(pt.Q2*
      (1 + eps2)**2)
    )
    
    CINT['unpV', (1, 1), (0)] = CINTunpV110

    def CINTunpV111(self, pt):
        """Same as CINT["unpV", (1, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
    (16*pt.K*t*xB*((2 - y)**2*(1 - (t*(1 - 2*xB))/pt.Q2) + 
       ((t - self.tmin(Q2, xB, eps2))*(1 - y - (y**2*eps2)/4)*(1 - 2*xB + 
          sqrt(1 + eps2)))/(2*pt.Q2)))/
     (pt.Q2*(1 + eps2)**(5/2.))
    )
    
    CINT['unpV', (1, 1), (1)] = CINTunpV111

    def CINTunpV112(self, pt):
        """Same as CINT["unpV", (1, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (8*t*xB*(2 - y)*(1 - y - (y**2*eps2)/4)*
      ((4*pt.tK2)/(pt.Q2*sqrt(1 + eps2)) + 
       ((t - self.tmin(Q2, xB, eps2))*(1 + t/pt.Q2)*(1 - 2*xB + 
          sqrt(1 + eps2)))/(2*pt.Q2)))/
     (pt.Q2*(1 + eps2)**2)
    )
    
    CINT['unpV', (1, 1), (2)] = CINTunpV112

    def CINTunpV113(self, pt):
        """Same as CINT["unpV", (1, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-8*pt.K*t*xB*(1 - y - (y**2*eps2)/4)*
      (-1 + sqrt(1 + eps2) + (t*(1 - 2*xB + sqrt(1 + eps2)))/
        pt.Q2))/(pt.Q2*(1 + eps2)**(5/2.))
    )
    
    CINT['unpV', (1, 1), (3)] = CINTunpV113

    def SINTLP211(self, pt):
        """Same as SINT["LP", (-1, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( -((4*K)/(1 + eps2)**3)*
     ((-(2 - y)**2)*(-1 - 2*eps2 + sqrt(1 + eps2) + 
        (t*(-1 + 2*xB + sqrt(1 + eps2)))/pt.Q2) + 
      (1 - y - (y**2*eps2)/4)*(-2 - eps2 + 
        2*sqrt(1 + eps2) + 
        (t*(-eps2 + 4*sqrt(1 + eps2) - 
           2*xB*(1 + sqrt(1 + eps2))))/pt.Q2))
    )
    
    SINT['LP', (-1, 1), (1)] = SINTLP211

    def SINTLP212(self, pt):
        """Same as SINT["LP", (-1, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
     -((4*(2 - y)*(1 - y - (y**2*eps2)/4))/(1 + eps2)**3)*
     (eps2*(1 + sqrt(1 + eps2)) + 
      (t*(2 + 2*sqrt(1 + eps2) + eps2*
          sqrt(1 + eps2) + xB*(-3 + eps2 - 
           3*sqrt(1 + eps2))))/pt.Q2 + 
      (t**2*(eps2 - 2*xB**2*(2 + sqrt(1 + eps2)) + 
         xB*(3 - eps2 + sqrt(1 + eps2))))/pt.Q2**2)
    )
    
    SINT['LP', (-1, 1), (2)] = SINTLP212

    def SINTLP213(self, pt):
        """Same as SINT["LP", (-1, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
     ((4*pt.K*(1 - y - (y**2*eps2)/4))/(1 + eps2)**3)*
     (2 + eps2 + 2*sqrt(1 + eps2) + 
      (t*(eps2 + 2*xB*(1 + sqrt(1 + eps2))))/
       pt.Q2)
    )
    
    SINT['LP', (-1, 1), (3)] = SINTLP213

    def SINTLP011(self, pt):
        """Same as SINT["LP", (0, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
     ((8*sqrt(2)*sqrt(1 - y - (y**2*eps2)/4))/
      (1 + eps2)**(5/2.))*((pt.tK2*(2 - y)**2)/pt.Q2 + 
      (1 + t/pt.Q2)*(1 - y - (y**2*eps2)/4)*
       ((2*xB*t)/pt.Q2 - (1 - t/pt.Q2)*
         eps2))
    )
    
    SINT['LP', (0, 1), (1)] = SINTLP011

    def SINTLP012(self, pt):
        """Same as SINT["LP", (0, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
     ((8*sqrt(2)*pt.K*(2 - y)*sqrt(1 - y - (y**2*eps2)/4))/
      (1 + eps2)**(5/2.))*(1 + (xB*t)/pt.Q2)
    )
    
    SINT['LP', (0, 1), (2)] = SINTLP012

    def SINTLP013(self, pt):
        """Same as SINT["LP", (0, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 0
    )
    
    SINT['LP', (0, 1), (3)] = SINTLP013

    def SINTLP111(self, pt):
        """Same as SINT["LP", (1, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
    ((8*pt.K*(2 - 2*y + y**2 + (y**2*eps2)/2))/
       (1 + eps2)**3)*((1 + sqrt(1 + eps2))/2)*
      (2*sqrt(1 + eps2) - 1 + (t/pt.Q2)*
        ((1 + sqrt(1 + eps2) - 2*xB)/(1 + sqrt(1 + eps2)))) + 
     ((8*pt.K*(1 - y - (y**2*eps2)/4))/
       (1 + eps2)**3)*((3*eps2)/2 + 
       (1 - sqrt(1 + eps2) - eps2/2. - 
         xB*(3 - sqrt(1 + eps2)))*(t/pt.Q2))
    )
    
    SINT['LP', (1, 1), (1)] = SINTLP111

    def SINTLP112(self, pt):
        """Same as SINT["LP", (1, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( ((-8*(2 - y)*
       (1 - y - (y**2*eps2)/4))/(1 + eps2)**(5/2.))*
     ((2*pt.tK2)/(pt.Q2*sqrt(1 + eps2)) + 
      ((1 + sqrt(1 + eps2) - 2*xB)/2)*(1 + sqrt(1 + eps2) + 
        (xB*t)/pt.Q2)*((t - self.tmin(Q2, xB, eps2))/pt.Q2))
    )
    
    SINT['LP', (1, 1), (2)] = SINTLP112

    def SINTLP113(self, pt):
        """Same as SINT["LP", (1, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
    ((-8*pt.K*(1 - y - (y**2*eps2)/4))/
      (1 + eps2)**3)*((1 + sqrt(1 + eps2) - 2*xB)/
      (2*(sqrt(1 + eps2) + 1)))*eps2*
     ((t - self.tmin(Q2, xB, eps2))/pt.Q2)
    )
    
    SINT['LP', (1, 1), (3)] = SINTLP113

    def SINTLPA211(self, pt):
        """Same as SINT["LPA", (-1, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-8*pt.K*t*xB*(1 - y - (y**2*eps2)/4)*
       (3 - sqrt(1 + eps2) - (t*(3 - 6*xB + sqrt(1 + eps2)))/
         pt.Q2))/(pt.Q2*
       (1 + eps2)**3) - 
     (8*pt.K*t*xB*(2 - 2*y + y**2 + (y**2*eps2)/2)*
       (1 + sqrt(1 + eps2))*
       (1 + (t*(-1 + 2*xB + sqrt(1 + eps2)))/(pt.Q2*
          (1 + sqrt(1 + eps2)))))/
      (pt.Q2*(1 + eps2)**3)
    )
    
    SINT['LPA', (-1, 1), (1)] = SINTLPA211

    def SINTLPA212(self, pt):
        """Same as SINT["LPA", (-1, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-4*t*xB*(2 - y)*(1 - y - (y**2*eps2)/4)*
      (1 + (4*pt.tK2)/pt.Q2 + sqrt(1 + eps2) - 
       (t*(2 - (t*(1 - 2*xB)*(1 - 2*xB - sqrt(1 + eps2)))/
           pt.Q2 - 2*xB*(2 + sqrt(1 + eps2))))/
        pt.Q2))/(pt.Q2*
      (1 + eps2)**3)
    )
    
    SINT['LPA', (-1, 1), (2)] = SINTLPA212

    def SINTLPA213(self, pt):
        """Same as SINT["LPA", (-1, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-8*pt.K*t*xB*(1 - y - (y**2*eps2)/4)*
      (1 + sqrt(1 + eps2) - (t*(1 - 2*xB - sqrt(1 + eps2)))/
        pt.Q2))/(pt.Q2*
      (1 + eps2)**3)
    )
    
    SINT['LPA', (-1, 1), (3)] = SINTLPA213

    def SINTLPA011(self, pt):
        """Same as SINT["LPA", (0, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-16*sqrt(2)*t*xB*(1 + t/pt.Q2)*
      (1 - (t*(1 - 2*xB))/pt.Q2)*(1 - y - (y**2*eps2)/4)**
       (3/2.))/(pt.Q2*(1 + eps2)**(5/2.))
    )
    
    SINT['LPA', (0, 1), (1)] = SINTLPA011

    def SINTLPA012(self, pt):
        """Same as SINT["LPA", (0, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-8*sqrt(2)*pt.K*t*xB*(2 - y)*
      (1 + t/pt.Q2)*sqrt(1 - y - (y**2*eps2)/4)
      )/(pt.Q2*(1 + eps2)**(5/2.))
    )
    
    SINT['LPA', (0, 1), (2)] = SINTLPA012

    def SINTLPA013(self, pt):
        """Same as SINT["LPA", (0, 1), (3)] """

        return 0
    
    SINT['LPA', (0, 1), (3)] = SINTLPA013

    def SINTLPA111(self, pt):
        """Same as SINT["LPA", (1, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (8*pt.K*t*xB*(1 - y - (y**2*eps2)/4)*
       (3 + sqrt(1 + eps2))*
       (1 - (t*(3 - 6*xB - sqrt(1 + eps2)))/(pt.Q2*
          (3 + sqrt(1 + eps2)))))/
      (pt.Q2*(1 + eps2)**3) - 
     (8*pt.K*t*xB*(2 - 2*y + y**2 + (y**2*eps2)/2)*
       (-1 + sqrt(1 + eps2) + (t*(1 - 2*xB + sqrt(1 + eps2)))/
         pt.Q2))/(pt.Q2*
       (1 + eps2)**3)
    )
    
    SINT['LPA', (1, 1), (1)] = SINTLPA111

    def SINTLPA112(self, pt):
        """Same as SINT["LPA", (1, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (8*t*xB*(2 - y)*(1 - y - (y**2*eps2)/4)*
      ((2*pt.tK2)/pt.Q2 - 
       ((t - self.tmin(Q2, xB, eps2))*(1 - (t*(1 - 2*xB))/pt.Q2)*
         (1 - 2*xB + sqrt(1 + eps2)))/(2*pt.Q2))
      )/(pt.Q2*(1 + eps2)**3)
    )
    
    SINT['LPA', (1, 1), (2)] = SINTLPA112

    def SINTLPA113(self, pt):
        """Same as SINT["LPA", (1, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-8*pt.K*t*(t - self.tmin(Q2, xB, eps2))*xB*
      (1 - y - (y**2*eps2)/4)*(1 - 2*xB + sqrt(1 + eps2))
      )/(pt.Q2**2*(1 + eps2)**3)
    )
    
    SINT['LPA', (1, 1), (3)] = SINTLPA113

    def SINTLPV211(self, pt):
        """Same as SINT["LPV", (-1, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
    (-4*pt.K*t*((2 - 2*y + y**2 + (y**2*eps2)/2)*(3 + 2*eps2 + 
         sqrt(1 + eps2) - 2*xB*(1 + sqrt(1 + eps2)) + 
         (t*(1 - 2*xB)*(-1 + 2*xB + sqrt(1 + eps2)))/
          pt.Q2) + (1 - y - (y**2*eps2)/4)*
        (8 + 5*eps2 + 2*xB*(-3 + sqrt(1 + eps2)) - 
         (t*(2 - eps2 + 2*sqrt(1 + eps2) - 
            4*xB*(3*(1 - xB) + sqrt(1 + eps2))))/sqrt(pt.Q2)**
           2)))/(pt.Q2*(1 + eps2)**3)
    )
    
    SINT['LPV', (-1, 1), (1)] = SINTLPV211

    def SINTLPV212(self, pt):
        """Same as SINT["LPV", (-1, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-4*t*(2 - y)*(1 - y - (y**2*eps2)/4)*
      (2 + eps2 + (4*pt.tK2*(1 - 2*xB))/(pt.Q2*
         sqrt(1 + eps2)) + 2*sqrt(1 + eps2) - 
       xB*(1 + sqrt(1 + eps2)) + 
       (t*(eps2 + xB*(3 - 2*xB + sqrt(1 + eps2))))/
        pt.Q2))/(pt.Q2*
      (1 + eps2)**(5/2.))
    )
    
    SINT['LPV', (-1, 1), (2)] = SINTLPV212

    def SINTLPV213(self, pt):
        """Same as SINT["LPV", (-1, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-16*pt.K*t*(1 - y - (y**2*eps2)/4)*
      (1 - xB + ((t - self.tmin(Q2, xB, eps2))*((1 - xB)*xB + eps2/4))/
        (pt.Q2*sqrt(1 + eps2))))/
     (pt.Q2*(1 + eps2)**(5/2.))
    )
    
    SINT['LPV', (-1, 1), (3)] = SINTLPV213

    def SINTLPV011(self, pt):
        """Same as SINT["LPV", (0, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-8*sqrt(2)*t*sqrt(1 - y - (y**2*eps2)/4)*
      ((pt.tK2*(2 - y)**2)/pt.Q2 + (1 + t/pt.Q2)*
        (1 - y - (y**2*eps2)/4)*(4 - 2*xB + 3*eps2 + 
         (t*(4*xB - 4*xB**2 + eps2))/pt.Q2))
      )/(pt.Q2*(1 + eps2)**(5/2.))
    )
    
    SINT['LPV', (0, 1), (1)] = SINTLPV011

    def SINTLPV012(self, pt):
        """Same as SINT["LPV", (0, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-8*sqrt(2)*pt.K*t*(1 - xB)*(2 - y)*
      sqrt(1 - y - (y**2*eps2)/4))/
     (pt.Q2*(1 + eps2)**(5/2.))
    )
    
    SINT['LPV', (0, 1), (2)] = SINTLPV012

    def SINTLPV013(self, pt):
        """Same as SINT["LPV", (0, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 0
    )
    
    SINT['LPV', (0, 1), (3)] = SINTLPV013

    def SINTLPV111(self, pt):
        """Same as SINT["LPV", (1, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (8*pt.K*t*(2 - 2*y + y**2 + (y**2*eps2)/2)*
       (1 - ((t - self.tmin(Q2, xB, eps2))*(1 - 2*xB)*(1 - 2*xB + sqrt(1 + eps2)))/
         (2*pt.Q2*(1 + eps2))))/
      (pt.Q2*(1 + eps2)**2) + 
     (32*pt.K*t*(1 - y - (y**2*eps2)/4)*(1 + (5*eps2)/8 - 
        (xB*(3 + sqrt(1 + eps2)))/4)*
       (1 - (t*(1 - eps2/2. - sqrt(1 + eps2) - 
           2*xB*(3*(1 - xB) - sqrt(1 + eps2))))/(pt.Q2*
          (4 + (5*eps2)/2 - xB*(3 + sqrt(1 + eps2)))))
       )/(pt.Q2*(1 + eps2)**3)
    )
    
    SINT['LPV', (1, 1), (1)] = SINTLPV111

    def SINTLPV112(self, pt):
        """Same as SINT["LPV", (1, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (4*t*(2 - y)*(1 - y - (y**2*eps2)/4)*
      ((4*pt.tK2*(1 - 2*xB))/(pt.Q2*sqrt(1 + eps2)) - 
       ((t - self.tmin(Q2, xB, eps2))*(-2*xB**2 + eps2 + 
          xB*(3 - sqrt(1 + eps2))))/pt.Q2)
      )/(pt.Q2*(1 + eps2)**(5/2.))
    )
    
    SINT['LPV', (1, 1), (2)] = SINTLPV112

    def SINTLPV113(self, pt):
        """Same as SINT["LPV", (1, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (16*pt.K*t*(t - self.tmin(Q2, xB, eps2))*((1 - xB)*xB + eps2/4)*
      (1 - y - (y**2*eps2)/4))/
     (pt.Q2**2*(1 + eps2)**3)
    )
    
    SINT['LPV', (1, 1), (3)] = SINTLPV113

    def SINTunp211(self, pt):
        """Same as SINT["unp", (-1, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (4*pt.K*(2 - y)*y*(1 + 2*eps2 - 
       sqrt(1 + eps2) - (2*t*xB*(1 + (-1 + sqrt(1 + eps2))/
           (2*xB)))/pt.Q2)*pt.in1polarization)/(1 + eps2)**2
    )
    
    SINT['unp', (-1, 1), (1)] = SINTunp211

    def SINTunp212(self, pt):
        """Same as SINT["unp", (-1, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (2*y*(1 - y - (y**2*eps2)/4)*
      (1 + sqrt(1 + eps2))*(eps2 - 
       (2*t*xB*(1 + eps2/(2*xB)))/pt.Q2)*
      (1 + (t*(-1 + 2*xB + sqrt(1 + eps2)))/(pt.Q2*
         (1 + sqrt(1 + eps2))))*pt.in1polarization)/(1 + eps2)**2
    )
    
    SINT['unp', (-1, 1), (2)] = SINTunp212

    def SINTunp213(self, pt):
        """Same as SINT["unp", (-1, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 0
    )
    
    SINT['unp', (-1, 1), (3)] = SINTunp213

    def SINTunp011(self, pt):
        """Same as SINT["unp", (0, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (8*sqrt(2)*pt.tK2*(2 - y)*y*
      sqrt(1 - y - (y**2*eps2)/4)*pt.in1polarization)/
     (pt.Q2*(1 + eps2)**2)
    )
    
    SINT['unp', (0, 1), (1)] = SINTunp011

    def SINTunp012(self, pt):
        """Same as SINT["unp", (0, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (8*sqrt(2)*pt.K*y*(1 + eps2/2.)*
      sqrt(1 - y - (y**2*eps2)/4)*
      (1 + (t*xB*(1 + eps2/(2*xB)))/(pt.Q2*
         (1 + eps2/2.)))*pt.in1polarization)/(1 + eps2)**2
    )
    
    SINT['unp', (0, 1), (2)] = SINTunp012

    def SINTunp013(self, pt):
        """Same as SINT["unp", (0, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 0
    )
    
    SINT['unp', (0, 1), (3)] = SINTunp013

    def SINTunp111(self, pt):
        """Same as SINT["unp", (1, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
    (8*pt.K*(2 - y)*y*(1 + ((t - self.tmin(Q2, xB, eps2))*(1 - xB + (-1 + sqrt(1 + eps2))/
           2))/(pt.Q2*(1 + eps2)))*pt.in1polarization)/
     (1 + eps2)
    )
    
    SINT['unp', (1, 1), (1)] = SINTunp111

    def SINTunp112(self, pt):
        """Same as SINT["unp", (1, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-4*(t - self.tmin(Q2, xB, eps2))*y*(1 - y - (y**2*eps2)/4)*
      (1 - 2*xB + sqrt(1 + eps2))*
      (-((t - self.tmin(Q2, xB, eps2))*(2*xB + eps2))/(2*pt.Q2*
         sqrt(1 + eps2)) + (eps2 - 
         xB*(-1 + sqrt(1 + eps2)))/(1 - 2*xB + 
         sqrt(1 + eps2)))*pt.in1polarization)/(pt.Q2*
      (1 + eps2)**(3/2.))
    )
    
    SINT['unp', (1, 1), (2)] = SINTunp112

    def SINTunp113(self, pt):
        """Same as SINT["unp", (1, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 0
    )
    
    SINT['unp', (1, 1), (3)] = SINTunp113

    def SINTunpA211(self, pt):
        """Same as SINT["unpA", (-1, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (4*pt.K*t*(2 - y)*y*(3 + 2*eps2 + 
       sqrt(1 + eps2) - 2*xB*(1 + sqrt(1 + eps2)) + 
       (t*(1 - 2*xB)*(-1 + 2*xB + sqrt(1 + eps2)))/
        pt.Q2)*pt.in1polarization)/(pt.Q2*
      (1 + eps2)**2)
    )
    
    SINT['unpA', (-1, 1), (1)] = SINTunpA211

    def SINTunpA212(self, pt):
        """Same as SINT["unpA", (-1, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (2*t*y*(1 - y - (y**2*eps2)/4)*
      (4 - 2*xB + 3*eps2 + (t*(4*xB - 4*xB**2 + eps2))/
        pt.Q2)*(1 + sqrt(1 + eps2) + 
       (t*(-1 + 2*xB + sqrt(1 + eps2)))/pt.Q2)*
      pt.in1polarization)/(pt.Q2*(1 + eps2)**2)
    )
    
    SINT['unpA', (-1, 1), (2)] = SINTunpA212

    def SINTunpA213(self, pt):
        """Same as SINT["unpA", (-1, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 0
    )
    
    SINT['unpA', (-1, 1), (3)] = SINTunpA213

    def SINTunpA011(self, pt):
        """Same as SINT["unpA", (0, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-8*sqrt(2)*t*pt.tK2*(1 - 2*xB)*(2 - y)*y*
      sqrt(1 - y - (y**2*eps2)/4)*pt.in1polarization)/
     (pt.Q2**2*(1 + eps2)**2)
    )
    
    SINT['unpA', (0, 1), (1)] = SINTunpA011

    def SINTunpA012(self, pt):
        """Same as SINT["unpA", (0, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-2*sqrt(2)*pt.K*t*y*
      sqrt(1 - y - (y**2*eps2)/4)*(4 - 4*xB + 2*eps2 + 
       (2*t*(4*xB - 4*xB**2 + eps2))/pt.Q2)*pt.in1polarization)/
     (pt.Q2*(1 + eps2)**2)
    )
    
    SINT['unpA', (0, 1), (2)] = SINTunpA012

    def SINTunpA013(self, pt):
        """Same as SINT["unpA", (0, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 0
    )
    
    SINT['unpA', (0, 1), (3)] = SINTunpA013

    def SINTunpA111(self, pt):
        """Same as SINT["unpA", (1, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( ((8*pt.in1polarization*pt.K*(2 - y)*y)/(1 + eps2))*
     (t/pt.Q2)*(1 - (1 - 2*xB)*((1 + sqrt(1 + eps2) - 
         2*xB)/(2*(1 + eps2)))*((t - self.tmin(Q2, xB, eps2))/pt.Q2))
    )
    
    SINT['unpA', (1, 1), (1)] = SINTunpA111

    def SINTunpA112(self, pt):
        """Same as SINT["unpA", (1, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
    (-((16*pt.in1polarization*(1 - y - (y**2*eps2)/4)*y)/(1 + eps2)**2))*
     (t/pt.Q2)*((t - self.tmin(Q2, xB, eps2))/pt.Q2)*
     (1 - xB/2 + (3*eps2)/4)*((1 + sqrt(1 + eps2) - 2*xB)/2)*
     (1 + ((4*(1 - xB)*xB + eps2)/(4 - 2*xB + 3*eps2))*
       (t/pt.Q2))
    )
    
    SINT['unpA', (1, 1), (2)] = SINTunpA112

    def SINTunpA113(self, pt):
        """Same as SINT["unpA", (1, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 0
    )
    
    SINT['unpA', (1, 1), (3)] = SINTunpA113

    def SINTunpV211(self, pt):
        """Same as SINT["unpV", (-1, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (8*pt.K*t*xB*(2 - y)*y*(1 + sqrt(1 + eps2))*
      (1 + (t*(-1 + 2*xB + sqrt(1 + eps2)))/(pt.Q2*
         (1 + sqrt(1 + eps2))))*pt.in1polarization)/(pt.Q2*
      (1 + eps2)**2)
    )
    
    SINT['unpV', (-1, 1), (1)] = SINTunpV211

    def SINTunpV212(self, pt):
        """Same as SINT["unpV", (-1, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
    (4*t*xB*y*(1 - (t*(1 - 2*xB))/pt.Q2)*
      (1 - y - (y**2*eps2)/4)*(1 + sqrt(1 + eps2))*
      (1 + (t*(-1 + 2*xB + sqrt(1 + eps2)))/(pt.Q2*
         (1 + sqrt(1 + eps2))))*pt.in1polarization)/(pt.Q2*
      (1 + eps2)**2)
    )
    
    SINT['unpV', (-1, 1), (2)] = SINTunpV212

    def SINTunpV213(self, pt):
        """Same as SINT["unpV", (-1, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 0
    )
    
    SINT['unpV', (-1, 1), (3)] = SINTunpV213

    def SINTunpV011(self, pt):
        """Same as SINT["unpV", (0, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (4*sqrt(2)*t*xB*(2 - y)*y*
      sqrt(1 - y - (y**2*eps2)/4)*
      ((4*t*(1 - xB)*(1 + (t*xB)/pt.Q2))/pt.Q2 + 
       (1 + t/pt.Q2)**2*eps2)*pt.in1polarization)/
     (pt.Q2*(1 + eps2)**2)
    )
    
    SINT['unpV', (0, 1), (1)] = SINTunpV011

    def SINTunpV012(self, pt):
        """Same as SINT["unpV", (0, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-8*sqrt(2)*pt.K*t*xB*y*
      (1 - (t*(1 - 2*xB))/pt.Q2)*
      sqrt(1 - y - (y**2*eps2)/4)*pt.in1polarization)/
     (pt.Q2*(1 + eps2)**2)
    )
    
    SINT['unpV', (0, 1), (2)] = SINTunpV012

    def SINTunpV013(self, pt):
        """Same as SINT["unpV", (0, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 0
    )
    
    SINT['unpV', (0, 1), (3)] = SINTunpV013

    def SINTunpV111(self, pt):
        """Same as SINT["unpV", (1, 1), (1)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-((8*pt.K*(2 - y)*y*xB*pt.in1polarization)/
       (1 + eps2)**2))*(t/pt.Q2)*
     (sqrt(1 + eps2) - 1 + (1 + sqrt(1 + eps2) - 2*xB)*
       (t/pt.Q2))
    )
    
    SINT['unpV', (1, 1), (1)] = SINTunpV111

    def SINTunpV112(self, pt):
        """Same as SINT["unpV", (1, 1), (2)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 
    (-((4*(1 - y - (y**2*eps2)/4)*y*xB*pt.in1polarization)/
       (1 + eps2)**2))*(t/pt.Q2)*
     (1 - (1 - 2*xB)*(t/pt.Q2))*(sqrt(1 + eps2) - 1 + 
      (1 + sqrt(1 + eps2) - 2*xB)*(t/pt.Q2))
    )
    
    SINT['unpV', (1, 1), (2)] = SINTunpV112

    def SINTunpV113(self, pt):
        """Same as SINT["unpV", (1, 1), (3)] """

        
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( 0
    )
    
    SINT['unpV', (1, 1), (3)] = SINTunpV113

    def cINTunp(self, pt, n):
        """BM10 (2.35) and (2.36)."""

        return (
         self.CINT["unp", (1, 1), n](self, pt) * self.CCALINTunp(pt) +
         self.CINT["unpV", (1, 1), n](self, pt) * self.CCALINTunpV(pt) +
         self.CINT["unpA", (1, 1), n](self, pt) * self.CCALINTunpA(pt) +
         sqrt(2.)/(2.-pt.xB+pt.xB*pt.t/pt.Q2)*pt.tK/sqrt(pt.Q2)*(
         self.CINT["unp", (0, 1), n](self, pt) * self.CCALINTunp(pt, eff=1) +
         self.CINT["unpV", (0, 1), n](self, pt) * self.CCALINTunpV(pt, eff=1) +
         self.CINT["unpA", (0, 1), n](self, pt) * self.CCALINTunpV(pt, eff=1) )
         )

    def sINTunp(self, pt, n):
        """BM10 (2.35) and (2.36)."""

        return (
         self.SINT["unp", (1, 1), n](self, pt) * self.CCALINTunp(pt, im=1) +
         self.SINT["unpV", (1, 1), n](self, pt) * self.CCALINTunpV(pt, im=1) +
         self.SINT["unpA", (1, 1), n](self, pt) * self.CCALINTunpA(pt, im=1) +
         sqrt(2.)/(2.-pt.xB+pt.xB*pt.t/pt.Q2)*pt.tK/sqrt(pt.Q2)*(
         self.SINT["unp", (0, 1), n](self, pt) * self.CCALINTunp(pt, im=1, eff=1) +
         self.SINT["unpV", (0, 1), n](self, pt) * self.CCALINTunpV(pt, im=1, eff=1) +
         self.SINT["unpA", (0, 1), n](self, pt) * self.CCALINTunpA(pt, im=1, eff=1) )
         )

    def cINT0unp(self, pt): return self.cINTunp(pt, 0)
    def cINT1unp(self, pt): return self.cINTunp(pt, 1)
    def cINT2unp(self, pt): return self.cINTunp(pt, 2)
    def cINT3unp(self, pt): return self.cINTunp(pt, 3)
    def sINT1unp(self, pt): return self.sINTunp(pt, 1)
    def sINT2unp(self, pt): return self.sINTunp(pt, 2)
    def sINT3unp(self, pt): return self.sINTunp(pt, 3)

    def cINTLP(self, pt, n):
        """BM10 (2.35) and (2.36)."""

        return (
         self.CINT["LP", (1, 1), n](self, pt) * self.CCALINTLP(pt) +
         self.CINT["LPV", (1, 1), n](self, pt) * self.CCALINTLPV(pt) +
         self.CINT["LPA", (1, 1), n](self, pt) * self.CCALINTLPA(pt) +
         sqrt(2.)/(2.-pt.xB+pt.xB*pt.t/pt.Q2)*pt.tK/sqrt(pt.Q2)*(
         self.CINT["LP", (0, 1), n](self, pt) * self.CCALINTLP(pt, eff=1) +
         self.CINT["LPV", (0, 1), n](self, pt) * self.CCALINTLPV(pt, eff=1) +
         self.CINT["LPA", (0, 1), n](self, pt) * self.CCALINTLPV(pt, eff=1) )
         )

    def sINTLP(self, pt, n):
        """BM10 (2.35) and (2.36)."""

        return (
         self.SINT["LP", (1, 1), n](self, pt) * self.CCALINTLP(pt, im=1) +
         self.SINT["LPV", (1, 1), n](self, pt) * self.CCALINTLPV(pt, im=1) +
         self.SINT["LPA", (1, 1), n](self, pt) * self.CCALINTLPA(pt, im=1) +
         sqrt(2.)/(2.-pt.xB+pt.xB*pt.t/pt.Q2)*pt.tK/sqrt(pt.Q2)*(
         self.SINT["LP", (0, 1), n](self, pt) * self.CCALINTLP(pt, im=1, eff=1) +
         self.SINT["LPV", (0, 1), n](self, pt) * self.CCALINTLPV(pt, im=1, eff=1) +
         self.SINT["LPA", (0, 1), n](self, pt) * self.CCALINTLPV(pt, im=1, eff=1) )
         )

    def cINT0LP(self, pt): return self.cINTLP(pt, 0)
    def cINT1LP(self, pt): return self.cINTLP(pt, 1)
    def cINT2LP(self, pt): return self.cINTLP(pt, 2)
    def cINT3LP(self, pt): return self.cINTLP(pt, 3)
    def sINT1LP(self, pt): return self.sINTLP(pt, 1)
    def sINT2LP(self, pt): return self.sINTLP(pt, 2)
    def sINT3LP(self, pt): return self.sINTLP(pt, 3)

    def TINTunp(self, pt):
        """ BH-DVCS interference. BM10 Eq. (2.34)"""

        return  - pt.in1charge * self.PreFacINT(pt) * ( self.cINT0unp(pt)  
                + self.cINT1unp(pt) * cos(pt.phi)
                + self.cINT2unp(pt) * cos(2.*pt.phi) 
                + self.cINT3unp(pt) * cos(3.*pt.phi) 
                + self.sINT1unp(pt) * sin(pt.phi)
                + self.sINT2unp(pt) * sin(2.*pt.phi)
                + self.sINT3unp(pt) * sin(3.*pt.phi)
                )

    def TINTLP(self, pt):
        """ BH-DVCS interference. BM10 Eq. (2.34)"""

        return  - pt.in1charge * self.PreFacINT(pt) * (self.cINT0LP(pt)  
                + self.cINT1LP(pt) * cos(pt.phi)
                + self.cINT2LP(pt) * cos(2.*pt.phi) 
                + self.cINT3LP(pt) * cos(3.*pt.phi) 
                + self.sINT1LP(pt) * sin(pt.phi)
                + self.sINT2LP(pt) * sin(2.*pt.phi)
                + self.sINT3LP(pt) * sin(3.*pt.phi)
                )

class BM10(BM10ex):
    """According to BM arXiv:1005.5209 [hep-ph]
    
    This is BM10ex, but with Q2-suppressed terms in CCAL coefs. removed
    """

    def CCALDVCSunp(self, pt, im=0, leff=0, reff=0): 
        """ BM10 (2.22), with 1/Q2 suppressed terms removed """

        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        if leff:
            H = (self.m.ReHeff(pt)) + 1j * ( self.m.ImHeff(pt))
            EE = (self.m.ReEeff(pt)) + 1j * ( self.m.ImEeff(pt))
            tH = (self.m.ReHteff(pt)) + 1j * ( self.m.ImHteff(pt))
            tE = (self.m.ReEteff(pt)) + 1j * ( self.m.ImEteff(pt))
        else:
            H = (self.m.ReH(pt)) + 1j * ( self.m.ImH(pt))
            EE = (self.m.ReE(pt)) + 1j * ( self.m.ImE(pt))
            tH = (self.m.ReHt(pt)) + 1j * ( self.m.ImHt(pt))
            tE = (self.m.ReEt(pt)) + 1j * ( self.m.ImEt(pt))
        if reff:
            HCC = (self.m.ReHeff(pt)) + 1j * ( -self.m.ImHeff(pt))
            EECC = (self.m.ReEeff(pt)) + 1j * ( -self.m.ImEeff(pt))
            tHCC = (self.m.ReHteff(pt)) + 1j * ( -self.m.ImHteff(pt))
            tECC = (self.m.ReEteff(pt)) + 1j * ( -self.m.ImEteff(pt))
        else:
            HCC = (self.m.ReH(pt)) + 1j * ( -self.m.ImH(pt))
            EECC = (self.m.ReE(pt)) + 1j * ( -self.m.ImE(pt))
            tHCC = (self.m.ReHt(pt)) + 1j * ( -self.m.ImHt(pt))
            tECC = (self.m.ReEt(pt)) + 1j * ( -self.m.ImEt(pt))
        res = (4.*(1.-xB)*(H*HCC+tH*tHCC) 
                - xB**2*(H*EECC+EE*HCC+tH*tECC+tE*tHCC)
                - (xB**2+(2.-xB)**2*t/(4.*Mp2))*EE*EECC
                - xB**2*t/(4.*Mp2)*tE*tECC) / (2.-xB)**2
        if im:
            return res.imag
        else:
            return res.real

    def CCALDVCSLP(self, pt, im=0, leff=0, reff=0): 
        """ BM10 (2.22), with 1/Q2 suppressed terms removed """

        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        if leff:
            H = (self.m.ReHeff(pt)) + 1j * ( self.m.ImHeff(pt))
            EE = (self.m.ReEeff(pt)) + 1j * ( self.m.ImEeff(pt))
            tH = (self.m.ReHteff(pt)) + 1j * ( self.m.ImHteff(pt))
            tE = (self.m.ReEteff(pt)) + 1j * ( self.m.ImEteff(pt))
        else:
            H = (self.m.ReH(pt)) + 1j * ( self.m.ImH(pt))
            EE = (self.m.ReE(pt)) + 1j * ( self.m.ImE(pt))
            tH = (self.m.ReHt(pt)) + 1j * ( self.m.ImHt(pt))
            tE = (self.m.ReEt(pt)) + 1j * ( self.m.ImEt(pt))
        if reff:
            HCC = (self.m.ReHeff(pt)) + 1j * ( -self.m.ImHeff(pt))
            EECC = (self.m.ReEeff(pt)) + 1j * ( -self.m.ImEeff(pt))
            tHCC = (self.m.ReHteff(pt)) + 1j * ( -self.m.ImHteff(pt))
            tECC = (self.m.ReEteff(pt)) + 1j * ( -self.m.ImEteff(pt))
        else:
            HCC = (self.m.ReH(pt)) + 1j * ( -self.m.ImH(pt))
            EECC = (self.m.ReE(pt)) + 1j * ( -self.m.ImE(pt))
            tHCC = (self.m.ReHt(pt)) + 1j * ( -self.m.ImHt(pt))
            tECC = (self.m.ReEt(pt)) + 1j * ( -self.m.ImEt(pt))
        res = (4.*(1.-xB)*(H*tHCC+tH*HCC) 
                - xB**2*(H*tECC+tE*HCC+tH*EECC+EE*tHCC)
                - xB*(xB**2/2.+(2.-xB)*t/(4.*Mp2))*(EE*tECC+tE*EECC)
                ) / (2.-xB)**2
        if im:
            return res.imag
        else:
            return res.real

    def CCALINTunp(self, pt, im=0, eff=0):
        """ BM10 Eq. (2.28) with 1/Q2 suppressed terms removed """
    
        if eff:
            CFFH = (self.m.ReHeff(pt)) + 1j * ( self.m.ImHeff(pt))
            CFFE = (self.m.ReEeff(pt)) + 1j * ( self.m.ImEeff(pt))
            CFFHt = (self.m.ReHteff(pt)) + 1j * ( self.m.ImHteff(pt))
            CFFEt = (self.m.ReEteff(pt)) + 1j * ( self.m.ImEteff(pt))
        else:
            CFFH = (self.m.ReH(pt)) + 1j * ( self.m.ImH(pt))
            CFFE = (self.m.ReE(pt)) + 1j * ( self.m.ImE(pt))
            CFFHt = (self.m.ReHt(pt)) + 1j * ( self.m.ImHt(pt))
            CFFEt = (self.m.ReEt(pt)) + 1j * ( self.m.ImEt(pt))
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        res = ( 
            self.m.F1(t)*CFFH - (t/(4*Mp2))*self.m.F2(t)*CFFE + 
              xB/(2 - xB)*(self.m.F1(t) + self.m.F2(t))*CFFHt
              )
        if im:
            return res.imag
        else:
            return res.real
    
    def CCALINTunpV(self, pt, im=0, eff=0):
        """ BM10 Eq. (2.29) with 1/Q2 suppressed terms removed """
    
        if eff:
            CFFH = (self.m.ReHeff(pt)) + 1j * ( self.m.ImHeff(pt))
            CFFE = (self.m.ReEeff(pt)) + 1j * ( self.m.ImEeff(pt))
            CFFHt = (self.m.ReHteff(pt)) + 1j * ( self.m.ImHteff(pt))
            CFFEt = (self.m.ReEteff(pt)) + 1j * ( self.m.ImEteff(pt))
        else:
            CFFH = (self.m.ReH(pt)) + 1j * ( self.m.ImH(pt))
            CFFE = (self.m.ReE(pt)) + 1j * ( self.m.ImE(pt))
            CFFHt = (self.m.ReHt(pt)) + 1j * ( self.m.ImHt(pt))
            CFFEt = (self.m.ReEt(pt)) + 1j * ( self.m.ImEt(pt))
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        res = ( 
        xB/(2 - xB)*(self.m.F1(t) + self.m.F2(t)) * (CFFH + CFFE)
              )
        if im:
            return res.imag
        else:
            return res.real

    def CCALINTunpA(self, pt, im=0, eff=0):
        """ BM10 ... with 1/Q2 suppressed terms removed """
    
        if eff:
            CFFH = (self.m.ReHeff(pt)) + 1j * ( self.m.ImHeff(pt))
            CFFE = (self.m.ReEeff(pt)) + 1j * ( self.m.ImEeff(pt))
            CFFHt = (self.m.ReHteff(pt)) + 1j * ( self.m.ImHteff(pt))
            CFFEt = (self.m.ReEteff(pt)) + 1j * ( self.m.ImEteff(pt))
        else:
            CFFH = (self.m.ReH(pt)) + 1j * ( self.m.ImH(pt))
            CFFE = (self.m.ReE(pt)) + 1j * ( self.m.ImE(pt))
            CFFHt = (self.m.ReHt(pt)) + 1j * ( self.m.ImHt(pt))
            CFFEt = (self.m.ReEt(pt)) + 1j * ( self.m.ImEt(pt))
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        res = ( 
    (xB/(2 - xB))*(self.m.F1(t) + self.m.F2(t))*CFFHt
    )
        if im:
            return res.imag
        else:
            return res.real

    def CCALINTLP(self, pt, im=0, eff=0):
        """ BM10 ... with 1/Q2 suppressed terms removed """
    
        if eff:
            CFFH = (self.m.ReHeff(pt)) + 1j * ( self.m.ImHeff(pt))
            CFFE = (self.m.ReEeff(pt)) + 1j * ( self.m.ImEeff(pt))
            CFFHt = (self.m.ReHteff(pt)) + 1j * ( self.m.ImHteff(pt))
            CFFEt = (self.m.ReEteff(pt)) + 1j * ( self.m.ImEteff(pt))
        else:
            CFFH = (self.m.ReH(pt)) + 1j * ( self.m.ImH(pt))
            CFFE = (self.m.ReE(pt)) + 1j * ( self.m.ImE(pt))
            CFFHt = (self.m.ReHt(pt)) + 1j * ( self.m.ImHt(pt))
            CFFEt = (self.m.ReEt(pt)) + 1j * ( self.m.ImEt(pt))
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        res = ( 
          xB/(2 - xB)*(self.m.F1(t) + self.m.F2(t))*(CFFH + (xB/2)*CFFE) 
          + self.m.F1(t)*CFFHt 
          - xB/(2 - xB)*((xB/2)*self.m.F1(t) 
                   + (t/(4*Mp2))*self.m.F2(t))*CFFEt
              )
        if im:
            return res.imag
        else:
            return res.real
    
    def CCALINTLPV(self, pt, im=0, eff=0):
        """ BM10 ... with 1/Q2 suppressed terms removed """
    
        if eff:
            CFFH = (self.m.ReHeff(pt)) + 1j * ( self.m.ImHeff(pt))
            CFFE = (self.m.ReEeff(pt)) + 1j * ( self.m.ImEeff(pt))
            CFFHt = (self.m.ReHteff(pt)) + 1j * ( self.m.ImHteff(pt))
            CFFEt = (self.m.ReEteff(pt)) + 1j * ( self.m.ImEteff(pt))
        else:
            CFFH = (self.m.ReH(pt)) + 1j * ( self.m.ImH(pt))
            CFFE = (self.m.ReE(pt)) + 1j * ( self.m.ImE(pt))
            CFFHt = (self.m.ReHt(pt)) + 1j * ( self.m.ImHt(pt))
            CFFEt = (self.m.ReEt(pt)) + 1j * ( self.m.ImEt(pt))
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        res = ( 
          xB/(2 - xB)*(self.m.F1(t) + self.m.F2(t))*(CFFH + (xB/2)*CFFE)
              )
        if im:
            return res.imag
        else:
            return res.real
    
    def CCALINTLPA(self, pt, im=0, eff=0):
        """ BM10 ... with 1/Q2 suppressed terms removed """
    
        if eff:
            CFFH = (self.m.ReHeff(pt)) + 1j * ( self.m.ImHeff(pt))
            CFFE = (self.m.ReEeff(pt)) + 1j * ( self.m.ImEeff(pt))
            CFFHt = (self.m.ReHteff(pt)) + 1j * ( self.m.ImHteff(pt))
            CFFEt = (self.m.ReEteff(pt)) + 1j * ( self.m.ImEteff(pt))
        else:
            CFFH = (self.m.ReH(pt)) + 1j * ( self.m.ImH(pt))
            CFFE = (self.m.ReE(pt)) + 1j * ( self.m.ImE(pt))
            CFFHt = (self.m.ReHt(pt)) + 1j * ( self.m.ImHt(pt))
            CFFEt = (self.m.ReEt(pt)) + 1j * ( self.m.ImEt(pt))
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        res = ( 
         xB/(2 - xB)*(self.m.F1(t) + self.m.F2(t))*(CFFHt + (xB/2)*CFFEt)
              )
        if im:
            return res.imag
        else:
            return res.real
    
# 
# class Trento(BMK):
#     """ This is just for data loading when one wants to stay in Trento frame
#     Needed for plotting.
# 
#     """
#     def to_conventions(pt):
#         pass
