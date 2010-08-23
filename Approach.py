#from IPython.Debugger import Tracer; debug_here = Tracer()

import copy, sys

from numpy import sin, cos, pi, sqrt, array, linspace
from scipy.special import gammainc

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

    def chisq(self, points, sigmas=False):
        """Return tuple (chi-square, d.o.f., probability). If the approach and model
           provide uncertainties, they are ignored - only experimental uncertainties
           are taken into account."""
        nfreepars=utils.npars(self.model)
        dof = len(points) - nfreepars
        allsigmas = [(self.predict(pt, observable=pt.yaxis) - pt.val) / pt.err for
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
        chis = []
        for val in linspace(0.5*mem, 1.5 * mem, npoints):
            self.m.parameters[parname] = val
            chi, dof, fitprob = self.chisq(points)
            print '%s  ->  %s' % (val, chi)
        self.m.parameters[parname] = mem  # restore original value


    def print_chisq(self, points, sigmas=False):
        """Pretty-print the chi-square."""
        if sigmas:
            print self.chisq(points, sigmas=True)
        print 'P(chi-square, d.o.f) = P(%1.2f, %2d) = %5.4f' % self.chisq(points)

    def predict(self, pt, error=False, **kwargs):
        """Give prediction for DataPoint pt.

        Keyword arguments:
        error - if available, produce tuple (mean, error)
        observable - string. Default is pt.yaxis
        parameters - dictionary which will update model's one

        """
        if kwargs.has_key('observable'):
            obs = kwargs['observable']
        else:
            obs = pt.yaxis

        if kwargs.has_key('parameters'):
            old = self.model.parameters.copy()
            self.model.parameters.update(kwargs['parameters'])
        #elif isinstance(self.model, Model.ComptonNeuralNets):
        #    # It is not training (which always uses 'parameters'), and
        #    # we are not asked for particular net (call would again come
        #    # with 'parameters'), so we want mean of all nets
        #    self.model.parameters['nnet'] = 'ALL'
        #    result = getattr(self, obs)(pt)
        #    if error:
        #        return (result.mean(), result.std())
        #    else:
        #        return result.mean()

        result = getattr(self, obs)(pt)

        if kwargs.has_key('parameters'):
            # restore old values
            self.model.parameters.update(old)

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
        """Transform stuff into Approach's conventions."""
        ##  --- go to BMK conventions ----
        # C1. azimutal angle phi should be in radians ...
        if pt.has_key('phi'):
            if pt.units['phi'][:3]== 'deg': # deg, degree, degrees -> radians
                pt.phi = pt.phi * pi / 180.
                pt.newunits['phi'] = 'rad'
        # C2. ... and in BMK convention. `frame` attribute is
        # obligatory for phi-dependent data.
            if pt.frame == 'Trento':  # Trento -> BMK
                pt.phi = pi - pt.phi
                pt.newframe = 'BMK'
    to_conventions = staticmethod(to_conventions)

    def from_conventions(pt):
        """Transform stuff from Approach's conventions into original data's."""
        # C1. azimutal angle phi should be in radians ...
        if pt.has_key('phi'):
            if pt.units['phi'][:3]== 'deg': # deg, degree, degrees -> radians
                pt.phi = pt.phi * pi / 180.
                pt.units['phi'] = 'rad'
        # C2. ... and in BKM convention. `frame` attribute is
        # obligatory for phi-dependent data.
            if pt.frame == 'Trento':  # Trento -> BKM
                pt.phi = pi - pt.phi
    from_conventions = staticmethod(from_conventions)

    def prepare(pt):
        """Pre-calculate GPD-independent kinamatical constants and functions."""
        pt.y = (pt.W**2 + pt.Q2 - Mp2) / (pt.s - Mp2)
        pt.eps = 2. * pt.xB * Mp / sqrt(pt.Q2)
        pt.eps2 = pt.eps**2
        if pt.has_key('t'):
            pt.J = BMK.J(pt.Q2, pt.xB, pt.t, pt.y, pt.eps2)
            pt.K2 = BMK.K2(pt.Q2, pt.xB, pt.t, pt.y, pt.eps2)
            pt.K = sqrt(pt.K2)
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
        brace3 = 2.*eps2*(1.-t/4.*Mp2) * FE2 - xB**2*(1.-t/Q2)**2 * FM2
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
        """ Bethe-Heitler amplitude squared. BKM Eq. (25)  """
        return  self.PreFacBH(pt) * ( self.cBH0unp(pt) + 
                   self.cBH1unp(pt)*cos(pt.phi) + self.cBH2unp(pt)*cos(2.*pt.phi) )


    #### DVCS

    # {\cal C} coefficients

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
           
    # DVCS amplitude squared Fourier coefficients

    def CDVCSunpPP(self, pt):
        """BMK Eq. (43)"""

        return 2. * (2. - 2.*pt.y + pt.y**2)


    def cDVCS0unp(self, pt):
        """ BKM Eq. (43) """
        return self.CDVCSunpPP(pt) * self.CCALDVCSunp(pt)

    def TDVCS2unp(self, pt):
        """ DVCS amplitude squared. BKM Eq. (26) - FIXME: only twist two now """
        return  self.PreFacDVCS(pt) * self.cDVCS0unp(pt)

    #### Interference

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
        return  8. * pt.K * pt.y * (2.-pt.y) * self.ImCCALINTunp(pt)

    def cINT2unp(self, pt):
        """ BKM Eq. (55) """
        return  -16. * pt.K2 * (2. - pt.y) / (2.-pt.xB) * self.ReCCALINTunpEFF(pt)

    def sINT2unp(self, pt):
        """ BKM Eq. (55) """
        return  16. * pt.K2 * pt.y / (2.-pt.xB) * self.ImCCALINTunpEFF(pt)

    def TINTunp(self, pt):
        """ BH-DVCS interference. BKM Eq. (27) - FIXME: only twist two """
        return  - pt.in1charge * self.PreFacINT(pt) * ( self.cINT0unp(pt)  
                + self.cINT1unp(pt) * cos(pt.phi)
                #+ self.cINT2unp(pt) * cos(2.*pt.phi) 
                + pt.in1polarization * self.sINT1unp(pt) * sin(pt.phi)
                #+ lam * self.sINT2unp(pt) * sin(2.*pt.phi)
                )
           
    #def TINTunpd(self, pt):
    #    """ BH-DVCS interference. (Normalized) part surviving after taking difference 
    #    of two lepton longitudinal polarization states.
    #    BKM Eq. (27) - FIXME: only twist two """
    #    return  - pt.in1charge * self.PreFacINT(pt) * self.sINT1unp(pt) * sin(pt.phi)

## Observables ##

    def Xunp(self, pt, **kwargs):
        """ Calculate 4-fold differential cross section for unpolarized target. 

        lam is lepton polarization \lambda .
        FIXME: Is this 'phi' bussiness below ugly?
        
        """
        if kwargs.has_key('vars'):
            ptvars = Data.DummyPoint(init=kwargs['vars'])
            kin = utils.fill_kinematics(ptvars, old=pt)
            BMK.prepare(kin)
        else:
            # just copy everything from pt
            ptempty = Data.DummyPoint()
            kin = utils.fill_kinematics(ptempty, old=pt)
            BMK.prepare(kin)
            ## Nothing seems to be gained by following approach:
            #kin = dict((i, getattr(pt, i)) for i in 
            #        ['xB', 'Q2', 'W', 's', 't', 'mt', 'phi', 'in1charge',
            #            'in1polarization', 'in2particle'])

        # copy non-kinematical info
        for atr in ['in1charge', 'in1polarization', 'in2particle']:
            if pt.has_key(atr):
                setattr(kin, atr, getattr(pt, atr))

        if kwargs.has_key('zeropolarized') and kwargs['zeropolarized']:
            kin.in1polarization = 0

        if kwargs.has_key('flip') and kwargs['flip']:
            if isinstance(kwargs['flip'], list):
                for item in kwargs['flip']:
                    setattr(kin, item, - getattr(pt, item))
            else:
                setattr(kin, kwargs['flip'], - getattr(pt, kwargs['flip']))

        if kwargs.has_key('weighted') and kwargs['weighted']:
            wgh = self.w(kin)
        else:
            wgh = 1

        # Gepard needs resetting
        if self.model.__dict__.has_key('Gepard'): self.m.g.newcall = 1

        return wgh * self.PreFacSigma(kin) * ( self.TBH2unp(kin) 
                + self.TINTunp(kin) 
                + self.TDVCS2unp(kin) )

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

    def _BSA(self, pt, **kwargs):
        """Calculate beam spin asymmetry (BSA)."""
        return self.BSD(pt, **kwargs) / self.BSS(pt, **kwargs)

    def BSA(self, pt, **kwargs):
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
        ### Exact but slower:
        #elif pt.has_key('FTn') and pt.FTn == -1:
        #    res = quadrature.Hquadrature(lambda phi: 
        #            self._BSA(pt, {'phi':phi}) * sin(phi), 0, 2*pi)
        #    return  res / pi

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

    def BCA(self, pt):
        """Calculate beam charge asymmetry (BCA) or its harmonics."""
        if pt.has_key('phi'):
            return self._BCA(pt)
        elif pt.has_key('FTn') and pt.FTn == 0:
            res = quadrature.Hquadrature(lambda phi: 
                    self._BCA(pt, vars={'phi':phi}), 0, 2.0*pi)
            return res / (2.0*pi)
        elif pt.has_key('FTn') and pt.FTn == 1:
            res = quadrature.Hquadrature(lambda phi: 
                    self._BCA(pt, vars={'phi':phi}) * cos(phi), 0, 2*pi)
            return  - res / pi

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

    def _XDVCSt4int(self, t, pt):
        """Same as XDVCSt but with additional variable t 
        to facilitate integration over it.
        
        """
        aux = []
        for t_single in t:
            pt.t = t_single
            res = self.XDVCSt(pt)
            #if debug == 2: print "t = %s  =>  dsig/dt = %s" % (t_single, res)
            aux.append(res)

        return array(aux)


    def XDVCSt(self, pt):
        """Partial DVCS cross section w.r.t. Mandelstam t."""

        W2 = pt.W * pt.W
        self.m.g.newcall = 1
        # Simplified formula used also in Fortran gepard code
        res = 260.5633976788416 * W2 * ( 
                (self.m.ImH(pt)**2 + self.m.ReH(pt)**2)
                - pt.t/(4.*Mp2)*(self.m.ReE(pt)**2 + self.m.ImE(pt)**2)) / (
            (W2 + pt.Q2) * (2.0 * W2 + pt.Q2)**2 )
        #sys.stderr.write('%s\n' % (res,))
        return res


    def XDVCS(self, pt):
        """Total DVCS cross section."""

        res = quadrature.tquadrature(lambda t: self._XDVCSt4int(t, pt), -1, 0)
        return res

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

    def BSSw(self, pt):
        """Weighted BSS as defined by HALL A.

        """
        if pt.FTn == 0:
            return quadrature.Hquadrature(lambda phi: self.BSS(pt, vars={'phi':phi}, 
                weighted=True), 0, 2.0*pi) / (2.0*pi)
        elif pt.FTn == 1:
            return quadrature.Hquadrature(lambda phi: self.BSS(pt, vars={'phi':phi},
                weighted=True) * cos(phi), 0, 2.0*pi) / pi
        elif pt.FTn == 2:
            return quadrature.Hquadrature(lambda phi: self.BSS(pt, vars={'phi':phi},
                weighted=True) * cos(2.*phi), 0, 2.0*pi) / pi

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
        return self.SINTunpPP1(pt) * self.ImCCALINTunp(pt)

