

from numpy import sin, cos, pi, sqrt

from quadrature import Hquadrature
from constants import *
#from ansatz import *

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

#     def __init__(self, calH):
#         """Choice of {\cal H} CFF is presently the only variable part."""
# 
#         self.calH = calH

    def __init__(self, ff, optimization = False):
        global F1, F2 
        global ImcffH, RecffH, ImcffE, RecffE
        global ImcffHt, RecffHt, ImcffEt, RecffEt
        F1 = ff.F1
        F2 = ff.F2
        ImcffH = ff.ImH
        RecffH = ff.ReH
        ImcffE = ff.ImE
        RecffE = ff.ReE
        ImcffHt = ff.ImHt
        RecffHt = ff.ReHt
        ImcffEt = ff.ImEt
        RecffEt = ff.ReEt
        self.optimization = optimization
        pass



class BMK(Approach):
    """Implementation of formulas from hep-ph/0112108  (BMK)"""

    ### Kinematics ###

    def tmin(self, Q2, xB, eps2):
        """BKM Eq. (31)"""
        return -Q2 * ( 2. * (1.-xB)*(1. - sqrt(1.+eps2)) + eps2 ) / (
                4. * xB * (1.-xB) + eps2 )

    def K2(self, Q2, xB, t, y, eps2):
        """BKM Eq. (30)"""
        tm = self.tmin(Q2, xB, eps2)
        brace = sqrt(1.+eps2) + (4. * xB * (1.-xB) + eps2 ) / (
                4. * (1.-xB) ) * (t - tm) / Q2
        return -(t/Q2) * (1.-xB) * (1.-y-y*y*eps2/4.) * (
                1. - tm / t ) * brace

    def J(self, Q2, xB, t, y, eps2):
        """BKM below Eq. (32)"""
        return (1.-y-y*eps2/2.) * (1. + t/Q2) - (1.-xB)*(2.-y)*t/Q2

    def r(self, Q2, xB, t, y, eps2):
        """DM's fitting notes, below Eq. (13)"""
        K = sqrt(self.K2(Q2, xB, t, y, eps2))
        brace = (2.-y)**2 * K / (1.-y) + (1./K)*(t/Q2)*(1.-y)*(2.-xB)
        return - (2.-y) / (2.-2.*y+y**2) * brace

    def P1P2(self, pt, phi):
        """ Product of Bethe-Heitler propagators, Eq(32) """
        P1 = - ( self.J(pt.Q2, pt.xB, pt.t, pt.y, pt.eps2) + 2. * 
                sqrt(self.K2(pt.Q2, pt.xB, pt.t, pt.y, pt.eps2)) * cos(phi) ) / (
                        pt.y * (1. + pt.eps2) )
        P2 = 1. + pt.t / pt.Q2  - P1
        return P1 * P2

    def anintP1P2(self, pt):
        """ Analitical integral of \int P1 P2 d\phi """
        xB, Q2, t, y, eps2, K2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2, pt.K2
        brace = ( (1 - y - (1+eps2/2.) * y**2 * eps2/2.) * (1. + t/Q2)**2 +
                  2.*K2 - (1.-xB)*(2.-y)**2 * (1. + xB*t/Q2) * t/Q2 )
        return -2. * pi * brace / (1+eps2)**2 / y**2

    def PreFacSigma(self, pt):
        """ Prefactor of 4-fold xs. Take prefactor in Eq(22) times 2\pi because of proton Phi integration
        and times y/Q2 because of dy -> dQ2. Convert to nanobarns."""
        return alpha**3 * pt.xB * pt.y**2 / (8. * pi * pt.Q2**2 * sqrt(1.+pt.eps2)) * GeV2nb

    def PreFacBH(self, pt, phi):
        """ Prefactor from Eq. (25), without e^6 """
        return 1./(pt.xB**2 * pt.y**2 * (1.+pt.eps2)**2 * pt.t * self.P1P2(pt, phi))

    def PreFacDVCS(self, pt):
        """ Prefactor from Eq. (26), without e^6 """
        return 1./(pt.y**2 * pt.Q2 )

    def PreFacINT(self, pt, phi):
        """ Prefactor from Eq. (27), without e^6 """
        return 1./(pt.xB * pt.y**3 * pt.t * self.P1P2(pt, phi))

    def w(self, pt, phi):
        """ Weight factor removing BH propagators from INT and BH amplitudes. 
        It is normalized to \int_0^2pi w = 1, and not 2pi as in BMK. """
        return self.P1P2(pt, phi) / self.anintP1P2(pt)


    ################################################
    #                                              #
    ###  Terms of squared ep->epgamma amplitude  ###
    #                                              #
    ################################################


    #### Bethe-Heitler amplitude squared Fourier coefficients

    def cBH0unpSX(self, pt):
        """ BKM Eq. (35) - small-x approximation """
        return 16. * pt.K2 * (pt.Q2/pt.t) * ( 
                F1(pt.t)**2 - (pt.t/(4.0*Mp2)) * F2(pt.t)**2 
                  ) + 8. * (2. - pt.y)**2 * ( 
                F1(pt.t)**2 - (pt.t/(4.0*Mp2)) * F2(pt.t)**2 )

    def cBH0unp(self, pt):
        """ BKM Eq. (35) """
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        FE2 = F1(t)**2 - t * F2(t)**2 / (4.0 * Mp2) 
        FM2 = (F1(t) + F2(t))**2
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
        FE2 = F1(t)**2 - t * F2(t)**2 / (4.0 * Mp2) 
        FM2 = (F1(t) + F2(t))**2
        brace = ( (4.*xB**2*Mp2/t - 2.*xB - eps2) * FE2 + 
                   2.*xB**2*(1.-(1.-2.*xB)*t/Q2) * FM2 )
        return 8. * pt.K * (2.-y) * brace

    def cBH2unp(self, pt):
        """ BKM Eq. (37) """
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        FE2 = F1(t)**2 - t * F2(t)**2 / (4.0 * Mp2) 
        FM2 = (F1(t) + F2(t))**2
        brace = 4.*Mp2/t * FE2 + 2. * FM2
        return 8. * xB**2 * pt.K2 * brace

    def TBH2unp(self, pt, phi):
        """ Bethe-Heitler amplitude squared. BKM Eq. (25)  """
        return  self.PreFacBH(pt, phi) * ( self.cBH0unp(pt) + 
                   self.cBH1unp(pt)*cos(phi) + self.cBH2unp(pt)*cos(2.*phi) )


    #### DVCS

    # {\cal C} coefficients

    def CCALDVCSunp(self, pt, pars):
        """ BKM Eq. (66) """

        xB, t = pt.xB, pt.t
        parenHH = ( RecffH(pt, pars)**2 + ImcffH(pt, pars)**2 
                +  RecffHt(pt, pars)**2 + ImcffHt(pt, pars)**2 )
        parenEH = 2.*( RecffE(pt, pars)*RecffH(pt, pars) + ImcffE(pt, pars)*ImcffH(pt, pars) 
                +  RecffEt(pt, pars)*RecffHt(pt, pars) + ImcffEt(pt, pars)*ImcffHt(pt, pars) ) 
        parenEE =  RecffE(pt, pars)**2 + ImcffE(pt, pars)**2 
        parenEtEt = RecffEt(pt, pars)**2 + ImcffEt(pt, pars)**2
        brace = 4. * (1.-xB) * parenHH - xB**2 * parenEH - (xB**2 
                + (2.-xB)**2 * t/(4.*Mp2)) * parenEE - xB**2 * t/(4.*Mp2) * parenEtEt
        return brace / (2.-xB)**2
           
    # DVCS amplitude squared Fourier coefficients

    def CDVCSunpPP(self, pt):
        """BMK Eq. (43)"""

        return 2. * (2. - 2.*pt.y + pt.y**2)


    def cDVCS0unp(self, pt, pars):
        """ BKM Eq. (43) """
        return self.CDVCSunpPP(pt) * self.CCALDVCSunp(pt, pars)

    def TDVCS2unp(self, pt, phi, pars):
        """ DVCS amplitude squared. BKM Eq. (26) - FIXME: only twist two now """
        return  self.PreFacDVCS(pt) * self.cDVCS0unp(pt, pars)

    #### Interference

    def ReCCALINTunp(self, pt, pars):
        """ Real part of BKM Eq. (69) """

        return F1(pt.t)*RecffH(pt, pars) + pt.xB/(2.-pt.xB)*(F1(pt.t)+
                F2(pt.t))*RecffHt(pt, pars) - pt.t/(4.*Mp2)*F2(pt.t)*RecffE(pt, pars)

    def ImCCALINTunp(self, pt, pars):
        """ Imag part of BKM Eq. (69) """

        return F1(pt.t)*ImcffH(pt, pars) + pt.xB/(2.-pt.xB)*(F1(pt.t)+
                F2(pt.t))*ImcffHt(pt, pars) - pt.t/(4.*Mp2)*F2(pt.t)*ImcffE(pt, pars)

    def ReDELCCALINTunp(self, pt, pars):
        """ Real part of BKM Eq. (72) """

        fx = pt.xB / (2. - pt.xB)
        return - fx * (F1(pt.t)+F2(pt.t)) * ( fx *(RecffH(pt, pars) 
            + RecffE(pt, pars)) + RecffHt(pt, pars) )

    def ImDELCCALINTunp(self, pt, pars):
        """ Imag part of BKM Eq. (72) """

        fx = pt.xB / (2. - pt.xB)
        return - fx * (F1(pt.t)+F2(pt.t)) * ( fx *(ImcffH(pt, pars) 
            + ImcffE(pt, pars)) + ImcffHt(pt, pars) )

    def ReCCALINTunpEFF(self, pt, pars):
        return 0

    def ImCCALINTunpEFF(self, pt, pars):
        return 0

    def cINT0unpSX(self, pt, pars):
        """ BKM Eq. (53) - small-x approximation!! """
        return -8. * (2. - pt.y) * (2. - 2. * pt.y + pt.y**2) *  (
                -pt.t/pt.Q2) * self.ReCCALINTunp(pt, pars)

    def cINT0unp(self, pt, pars):
        """ BKM Eq. (53) """
        return -8. * (2. - pt.y) * (
                   (2.-pt.y)**2 * pt.K2 * self.ReCCALINTunp(pt, pars) / (1.-pt.y) + 
                   (pt.t/pt.Q2) * (1.-pt.y) * (2.-pt.xB) * 
                      ( self.ReCCALINTunp(pt, pars) + self.ReDELCCALINTunp(pt, pars) ) )

    def cINT1unp(self, pt, pars):
        """ BKM Eq. (54) """
        return -8. * pt.K * (2. - 2. * pt.y + pt.y**2) * self.ReCCALINTunp(pt, pars)

    def sINT1unp(self, pt, pars):
        """ BKM Eq. (54) """
        return  8. * pt.K * pt.y * (2.-pt.y) * self.ImCCALINTunp(pt, pars)

    def cINT2unp(self, pt, pars):
        """ BKM Eq. (55) """
        return  -16. * pt.K2 * (2. - pt.y) / (2.-pt.xB) * self.ReCCALINTunpEFF(pt, pars)

    def sINT2unp(self, pt, pars):
        """ BKM Eq. (55) """
        return  16. * pt.K2 * pt.y / (2.-pt.xB) * self.ImCCALINTunpEFF(pt, pars)

    def TINTunp(self, pt, phi, lam, charge, pars):
        """ BH-DVCS interference. BKM Eq. (27) - FIXME: only twist two """
        return  - charge * self.PreFacINT(pt, phi) * ( self.cINT0unp(pt, pars)  
                + self.cINT1unp(pt, pars) * cos(phi)
                #+ self.cINT2unp(pt, pars) * cos(2.*phi) 
                + lam * self.sINT1unp(pt, pars) * sin(phi)
                #+ lam * self.sINT2unp(pt, pars) * sin(2.*phi)
                )
           
    def TINTunpd(self, pt, phi, charge, pars):
        """ BH-DVCS interference. (Normalized) part surviving after taking difference 
        of two lepton longitudinal polarization states.
        BKM Eq. (27) - FIXME: only twist two """
        return  - charge * self.PreFacINT(pt, phi) * self.sINT1unp(pt, pars) * sin(phi)

       

    def prepare(self, pt):
        """Precalculate everything that is known and attach it to data point.
        FIXME: some of the stuff here, like completing kinematics info should
        be in parent class or somewhere else, because it is not BMK-specific.

        """

        ##  --- go to BMK conventions ----
        # C1. azimutal angle phi should be in radians ...
        if pt.has('phi'):
            if pt.units['phi'][:3]== 'deg': # deg, degree, degrees -> radians
                pt.phi = pt.phi * pi / 180.
                pt.newunits['phi'] = 'rad'
        # C2. ... and in BMK convention. `frame` attribute is
        # obligatory for phi-dependent data.
            if pt.frame == 'Trento':  # Trento -> BMK
                pt.phi = pi - pt.phi
                pt.newframe = 'BMK'
        ## --- calculate standard kinematical variables ---
        # completing kinematics info
        #if pt.has('W') and pt.has('Q2'):  # it never happens
        #    pt.xB = pt.Q2 / (pt.W**2 + pt.Q2 - Mp2)
        if pt.has('xB') and pt.has('Q2'):
            pt.W = sqrt(pt.Q2 / pt.xB - pt.Q2 + Mp2)
        if pt.has('xB'):
            # There are t/Q2 corrections, cf. BMK Eq. (4), but they are 
            # formally higher twist and it is maybe sensible to DEFINE xi, 
            # the argument of CFF, as follows:
            pt.xi = pt.xB / (2. - pt.xB)
        if pt.has('t'):
            pt.mt = - pt.t
        elif pt.has('mt'):
            pt.t = - pt.mt
        ##  --- pre-calculate BMK kinematical constants ----
        # derived kinematics
        pt.y = (pt.W**2 + pt.Q2 - Mp2) / (pt.s - Mp2)
        pt.eps = 2. * pt.xB * Mp / sqrt(pt.Q2)
        pt.eps2 = pt.eps**2
        if pt.has('t'):
            pt.J = self.J(pt.Q2, pt.xB, pt.t, pt.y, pt.eps2)
            pt.K2 = self.K2(pt.Q2, pt.xB, pt.t, pt.y, pt.eps2)
            pt.K = sqrt(pt.K2)
            pt.r = self.r(pt.Q2, pt.xB, pt.t, pt.y, pt.eps2)
            # First option is numerical, second is faster
            #pt.intP1P2 = Hquadrature(lambda phi: P1P2(pt, phi), 0, 2.0*pi)
            pt.intP1P2 = self.anintP1P2(pt)

    def antiprepare(self, pt):
        """Return """

        # C1. azimutal angle phi should be in radians ...
        if pt.has('phi'):
            if pt.units['phi'][:3]== 'deg': # deg, degree, degrees -> radians
                pt.phi = pt.phi * pi / 180.
                pt.units['phi'] = 'rad'
        # C2. ... and in BKM convention. `frame` attribute is
        # obligatory for phi-dependent data.
            if pt.frame == 'Trento':  # Trento -> BKM
                pt.phi = pi - pt.phi

## Observables ##

    def Xunp(self, pt, lam, charge, pars, vars={}):
        """ 4-fold differential cross section for unpolarized target. 

        lam is lepton polarization \lambda .
        FIXME: Is this 'phi' bussiness below ugly?
        
        """
        if 'phi' in vars:
            phi = vars['phi']
        else:
            phi = pt.phi
        return self.PreFacSigma(pt) * ( self.TBH2unp(pt, phi) 
                + self.TINTunp(pt, phi, lam, charge, pars) 
                + self.TDVCS2unp(pt, phi, pars) )

    def XLU(self, pt, pars, vars={}):
        """4-fold helicity-dependent cross section measured by HALL A """

        return ( self.Xunp(pt, 1, pt.charge, pars, vars) 
                - self.Xunp(pt, -1, pt.charge, pars, vars) ) / 2.

    def XUU(self, pt, pars, vars={}):
        """4-fold helicity-independent cross section measured by HALL A """

        return ( self.Xunp(pt, 1, pt.charge, pars, vars) 
                + self.Xunp(pt, -1, pt.charge, pars, vars) ) / 2.

    def BCA(self, pt, phi, pars):
        """Beam charge asymmetry. """

        if not self.optimization:
            # use defining formula:  (sigma+ - sigma-)/(sigma+ + sigma-)
            # i.e. charge is not read from datafile!
            return (
               self.Xunp(pt, 0, 1, pars, {'phi':phi}) - self.Xunp(pt, 0, -1, pars, {'phi':phi}) )/(
               self.Xunp(pt, 0, 1, pars, {'phi':phi}) + self.Xunp(pt, 0, -1, pars, {'phi':phi}) )
        else:
            # optimized formula (remove parts which cancel anyway)
            return  self.TINTunp(pt, phi, 0, 1, pars) / ( 
                           self.TBH2unp(pt, phi) + self.TDVCS2unp(pt, phi, pars) )

    def BCSD(self, pt, lam, pars, vars={}):
        """4-fold beam charge-spin cross section difference measured by COMPASS """

        # charge is not read from datafile!
        return (self.Xunp(pt, lam, 1, pars, vars) - self.Xunp(pt, -lam, -1, pars, vars))/2.

    def BCSS(self, pt, lam, pars, vars={}):
        """4-fold beam charge-spin cross section sum measured by COMPASS. """

        # charge is not read from datafile!
        return (self.Xunp(pt, lam, 1, pars, vars) + self.Xunp(pt, -lam, -1, pars, vars))/2.

    def BCSA(self, pt, lam, pars, vars={}):
        """Beam charge-spin asymmetry as measured by COMPASS. """

        if not self.optimization:
            return  self.BCSD(pt, lam, pars, vars) / self.BCSS(pt, lam, pars, vars)
        else:
            if 'phi' in vars:
                phi = vars['phi']
            else:
                phi = pt.phi
            return  self.TINTunp(pt, phi, 0, 1, pars) / ( 
                  self.TBH2unp(pt, phi) + self.TDVCS2unp(pt, phi, pars) 
               + (self.TINTunp(pt, phi, lam, 1, pars) + self.TINTunp(pt, phi, -lam, -1, pars))/2.)

    def BSA(self, pt, pars, vars={}):
        """Beam spin asymmetry. HERMES uses positron convention (charge=+1)
           CLAS uses electron (charge=-1), so datafile has to provide this
           information in the name of 'in1' particle."""

        if 'phi' in vars:
            phi = vars['phi']
        else:
            phi = pt.phi
        if not self.optimization:
            # use defining formula:  
            return (
               self.Xunp(pt, 1, pt.charge, pars, {'phi':phi})  
                     - self.Xunp(pt, -1, pt.charge, pars, {'phi':phi}) )/(
               self.Xunp(pt, 1, pt.charge, pars, {'phi':phi}) 
                     + self.Xunp(pt, -1, pt.charge, pars, {'phi':phi}) )
        else:
            # optimized formula (by removing parts which cancel anyway)
            return  self.TINTunpd(pt, phi, pt.charge, pars) / ( self.TBH2unp(pt, phi) 
                + self.TDVCS2unp(pt, phi, pars) + self.TINTunp(pt, phi, 0, pt.charge, pars) )

    def PartialCrossSection4int(self, t, pt, pars):
        """Same as PartialCrossSection but with additional variable t 
        to facilitate integration over it.
        
        """
        pt.t = t
        return PartialCrossSection(pt, pars)

    def PartialCrossSection(self, pt, pars):
        """Partial DVCS cross section w.r.t. Mandelstam t."""

        W2 = pt.W * pt.W
        return 260.5633976788416 * W2 * (ImcffH(pt, pars)**2 
                + RecffH(pt, pars)**2)  / (
            (W2 + pt.Q2) * (2.0 * W2 + pt.Q2)**2 )

    def TotalCrossSection(self, pt, pars):
        """Total DVCS cross section."""

        #res = Hquadrature(lambda t: self.PartialCrossSection4int(t, pt, pars), 0, 1)
        res = PartialCrossSection4int(-0.22, pt, pars)
        return res

    def BCAcos0(self, pt, pars):
        res = Hquadrature(lambda phi: self.BCA(pt, phi, pars), 0, 2.0*pi)
        return res / (2.0*pi)

    def BCAcos1(self, pt, pars):
        res = Hquadrature(lambda phi: self.BCA(pt, phi, pars) * cos(phi), 0, 2*pi)
        return  - res / pi

    def ALUIsin1orig(self, pt, pars):
        res = Hquadrature(lambda phi: self.BSA(pt, pars, {'phi':phi}) * sin(phi), 0, 2*pi)
        return  res / pi

    def ALUIsin1(self, pt, pars):
        """FIXME: shortcut, charge conventions, ... for HERMES - negative. FIXED?"""
        return  self.BSA(pt, pars, {'phi':pi/2.}) 

    def ALUa(self, pt, pars):
        """FIXME: charge conventions ... for CLAS - positive. FIXED?"""
        return  self.BSA(pt, pars, {'phi':pi/2.}) 

    def BCA0minusr1(self, pt, pars):
        return self.BCAcos0(pt, pars) - pt.r * self.BCAcos1(pt, pars)

    def ImCI(self, pt, pars):
        """ As defined by HALL A """
        return self.ImCCALINTunp(pt, pars)

    def ReCI(self, pt, pars):
        """ As defined by HALL A """
        return self.ReCCALINTunp(pt, pars)

    def ReCpDCI(self, pt, pars):
        """ As defined by HALL A """
        return self.ReCCALINTunp(pt, pars) + self.ReDELCCALINTunp(pt, pars)

    def b1ovb0(self, pt, pars):
        """Ratio of first two cos harmonics of w-weighted cross section. In BMK, not Trento??"""

        b0 = Hquadrature(lambda phi: self.w(pt, phi) * self.XUU(pt, pars, {'phi':phi}), 
                0, 2.0*pi) / (2.0*pi)
        b1 = Hquadrature(lambda phi: self.w(pt, phi) * self.XUU(pt, pars, {'phi':phi}) * cos(phi), 
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
                (1. + eps2)**2*((8.*pt.K*(2. - y)*y)/(1. + eps2))*
         (1. + ((t - self.tmin(Q2, xB, eps2))*(1. - xB + 0.5*(-1. + sqrt(1. + eps2))))/
           (Q2*(1. + eps2)))
         )


    def cINT0unp(self, pt, pars):
        """ hotfixed BKM Eq. (53) """
        return (self.CINTunpPP0(pt) * self.ReCCALINTunp(pt, pars)
                + self.CINTunpPP0q(pt) * self.ReDELCCALINTunp(pt, pars))

    def cINT1unp(self, pt, pars):
        """ hotfixed BKM Eq. (54) """
        return self.CINTunpPP1(pt) * self.ReCCALINTunp(pt, pars)

    def sINT1unp(self, pt, pars):
        """ hotfixed BKM Eq. (54) """
        return self.SINTunpPP1(pt) * self.ImCCALINTunp(pt, pars)

