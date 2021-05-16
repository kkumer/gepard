"""Definition of theory frameworks."""

from joblib import Parallel, delayed
from numpy import array, cos, linspace, ndarray, pi, sin, sqrt, transpose
from scipy.stats import scoreatpercentile

import gepard.data
import gepard.model
import gepard.quadrature
from gepard.constants import GeV2nb, Mp, Mp2, alpha

NCPU = 23  # how many CPUs to use in parallel


class Theory(object):
    """Class of theory frameworks for calculation of observables.

    This is a base class and subclasses implement various observables as
    methods. Instances are instantiated with specific model (instance of
    `gepard.model.Model`) for structure functions such as CFFs and elastic FFs,
    so theory object finally contains complete algorithm for calculation of
    specific observable.

    Implemented subclasses are (FIXME: not yet transfered):

     - BMK - From Belitsky, Mueller and Kirchner arXiv:hep-ph/0112108
     - hotfixedBMK  - simple valence xB region improvement by Dieter
     - BM10ex
     - BM10 - Belitsky, Mueller and Ji
     - BM10tw2 - BM10 with higher twists set to zero
    """

    def __init__(self, model: gepard.model.Model) -> None:
        """Initialize with m as instance of `g.Model` class."""
        self.model = model
        self.m = self.model  # shortcut
        self.name = model.name
        self.texname = model.texname
        self.description = model.description

    def chisq_single(self, points: gepard.data.DataSet, asym: bool = False,
                     **kwargs) -> float:
        """Return total chi-square.

        Args:
            points: measurements with uncertainties, observable is
                named in `yaxis` attribute of each DataPoint,
                value is `val` and uncertainty is `err`.
            asym: if measurements provide asymmetric uncertainties
                `errplus` and `errminus`, this enables their usage
            observable (str): overrides `yaxis` of DataPoints

        Notes:
            If the theory or model provide uncertainties, they are ignored -
            only experimental uncertainties are taken into account.
           """
        allpulls = []
        for pt in points:
            diff = (self.predict(pt, observable=pt.yaxis, **kwargs) 
                          - pt.val)
            if asym:
                if diff > 0:
                    allpulls.append(diff/pt.errplus)
                else:
                    allpulls.append(diff/pt.errminus)
            else:
                allpulls.append(diff/pt.err)
        chi = sum(p*p for p in allpulls)  # equal to m.fval if minuit fit is done
        return chi

    def pull(self, pt: gepard.data.DataPoint):
        """Return pull of a single Datapoint."""
        return (self.predict(pt, observable=pt.yaxis) - pt.val) / pt.err

    def chisq_para(self, points: gepard.data.DataSet, asym: bool = False,
                   **kwargs) -> float:
        """Return total chi-square - parallel version.

        Warning:
            Cannot be used until underlying global Fortran variables can change
            during session. (Like different kinematics of data points.)
        """
        allpulls = Parallel(n_jobs=NCPU)(delayed(self.pull)(pt) for pt in points)
        chi = sum(p*p for p in allpulls)  # equal to m.fval if minuit fit is done
        return chi

    chisq = chisq_para

    def predict(self, pt, error=False, CL=False, **kwargs):
        """Give prediction for DataPoint pt.

        Args:
            pt: instance of DataPoint
            error: if available, produce tuple (mean, error)
            CL: (NNet only) error is not std.dev.,
                but 68% C.L. (mean, errplus, errminus)
            observable: string. Default is pt.yaxis. It is acceptable also
                        to pass CFF as observable, e.g., observable = 'ImH'
            parameters: dictionary which will temporarily update model's one
            orig_conventions: give prediction using original conventions of
                              the given DataPoint (e.g. for plotting)
        """
        m = self.model
        if 'observable' in kwargs:
            obs = kwargs['observable']
        else:
            obs = pt.yaxis

        if 'parameters' in kwargs:
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

        # if obs in m.allCFFs or obs in m.allGPDs:
            # # we want GPD/CFF, which are model attributes
            # fun = getattr(m, obs)
        # else:
            # # we want a "real" observable, which is theory attribute
            # fun = getattr(self, obs)

        fun = getattr(self, obs)

        if error:
            try:
                # We now do standard propagation of error from m to observable
                pars = [p for p in m.parameters if not m.parameters_fix['p']]
                var = 0
                dfdp = {}
                for p in pars:
                    # calculating dfdp = derivative of observable w.r.t. parameter:
                    h=sqrt(m.covariance[p,p])
                    mem = m.parameters[p]
                    m.parameters[p] = mem+h/2.
                    if 'g' in self.model.__dict__: self.m.g.newcall = 1
                    up = fun(pt)
                    m.parameters[p] = mem-h/2.
                    if 'g' in self.model.__dict__: self.m.g.newcall = 1
                    down = fun(pt)
                    m.parameters[p] = mem
                    dfdp[p] = (up-down)/h
                for p1 in pars:
                    for p2 in pars:
                        var += dfdp[p1]*m.covariance[p1,p2]*dfdp[p2]
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
            if 'g' in self.model.__dict__: self.m.g.newcall = 1
            result = fun(pt)
            if isinstance(result, ndarray):
                # we have neural net
                result = result.mean()

        if 'parameters' in kwargs:
            # restore old values
            self.model.parameters.update(old)

        if kwargs.pop('orig_conventions', False):
            # express result in conventions of original datapoint
            try:
                result = (self.orig_conventions(pt, result[0]),) + result[1:]
            except (IndexError, TypeError):
                result = self.orig_conventions(pt, result)
        return result


class BMK(Theory):
    """Implementation of formulas from hep-ph/0112108  (BMK)"""

    ### Kinematics ###
    # (implemented as static and class methods)

    def tmin(Q2, xB, eps2):
        """BMK Eq. (31)"""
        return -Q2 * ( 2. * (1.-xB)*(1. - sqrt(1.+eps2)) + eps2 ) / (
                4. * xB * (1.-xB) + eps2 )
    tmin = staticmethod(tmin)

    def tmax(Q2, xB, eps2):
        return -Q2 * ( 2. * (1.-xB)*(1. + sqrt(1.+eps2)) + eps2 ) / (
                4. * xB * (1.-xB) + eps2 )
    tmax = staticmethod(tmax)

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
        """ Analitical integral of  P1 P2  """
        xB, Q2, t, y, eps2, K2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2, pt.K2
        brace = ( (1 - y - (1+eps2/2.) * y**2 * eps2/2.) * (1. + t/Q2)**2 +
                  2.*K2 - (1.-xB)*(2.-y)**2 * (1. + xB*t/Q2) * t/Q2 )
        return -2. * pi * brace / (1+eps2)**2 / y**2
    anintP1P2 = staticmethod(anintP1P2)

    def PreFacSigma(self, pt):
        """ Prefactor of 4-fold xs. Take prefactor in Eq(22) times 2pi because of proton Phi integration
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
        It is normalized to int_0^2pi w  2pi as in BMK. """
        return 2.*pi*pt.P1P2 / pt.intP1P2

    def prepare(pt):
        """Pre-calculate GPD-independent kinamatical constants and functions."""
        if not hasattr(pt, "s"):
            #This is for variable beam energy;  code duplication
            if pt.process in ['ep2epgamma', 'en2engamma']:
                if pt.exptype == 'fixed target':
                    pt.s = 2 * Mp * pt.in1energy + Mp2
                elif pt.exptype == 'collider':
                    pt.s = 2 * pt.in1energy * (pt.in2energy + math.sqrt(
                        pt.in2energy**2 - Mp2)) + Mp2
                else:
                    pass # FIXME: raise error
            else:
                pass # FIXME: should I raise error here?
        pt.y = (pt.W**2 + pt.Q2 - Mp2) / (pt.s - Mp2)
        pt.eps = 2. * pt.xB * Mp / sqrt(pt.Q2)
        pt.eps2 = pt.eps**2
        if 't' in pt:
            pt.J = BMK.J(pt.Q2, pt.xB, pt.t, pt.y, pt.eps2)
            pt.K2 = BMK.K2(pt.Q2, pt.xB, pt.t, pt.y, pt.eps2)
            pt.K = sqrt(pt.K2)
            pt.tK2 = pt.K2*pt.Q2/(1-pt.y-pt.eps2*pt.y**2/4.)
            pt.tK = sqrt(pt.tK2)
            pt.r = BMK.r(pt.Q2, pt.xB, pt.t, pt.y, pt.eps2)
            # Needed for BMP higher-twist stuff:
            pt.chi0 = sqrt(2.*pt.Q2)*pt.tK/sqrt(1+pt.eps2)/(pt.Q2+pt.t)
            pt.chi = (pt.Q2-pt.t+2.*pt.xB*pt.t)/sqrt(1+pt.eps2)/(pt.Q2+pt.t) - 1.
            # First option is numerical, second is analytical and faster
            #pt.intP1P2 = g.quadrature.Hquadrature(lambda phi: P1P2(pt, phi), 0, 2.0*pi)
            pt.intP1P2 = BMK.anintP1P2(pt)
        if 'phi' in pt:
            pt.P1P2 = BMK.P1P2(pt)
    prepare = staticmethod(prepare)

    ### Some kinematical functions which are not in BKM paper but are convenient
    ### to define here

    def long2trans(self, pt):
        """ Ratio of longitudinal to transverse photon flux 1304.0077 Eq. (2.9) """
        return (1.-pt.y-pt.eps2*pt.y**2/4.)/(1-pt.y+pt.y**2/2+pt.eps2*pt.y**2/4.)

    def HandFlux(self, pt):
        """ Virtual photon flux (Hand convention) 1304.0077 Eq. (2.9) """
        return (alpha/2./pi)*(pt.y**2/(1.-self.long2trans(pt)))*(1-pt.xB)/pt.xB/pt.Q2

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
                self.m.F1(pt)**2 - (pt.t/(4.0*Mp2)) * self.m.F2(pt)**2
                  ) + 8. * (2. - pt.y)**2 * (
                self.m.F1(pt)**2 - (pt.t/(4.0*Mp2)) * self.m.F2(pt)**2 )

    def cBH0unp(self, pt):
        """ BKM Eq. (35) """
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        FE2 = self.m.F1(pt)**2 - t * self.m.F2(pt)**2 / (4.0 * Mp2)
        FM2 = (self.m.F1(pt) + self.m.F2(pt))**2
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
        FE2 = self.m.F1(pt)**2 - t * self.m.F2(pt)**2 / (4.0 * Mp2)
        FM2 = (self.m.F1(pt) + self.m.F2(pt))**2
        brace = ( (4.*xB**2*Mp2/t - 2.*xB - eps2) * FE2 +
                   2.*xB**2*(1.-(1.-2.*xB)*t/Q2) * FM2 )
        return 8. * pt.K * (2.-y) * brace

    def cBH2unp(self, pt):
        """ BKM Eq. (37) """
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        FE2 = self.m.F1(pt)**2 - t * self.m.F2(pt)**2 / (4.0 * Mp2)
        FM2 = (self.m.F1(pt) + self.m.F2(pt))**2
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
        F1, F2 = self.m.F1(pt), self.m.F2(pt)
        sqrt1yeps = sqrt(1-y-eps2*y**2/4.)
        brace = ( xB**3*Mp2/Q2*(1-t/Q2)*(F1+F2) +
                (1-(1-xB)*t/Q2)*(xB**2*Mp2/t*(1-t/Q2)*F1+xB/2.*F2) )
        return (-8*pt.in1polarization*cos(pt.varphi)*(2-y)*y*sqrt(Q2)/Mp*
                sqrt(1+eps2)*pt.K/sqrt1yeps*(F1+F2)*brace)

    def cBH1TP(self, pt):
        """ BKM Eq. (41) """
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        F1, F2 = self.m.F1(pt), self.m.F2(pt)
        sqrt1yeps = sqrt(1-y-eps2*y**2/4.)
        brace = ( 2*pt.K**2*Q2/t/sqrt1yeps**2 * (xB*(1-t/Q2)*F1 +
            t/4./Mp2*F2) + (1+eps2)*xB*(1-t/Q2)*(F1+t/4./Mp2*F2) )
        return (-16*pt.in1polarization*cos(pt.varphi)*xB*y*
                sqrt1yeps*Mp/sqrt(Q2)*sqrt(1+eps2)*(F1+F2)*brace)

    def sBH1TP(self, pt):
        """ BKM Eq. (42) """
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        F1, F2 = self.m.F1(pt), self.m.F2(pt)
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

        return self.m.F1(pt)*self.m.ReH(pt) + pt.xB/(2.-pt.xB)*(self.m.F1(pt)+
                self.m.F2(pt))*self.m.ReHt(pt) - pt.t/(4.*Mp2)*self.m.F2(pt)*self.m.ReE(pt)

    def ImCCALINTunp(self, pt):
        """ Imag part of BKM Eq. (69) """

        return self.m.F1(pt)*self.m.ImH(pt) + pt.xB/(2.-pt.xB)*(self.m.F1(pt)+
                self.m.F2(pt))*self.m.ImHt(pt) - pt.t/(4.*Mp2)*self.m.F2(pt)*self.m.ImE(pt)

    def ReDELCCALINTunp(self, pt):
        """ Real part of BKM Eq. (72) """

        fx = pt.xB / (2. - pt.xB)
        return - fx * (self.m.F1(pt)+self.m.F2(pt)) * ( fx *(self.m.ReH(pt)
            + self.m.ReE(pt)) + self.m.ReHt(pt) )

    def ImDELCCALINTunp(self, pt):
        """ Imag part of BKM Eq. (72) """

        fx = pt.xB / (2. - pt.xB)
        return - fx * (self.m.F1(pt)+self.m.F2(pt)) * ( fx *(self.m.ImH(pt)
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

        xB, t, F1, F2 = pt.xB, pt.t, self.m.F1(pt), self.m.F2(pt)
        H, E, Ht, Et = self.m.ReH(pt), self.m.ReE(pt), self.m.ReHt(pt), self.m.ReEt(pt)
        brace1 = xB**2/(2-xB)*(H+xB/2.*E) + xB*t/4./Mp2*E
        brace2 = 4*(1-xB)/(2-xB)*F2*Ht - (xB*F1+xB**2/(2-xB)*F2)*Et
        return (F1+F2)*brace1 - xB**2/(2-xB)*F1*(Ht+xB/2.*Et) + t/4./Mp2*brace2

    def ImCCALINTTPp(self, pt):
        """ Imag part of BKM Eq. (71) FIXME: code duplication """

        xB, t, F1, F2 = pt.xB, pt.t, self.m.F1(pt), self.m.F2(pt)
        H, E, Ht, Et = self.m.ImH(pt), self.m.ImE(pt), self.m.ImHt(pt), self.m.ImEt(pt)
        brace1 = xB**2/(2-xB)*(H+xB/2.*E) + xB*t/4./Mp2*E
        brace2 = 4*(1-xB)/(2-xB)*F2*Ht - (xB*F1+xB**2/(2-xB)*F2)*Et
        return (F1+F2)*brace1 - xB**2/(2-xB)*F1*(Ht+xB/2.*Et) + t/4./Mp2*brace2

    def ReCCALINTTPm(self, pt):
        """ Real part of BKM Eq. (71) """

        xB, t, F1, F2 = pt.xB, pt.t, self.m.F1(pt), self.m.F2(pt)
        H, E, Ht, Et = self.m.ReH(pt), self.m.ReE(pt), self.m.ReHt(pt), self.m.ReEt(pt)
        xBaux = xB**2/(2-xB)
        paren1 = xB**2*F1 - (1-xB)*t/Mp2*F2
        brace = t/4./Mp2*((2-xB)*F1 + xBaux*F2) + xBaux*F1
        return paren1*H/(2-xB) + brace*E - xBaux*(F1+F2)*(Ht+t/4./Mp2*Et)

    def ImCCALINTTPm(self, pt):
        """ Imag part of BKM Eq. (71)  FIXME: code duplication"""

        xB, t, F1, F2 = pt.xB, pt.t, self.m.F1(pt), self.m.F2(pt)
        H, E, Ht, Et = self.m.ImH(pt), self.m.ImE(pt), self.m.ImHt(pt), self.m.ImEt(pt)
        xBaux = xB**2/(2-xB)
        paren1 = xB**2*F1 - (1-xB)*t/Mp2*F2
        brace = t/4./Mp2*((2-xB)*F1 + xBaux*F2) + xBaux*F1
        return paren1*H/(2-xB) + brace*E - xBaux*(F1+F2)*(Ht+t/4./Mp2*Et)

    def ReDELCCALINTTPp(self, pt):
        """ Real part of BKM Eq. (74) """

        xB, t, F1, F2 = pt.xB, pt.t, self.m.F1(pt), self.m.F2(pt)
        Ht, Et = self.m.ReHt(pt), self.m.ReEt(pt)
        return -t/Mp2*(F2*Ht - xB/(2-xB)*(F1+xB*F2/2)*Et)

    def ImDELCCALINTTPp(self, pt):
        """ Imag part of BKM Eq. (74) FIXME: code duplication """

        xB, t, F1, F2 = pt.xB, pt.t, self.m.F1(pt), self.m.F2(pt)
        Ht, Et = self.m.ImHt(pt), self.m.ImEt(pt)
        return -t/Mp2*(F2*Ht - xB/(2-xB)*(F1+xB*F2/2)*Et)

    def ReDELCCALINTTPm(self, pt):
        """ Real part of BKM Eq. (75) """

        xB, t, F1, F2 = pt.xB, pt.t, self.m.F1(pt), self.m.F2(pt)
        H, E = self.m.ReH(pt), self.m.ReE(pt)
        return t/Mp2*(F2*H - F1*E)

    def ImDELCCALINTTPm(self, pt):
        """ Imag part of BKM Eq. (75) FIXME: code duplication """

        xB, t, F1, F2 = pt.xB, pt.t, self.m.F1(pt), self.m.F2(pt)
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
        if 'vars' in kwargs:
            ptvars = gepard.data.DummyPoint(init=kwargs['vars'])
            kin = gepard.data._fill_kinematics(ptvars, old=pt)
            BMK.prepare(kin)
        else:
            # just copy everything from pt
            ptempty = gepard.data.DummyPoint()
            kin = gepard.data._fill_kinematics(ptempty, old=pt)
            BMK.prepare(kin)
            ## Nothing seems to be gained by the following approach:
            #kin = dict((i, getattr(pt, i)) for i in
            #        ['xB', 'Q2', 'W', 's', 't', 'mt', 'phi', 'in1charge',
            #            'in1polarization', 'in2particle'])

        # Copy non-kinematical info
        for atr in ['in1charge', 'in1polarization', 'in2polarization', 'in2particle']:
            if atr in pt:
                setattr(kin, atr, getattr(pt, atr))

        # For efficient calculation of XS with unpolarized beam
        if 'zeropolarized' in kwargs and kwargs['zeropolarized']:
            kin.in1polarization = 0

        # Flipping spins and/or charges for asymmetries
        if 'flip' in kwargs and kwargs['flip']:
            if isinstance(kwargs['flip'], list):
                for item in kwargs['flip']:
                    setattr(kin, item, - getattr(pt, item))
            else:
                setattr(kin, kwargs['flip'], - getattr(pt, kwargs['flip']))

        # Weighting the integrand by BH propagators
        if 'weighted' in kwargs and kwargs['weighted']:
            wgh = self.w(kin)
        else:
            wgh = 1

        # Gepard may need resetting
        if 'g' in self.model.__dict__:
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

        # Gepard may need resetting
        if 'g' in self.model.__dict__:
            self.m.g.parint.pid = 0
            self.m.g.newcall = 1
        res = self.m.DISF2(pt)
        return res

    def _XDVCStApprox(self, pt):
        """Partial DVCS cross section w.r.t. Mandelstam t.
         Approx. formula used in NPB10 paper."""

        W2 = pt.W * pt.W
        # Simplified formula used also in Fortran gepard code
        ReH, ImH, ReE, ImE, ReHt, ImHt, ReEt, ImEt = self.m.cff(pt.xi, pt.t, pt.Q2)
        res = 260.5633976788416 * W2 * (
                (ImH**2 + ReH**2)
                - pt.t/(4.*Mp2)*(ReE**2 + ImE**2)) / (
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
        ReH, ImH, ReE, ImE, ReHt, ImHt, ReEt, ImEt = self.m.cff(pt.xi, pt.t, pt.Q2)
        # if pt.t < -0.9:
        #     print('t, Q2:  ReH, ImH = {:.2} {}: {:.3}, {:5.1f} '.format(pt.t, pt.Q2, ReH, ImH))
        res = 65.14079453579676 * (pt.xB**2 / pt.Q2**2 / (1-pt.xB) / (2-pt.xB)**2 /
                sqrt(1 + eps2) * (
                    4 * (1 - pt.xB) * (ImH**2 + ReH**2)
                - pt.xB**2 * (ReE**2+ImE**2 + 2*ReE*ReH + 2*ImE*ImH)
                - (2-pt.xB)**2 * pt.t/4/Mp2*(ImE**2+ReE**2)))
        return res

    # _XDVCSt = _XDVCStApprox
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
            aux.append(res)
        return array(aux)


    def X(self, pt):
        """Total DVCS or DVMP cross section. """
        if 't' in pt or 'tm' in pt:
            # partial XS w.r.t momentum transfer t
            if hasattr(pt, 'process') and pt.process == 'gammastarp2rho0p':
                return self._Xrhot(pt)
            else:
                return self._XDVCSt(pt)

        else:
            # total XS
            if 'tmmax' in pt:
                tmmax = pt.tmmax
            else:
                tmmax = 1.  # default -t cuttoff in GeV^2
            res = gepard.quadrature.tquadrature(lambda t: self._Xt4int(t, pt), -tmmax, 0)
            return res


## General assymetries
## TODO: Lot of code duplication here - this should be united in one clever function

    def _phiharmonic(self, fun, pt, **kwargs):
        """Return fun evaluated for phi=pt.phi, or harmonic of fun
        corresponding to pt.FTn.

        """
        if 'phi' in pt or ('vars' in kwargs
                and 'phi' in kwargs['vars']):
            return fun(pt, **kwargs)
        elif 'FTn' in pt:
            if pt.FTn < 0:
                res = gepard.quadrature.Hquadrature(lambda phi:
                        fun(pt, vars={'phi':phi}, **kwargs) * sin(-pt.FTn*phi), 0, 2*pi)
            elif pt.FTn > 0:
                res = gepard.quadrature.Hquadrature(lambda phi:
                        fun(pt, vars={'phi':phi}, **kwargs) * cos(pt.FTn*phi), 0, 2*pi)
            elif pt.FTn == 0:
                res = gepard.quadrature.Hquadrature(lambda phi:
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

    def _BSD(self, pt, **kwargs):
        """Calculate 4-fold helicity-dependent cross section"""

        R = kwargs.copy()
        R.update({'flip':'in1polarization'})
        return ( self.Xunp(pt, **kwargs)
                - self.Xunp(pt, **R) ) / 2.

    def BSD(self, pt, **kwargs):
        """Calculate beam spin difference (BSD) or its harmonics."""
        return self._phiharmonic(self._BSD, pt, **kwargs)

    def _BSS(self, pt, **kwargs):
        """4-fold helicity-independent cross section"""
        R = kwargs.copy()
        R.update({'flip':'in1polarization'})
        return ( self.Xunp(pt, **kwargs)
                + self.Xunp(pt, **R) ) / 2.

    def XSintphi(self, pt, **kwargs):
        """Return XS integrated over azimuthal angle"""
        mem = pt.__dict__.pop('FTn', None)
        pt.FTn = 0
        res = self.BSS(pt, **kwargs)
        # restore old value if needed
        if mem: pt.in2polarization = mem
        return 2*pi*self.BSS(pt, **kwargs)

    def BSS(self, pt, **kwargs):
        """Calculate beam spin sum (BSS) or its harmonics."""
        return self._phiharmonic(self._BSS, pt, **kwargs)

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
        if 'phi' in pt:
            return self._BSA(pt, **kwargs)
        elif 'FTn' in pt and pt.FTn == -1:
            # FIXME: faster shortcut (approximate!)
            if 'vars' in kwargs:
                kwargs['vars'].update({'phi':pi/2.})
            else:
                kwargs['vars'] = {'phi':pi/2.}
            return  self._BSA(pt, **kwargs)
        else:
            raise ValueError('[%s] has neither azimuthal angle phi\
 nor harmonic FTn = -1 defined!' % pt)
        ### Exact but slower:
            #res = g.quadrature.Hquadrature(lambda phi:
            #        self._BSA(pt, vars={'phi':phi}) * sin(phi), 0, 2*pi)
            #return  res / pi

    def BSAexact(self, pt, **kwargs):
        """Calculate beam spin asymmetry (BSA) or its harmonics."""
        if 'phi' in pt:
            return self._BSA(pt, **kwargs)
        elif 'FTn' in pt and pt.FTn == -1:
            res = gepard.quadrature.Hquadrature(lambda phi:
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
        b0 = gepard.quadrature.Hquadrature(lambda phi: self.BSS(pt, vars={'phi':phi}, weighted=True),
                0, 2.0*pi) / (2.0*pi)
        b1 = gepard.quadrature.Hquadrature(lambda phi: self.BSS(pt, vars={'phi':phi}, weighted=True) * cos(phi),
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
