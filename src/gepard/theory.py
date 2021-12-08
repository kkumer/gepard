"""Definition of theory frameworks."""

# from joblib import Parallel, delayed
from numpy import array, cos, linspace, ndarray, pi, sin, sqrt, transpose
from scipy.stats import scoreatpercentile

from . import data, model, quadrature
from .constants import GeV2nb, Mp, Mp2, alpha
from .kinematics import K2, tmin, weight_BH

NCPU = 23  # how many CPUs to use in parallel


class Theory(object):
    """Class of theory frameworks for calculation of observables.

    This is a base class and subclasses implement various observables as
    methods. Instances are instantiated with specific model (instance of
    `model.Model`) for structure functions such as CFFs and elastic FFs,
    so theory object finally contains complete algorithm for calculation of
    specific observable.
    """

    def __init__(self, model: model.Model) -> None:
        """Construct theory framework with specific model."""
        self.model = model
        self.m = self.model  # shortcut
        self.name = model.name
        # self.texname = model.texname
        # self.description = model.description

    def chisq_single(self, points: data.DataSet, asym: bool = False,
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
            diff = (self.predict(pt, observable=pt.yaxis, **kwargs) - pt.val)
            if asym:
                if diff > 0:
                    allpulls.append(diff/pt.errplus)
                else:
                    allpulls.append(diff/pt.errminus)
            else:
                allpulls.append(diff/pt.err)
        chi = sum(p*p for p in allpulls)  # equal to m.fval if minuit fit is done
        return chi

    def pull(self, pt: data.DataPoint):
        """Return pull of a single Datapoint."""
        return (self.predict(pt, observable=pt.yaxis) - pt.val) / pt.err

#     def chisq_para(self, points: DataSet, asym: bool = False,
#                    **kwargs) -> float:
#         """Return total chi-square - parallel version.
# 
#         Warning:
#             Cannot be used until underlying global Fortran variables can change
#             during session. (Like different kinematics of data points.)
#         """
#         allpulls = Parallel(n_jobs=NCPU)(delayed(self.pull)(pt) for pt in points)
#         chi = sum(p*p for p in allpulls)  # equal to m.fval if minuit fit is done
#         return chi

    chisq = chisq_single

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
                    up = fun(pt)
                    m.parameters[p] = mem-h/2.
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
                result = (pt.orig_conventions(result[0]),) + result[1:]
            except (IndexError, TypeError):
                result = pt.orig_conventions(result)
        return result


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
            # print("t = {}  =>  dsig/dt = {}".format(t_single, res))
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
            res = quadrature.tquadrature(lambda t: self._Xt4int(t, pt), -tmmax, 0)
            return res


class DVCS(Theory):
    """DVCS observables base class.

    Implements cross-section for electroproduction of real photon, and derived
    cross-sections, asymmetries, chi-squares etc. 
    Subclasses should, as their methods, implement expressions for squared
    BH, interference and squared DVCS amplitudes:
    - TBH2unp, TINTunp, TDVCS2unp  (for unpolarized target)
    - TBH2LP, TINTLP, TDVCS2LP  (for longitudinally polarized target)
    - TBH2TP, TINTTP, TDVCS2TP  (for transversally polarized target)

    Implemented subclasses are:

     - BMK - From Belitsky, Mueller and Kirchner arXiv:hep-ph/0112108
     - hotfixedBMK  - simple valence xB region improvement by Dieter
     - BM10ex
     - BM10 - Belitsky, Mueller and Ji
     - BM10tw2 - BM10 with higher twists set to zero
    """

    def PreFacSigma(self, pt):
        """ Prefactor of 4-fold XS. 

        Take prefactor in Eq(22) times 2pi because of proton Phi integration
        and times y/Q2 because of dy -> dQ2. Convert to nanobarns.

        """
        return alpha**3 * pt.xB * pt.y**2 / (8. * pi * pt.Q2**2 * 
                                             sqrt(1.+pt.eps2)) * GeV2nb


    def XS(self, pt, **kwargs):
        """Differential 4-fold e p --> e p gamma cross section.

        """
        # Overriding pt kinematics with those from kwargs
        if 'vars' in kwargs:
            ptvars = data.DataPoint(init=kwargs['vars'])
            data._fill_kinematics(ptvars, old=pt)
            kin = ptvars
        else:
            # just copy everything from pt
            ptempty = data.DataPoint()
            data._fill_kinematics(ptempty, old=pt)
            kin = ptempty
            # Nothing seems to be gained by the following approach:
            # kin = dict((i, getattr(pt, i)) for i in
            #        ['xB', 'Q2', 'W', 's', 't', 'mt', 'phi', 'in1charge',
            #            'in1polarization', 'in2particle'])

        # Copy non-kinematical info
        for atr in ['in1charge', 'in1polarization', 'in2polarization', 'in2particle',
                    'process', 'exptype', 'in1energy', 'in2energy', 'yaxis']:
            if atr in pt:
                setattr(kin, atr, getattr(pt, atr))

        kin.prepare()

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
            wgh = weight_BH(kin)
        else:
            wgh = 1

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

        This is a convenience function that calculates XS as if the target were
        unpolarized even if it is i.e. if pt.in2polarization has some value.

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
        pol.update({'flip': 'in2polarization'})
        o =  self.XS(pt, **kwargs)
        f =  self.XS(pt, **pol)
        return (o-f)/2.

    def _XDVCStApprox(self, pt):
        """Partial DVCS cross section w.r.t. Mandelstam t.

        Approx. formula used in NPB10 paper.
        """
        W2 = pt.W * pt.W
        # Simplified formula used also in Fortran gepard code
        ReH, ImH, ReE, ImE, ReHt, ImHt, ReEt, ImEt = self.m.cff(pt)
        res = 260.5633976788416 * W2 * (
                (ImH**2 + ReH**2)
                - pt.t/(4.*Mp2)*(ReE**2 + ImE**2)) / (
            (W2 + pt.Q2) * (2.0 * W2 + pt.Q2)**2 )
        return res


    def _XDVCStEx(self, pt):
        """Partial DVCS cross section w.r.t. Mandelstam t."""

        eps2 = 4. * pt.xB**2 * Mp2 / pt.Q2
        try:
            ReH, ImH, ReE, ImE, ReHt, ImHt, ReEt, ImEt = self.m.cff(pt)
        except AttributeError:
            ReH = self.m.ReH(pt)
            ImH = self.m.ImH(pt)
            ReE = self.m.ReE(pt)
            ImE = self.m.ImE(pt)
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

    def F2(self, pt):
        """Return DIS F2."""
        return self.m.F2(pt)

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
            res = quadrature.Hquadrature(lambda phi:
                    self._BSA(pt, vars={'phi':phi}) * sin(phi), 0, 2*pi)
        else:
            raise ValueError('[%s] has neither azimuthal angle phi\
 nor harmonic FTn == -1 defined!' % pt)
        return res / pi

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


class BMK(DVCS):
    """Implementation of formulas from hep-ph/0112108  (BMK)"""


    def PreFacBH(self, pt):
        """ Prefactor from Eq. (25), without e^6 """
        return 1./(pt.xB**2 * pt.y**2 * (1.+pt.eps2)**2 * pt.t * pt.P1P2)

    def PreFacDVCS(self, pt):
        """ Prefactor from Eq. (26), without e^6 """
        return 1./(pt.y**2 * pt.Q2 )

    def PreFacINT(self, pt):
        """ Prefactor from Eq. (27), without e^6 """
        return 1./(pt.xB * pt.y**3 * pt.t * pt.P1P2)


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
         (1. + ((t - tmin(Q2, xB, eps2))*(1. - xB + 0.5*(-1. + sqrt(1. + eps2))))/
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
    """According to BM arXiv:1005.5209 [hep-ph] - exact formulas.

    Only this Approach implements observables for longitudinally polarized target!
    """


    #### Bethe-Heitler - also longitudinally polarized target  (LP)

    def cBH0LP(self, pt):
        """ BKM Eq. (38) """
        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        FE = self.m.F1(pt) + t * self.m.F2(pt) / (4.0 * Mp2)
        FM = self.m.F1(pt) + self.m.F2(pt)
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
        FE = self.m.F1(pt) + t * self.m.F2(pt) / (4.0 * Mp2)
        FM = self.m.F1(pt) + self.m.F2(pt)
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
    self.m.F1(pt)*CFFH - (t/(4*Mp2))*self.m.F2(pt)*CFFE +
     (xB/((2 - xB) + (t*xB)/pt.Q2))*(self.m.F1(pt) + self.m.F2(pt))*CFFHt
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
    (xB/((2 - xB) + (t*xB)/pt.Q2))*(self.m.F1(pt) + self.m.F2(pt))*
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
    (xB/((2 - xB) + (t*xB)/pt.Q2))*(self.m.F1(pt) + self.m.F2(pt))*CFFHt
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
    (xB/((2 - xB) + (t*xB)/pt.Q2))*(self.m.F1(pt) + self.m.F2(pt))*
      (CFFH + (xB/2)*(1 - t/pt.Q2)*
        CFFE) + (1 + (Mp2/pt.Q2)*
        ((2*xB**2)/((2 - xB) + (t*xB)/pt.Q2))*
        (3 + t/pt.Q2))*self.m.F1(pt)*CFFHt -
     (t/pt.Q2)*((xB*(1 - 2*xB))/((2 - xB) +
        (t*xB)/pt.Q2))*self.m.F2(pt)*CFFHt -
     (xB/((2 - xB) + (t*xB)/pt.Q2))*
      ((xB/2)*(1 - t/pt.Q2)*self.m.F1(pt) + (t/(4*Mp2))*self.m.F2(pt))*
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
    (xB/((2 - xB) + (t*xB)/pt.Q2))*(self.m.F1(pt) + self.m.F2(pt))*
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
    (xB/((2 - xB) + (t*xB)/pt.Q2))*(self.m.F1(pt) + self.m.F2(pt))*
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
      (1 - ((t - tmin(Q2, xB, eps2))*(1 + (1 - eps2)/sqrt(1 + eps2) -
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
       (t*(t - tmin(Q2, xB, eps2))*xB*(1 - y - (y**2*eps2)/4)*
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
      (((t - tmin(Q2, xB, eps2))*(2 - y)**2*(1 - xB + ((t - tmin(Q2, xB, eps2))*((1 - xB)*xB +
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
          sqrt(1 + eps2))) + (t*(t - tmin(Q2, xB, eps2))*xB*
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
    (-4*t*(2 - y)*(((t - tmin(Q2, xB, eps2))*(1 - y - (y**2*eps2)/4)*
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
       ((t - tmin(Q2, xB, eps2))*(-2*xB**2 + eps2 +
          xB*(3 + sqrt(1 + eps2))))/(4*pt.Q2*
         sqrt(1 + eps2))))/(pt.Q2*
      (1 + eps2)**(3/2.))
    )

    CINT['unpA', (-1, 1), (2)] = CINTunpA212

    def CINTunpA213(self, pt):
        """Same as CINT["unpA", (-1, 1), (3)] """


        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (16*pt.K*t*(1 - y - (y**2*eps2)/4)*
      (1 - xB + ((t - tmin(Q2, xB, eps2))*((1 - xB)*xB + eps2/4))/
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
       (2*(t - tmin(Q2, xB, eps2))*((1 - xB)*xB + eps2/4))/(pt.Q2*
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
         ((t - tmin(Q2, xB, eps2))*(4*(1 - xB)*xB + eps2))/(4*pt.Q2*
           sqrt(1 + eps2))) - (2 - y)**2*(1 - xB/2 +
         ((t - tmin(Q2, xB, eps2))*(4*(1 - xB)*xB + eps2))/(2*pt.Q2*
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
       ((t - tmin(Q2, xB, eps2))*xB*(3 - 2*xB + eps2/xB - sqrt(1 + eps2)))/
        pt.Q2))/(pt.Q2*(1 + eps2)**2)
    )

    CINT['unpA', (1, 1), (2)] = CINTunpA112

    def CINTunpA113(self, pt):
        """Same as CINT["unpA", (1, 1), (3)] """


        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (16*pt.K*t*(t - tmin(Q2, xB, eps2))*((1 - xB)*xB + eps2/4)*
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
       ((t - tmin(Q2, xB, eps2))*(1 - y - (y**2*eps2)/4)*(1 - 2*xB +
          sqrt(1 + eps2)))/(2*pt.Q2)))/
     (pt.Q2*(1 + eps2)**(5/2.))
    )

    CINT['unpV', (1, 1), (1)] = CINTunpV111

    def CINTunpV112(self, pt):
        """Same as CINT["unpV", (1, 1), (2)] """


        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (8*t*xB*(2 - y)*(1 - y - (y**2*eps2)/4)*
      ((4*pt.tK2)/(pt.Q2*sqrt(1 + eps2)) +
       ((t - tmin(Q2, xB, eps2))*(1 + t/pt.Q2)*(1 - 2*xB +
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
        (xB*t)/pt.Q2)*((t - tmin(Q2, xB, eps2))/pt.Q2))
    )

    SINT['LP', (1, 1), (2)] = SINTLP112

    def SINTLP113(self, pt):
        """Same as SINT["LP", (1, 1), (3)] """


        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return (
    ((-8*pt.K*(1 - y - (y**2*eps2)/4))/
      (1 + eps2)**3)*((1 + sqrt(1 + eps2) - 2*xB)/
      (2*(sqrt(1 + eps2) + 1)))*eps2*
     ((t - tmin(Q2, xB, eps2))/pt.Q2)
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
       ((t - tmin(Q2, xB, eps2))*(1 - (t*(1 - 2*xB))/pt.Q2)*
         (1 - 2*xB + sqrt(1 + eps2)))/(2*pt.Q2))
      )/(pt.Q2*(1 + eps2)**3)
    )

    SINT['LPA', (1, 1), (2)] = SINTLPA112

    def SINTLPA113(self, pt):
        """Same as SINT["LPA", (1, 1), (3)] """


        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-8*pt.K*t*(t - tmin(Q2, xB, eps2))*xB*
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
      (1 - xB + ((t - tmin(Q2, xB, eps2))*((1 - xB)*xB + eps2/4))/
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
       (1 - ((t - tmin(Q2, xB, eps2))*(1 - 2*xB)*(1 - 2*xB + sqrt(1 + eps2)))/
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
       ((t - tmin(Q2, xB, eps2))*(-2*xB**2 + eps2 +
          xB*(3 - sqrt(1 + eps2))))/pt.Q2)
      )/(pt.Q2*(1 + eps2)**(5/2.))
    )

    SINT['LPV', (1, 1), (2)] = SINTLPV112

    def SINTLPV113(self, pt):
        """Same as SINT["LPV", (1, 1), (3)] """


        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (16*pt.K*t*(t - tmin(Q2, xB, eps2))*((1 - xB)*xB + eps2/4)*
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
    (8*pt.K*(2 - y)*y*(1 + ((t - tmin(Q2, xB, eps2))*(1 - xB + (-1 + sqrt(1 + eps2))/
           2))/(pt.Q2*(1 + eps2)))*pt.in1polarization)/
     (1 + eps2)
    )

    SINT['unp', (1, 1), (1)] = SINTunp111

    def SINTunp112(self, pt):
        """Same as SINT["unp", (1, 1), (2)] """


        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return ( (-4*(t - tmin(Q2, xB, eps2))*y*(1 - y - (y**2*eps2)/4)*
      (1 - 2*xB + sqrt(1 + eps2))*
      (-((t - tmin(Q2, xB, eps2))*(2*xB + eps2))/(2*pt.Q2*
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
         2*xB)/(2*(1 + eps2)))*((t - tmin(Q2, xB, eps2))/pt.Q2))
    )

    SINT['unpA', (1, 1), (1)] = SINTunpA111

    def SINTunpA112(self, pt):
        """Same as SINT["unpA", (1, 1), (2)] """


        xB, Q2, t, y, eps2  = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2
        return (
    (-((16*pt.in1polarization*(1 - y - (y**2*eps2)/4)*y)/(1 + eps2)**2))*
     (t/pt.Q2)*((t - tmin(Q2, xB, eps2))/pt.Q2)*
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
            self.m.F1(pt)*CFFH - (t/(4*Mp2))*self.m.F2(pt)*CFFE +
              xB/(2 - xB)*(self.m.F1(pt) + self.m.F2(pt))*CFFHt
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
        xB/(2 - xB)*(self.m.F1(pt) + self.m.F2(pt)) * (CFFH + CFFE)
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
    (xB/(2 - xB))*(self.m.F1(pt) + self.m.F2(pt))*CFFHt
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
          xB/(2 - xB)*(self.m.F1(pt) + self.m.F2(pt))*(CFFH + (xB/2)*CFFE)
          + self.m.F1(pt)*CFFHt
          - xB/(2 - xB)*((xB/2)*self.m.F1(pt)
                   + (t/(4*Mp2))*self.m.F2(pt))*CFFEt
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
          xB/(2 - xB)*(self.m.F1(pt) + self.m.F2(pt))*(CFFH + (xB/2)*CFFE)
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
         xB/(2 - xB)*(self.m.F1(pt) + self.m.F2(pt))*(CFFHt + (xB/2)*CFFEt)
              )
        if im:
            return res.imag
        else:
            return res.real

class BM10tw2(BM10):
    """According to BM arXiv:1005.5209 [hep-ph]

    This is BM10, but with C^{INT}_{2,3} set to zero for agreement
    with Dieter.
    """

    def cINT2unp(self, pt): return 0
    def cINT3unp(self, pt): return 0
    def sINT2unp(self, pt): return 0
    def sINT3unp(self, pt): return 0
    def cINT2LP(self, pt): return 0
    def cINT3LP(self, pt): return 0
    def sINT2LP(self, pt): return 0
    def sINT3LP(self, pt): return 0


class DVMP(Theory):
    """DVMP observables

    Implements cross-section for electroproduction of meson.
    """

    def _XrhotApprox(self, pt):
        """Partial DVrhoP cross section w.r.t. Mandelstam t.

        Approximate formula valid for small xB.

        """

        # 4 * pi**2 * alpha_em * GeV2nb = 112175.5
        res = 112175.5 * pt.xB**2 * (
                self.m.ImHrho(pt)**2 + self.m.ReHrho(pt)**2) / pt.Q2**2
        return res


    _Xrhot = _XrhotApprox


class Default(BMK, DVMP):
    """Default DVCS+DVMP class."""
