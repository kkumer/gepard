"""Deeply Virtual Compton Scattering (DVCS) observables."""

# Unlike DVMP and DIS, DVCS theory is separated into three files
#  dvcs.py  -  cross-section and derived observables
#   bmk.py  -  BMK formulas for interference, DVCS-squared and BH-squared terms
#   cff.py  -  models for Compton Form Factors

from numpy import cos, pi, sin, sqrt

from . import cff, data, quadrature, theory
from .constants import GeV2nb, Mp2, alpha
from .kinematics import weight_BH


class DVCS(theory.Theory):
    """DVCS observables base class.

    Implements cross-section for electroproduction of polarized real photon,
    on polarized nucleon target, as well as derived cross-sections and  asymmetries.
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

    def PreFacSigma(self, pt: data.DataPoint) -> float:
        """Overall prefactor of 4-fold XS.

        Args:
            pt: instance of DataPoint

        Returns:
            Overal prefactor from Eq. (22) of BMK,  times 2pi because of
            proton Phi integration and times y/Q2 because of dy -> dQ2.
            Converted to nanobarns.

        """
        return alpha**3 * pt.xB * pt.y**2 / (8. * pi * pt.Q2**2 *
                                             sqrt(1.+pt.eps2)) * GeV2nb

    def XS(self, pt: data.DataPoint, **kwargs) -> float:
        """Differential 4-fold e p --> e p gamma cross section.

        Args:
            pt: DataPoint instance
            **kwargs: keyword arguments

        Keyword args:
            vars (dict): dict which can override pt kinematics, e.g., vars={'phi': 0}
            flip (string, list): flip the sign of pt attribute,
                                 e.g., flip='in1polarization'

        Returns:
            Cross section.

        This is for calculation of fully beam and target polarized cross-section.
        Most of the usually measured cross-sections and asymmetries are derived
        from this one.

        Todo:
            * transversal target BMK/Trento switching

        """
        # We create new instance of data.Datapoint, which can be manipulated with
        # without messing up the original pt.

        # Overriding pt kinematics with those from 'vars' keyword attribute
        if 'vars' in kwargs:
            kin = data.DataPoint(kindict=kwargs['vars'])
            data._fill_kinematics(kin, old=pt)
        else:
            kin = pt.copy()

        # Copy non-kinematical info
        for atr in ['in1charge', 'in1polarizationvector', 'in1polarization',
                    'in2polarizationvector', 'in2polarization', 'in2particle',
                    'process', 'exptype', 'in1energy', 'in2energy', 'observable']:
            if atr in pt:
                setattr(kin, atr, getattr(pt, atr))

        # Pre-calculate some kinematical functions
        kin.prepare()

        # Flipping spins and/or charges for asymmetries, if asked for
        if 'flip' in kwargs and kwargs['flip']:
            if isinstance(kwargs['flip'], list):
                for item in kwargs['flip']:
                    setattr(kin, item, - getattr(pt, item))
            else:
                setattr(kin, kwargs['flip'], - getattr(pt, kwargs['flip']))

        # Weighting the integrand by BH propagators, if asked for
        if 'weighted' in kwargs and kwargs['weighted']:
            wgh = weight_BH(kin)
        else:
            wgh = 1

        # Finally, we build up the cross-section
        # 1. -- unpolarized target part --
        aux = self.TBH2unp(kin) + self.TINTunp(kin) + self.TDVCS2unp(kin)
        if hasattr(pt, 'in2polarizationvector'):
            # 2. -- longitudinally polarized target part --
            if pt.in2polarizationvector == 'L':
                aux += kin.in2polarization*(
                        self.TBH2LP(kin) + self.TINTLP(kin) + self.TDVCS2LP(kin))
            elif pt.in2polarizationvector == 'T':
                # 3. -- transversally polarized target part --
                # We directly take cos(varphi) or sin(varphi) terms depending if
                # varFTn is specified
                if hasattr(pt, 'varFTn'):
                    # FIXME: should deal properly with Trento/BKM differences
                    # Also, this gives result in Trento convention now,
                    # and is completely O.K. for Trento asymmetries
                    # kin.varphi = (3*pt.varFTn+1)*pi/4.  # pi for cos, -pi/2 for sin
                    # And this is for BMK
                    kin.varphi = (1-pt.varFTn)*pi/4.  # 0 for cos, pi/2 for sin
                aux += kin.in2polarization*(
                        self.TBH2TP(kin) + self.TINTTP(kin) + self.TDVCS2TP(kin))
            elif pt.in2polarizationvector == 'U':
                pass
            else:
                raise ValueError('in2polarizationvector must be L, T, or U!')
        return wgh * self.PreFacSigma(kin) * aux

    def _XGAMMA_DVCS_t_Approx(self, pt):
        """Partial DVCS (gamma* p -> gamma p) cross section differential in t.

        Args:
            pt: instance of DataPoint

        Returns:
            Cross section. Approx. formula used in NPB10 paper.

        """
        W2 = pt.W * pt.W
        # Simplified formula used also in Fortran gepard code
        ReH, ImH, ReE, ImE, ReHt, ImHt, ReEt, ImEt = self.m.cff(pt)
        res = 260.5633976788416 * W2 * (
                (ImH**2 + ReH**2)
                - pt.t/(4.*Mp2)*(ReE**2 + ImE**2)) / (
            (W2 + pt.Q2) * (2.0 * W2 + pt.Q2)**2)
        return res

    def _XGAMMA_DVCS_t_Ex(self, pt):
        """Partial DVCS (gamma* p -> gamma p) cross section differential in t."""
        eps2 = 4. * pt.xB**2 * Mp2 / pt.Q2
        if cff.HybridCFF in self.__class__.mro():
            # For hybrid models we cannot ask for cff() since self gets misinterpreted
            ReH = self.m.ReH(pt)
            ImH = self.m.ImH(pt)
            ReE = self.m.ReE(pt)
            ImE = self.m.ImE(pt)
        else:
            ReH, ImH, ReE, ImE, ReHt, ImHt, ReEt, ImEt = self.m.cff(pt)
        res = 65.14079453579676 * (pt.xB**2 / pt.Q2**2 / (1-pt.xB) / (2-pt.xB)**2 /
                                   sqrt(1 + eps2) * (
                    4 * (1 - pt.xB) * (ImH**2 + ReH**2)
                    - pt.xB**2 * (ReE**2+ImE**2 + 2*ReE*ReH + 2*ImE*ImH)
                    - (2-pt.xB)**2 * pt.t/4/Mp2*(ImE**2+ReE**2)))
        return res

    # Choice of what is actually used:
    _XGAMMA_DVCS_t = _XGAMMA_DVCS_t_Ex

    def XSintphi(self, pt: data.DataPoint, **kwargs) -> float:
        """XS integrated over azimuthal angle.

        Args:
            pt: instance of DataPoint
            **kwargs: keyword arguments

        Returns:
            Cross section (unpolarized)

        Todo:
            This is unnecessary function. Should be changed to XUU.

        """
        mem = pt.__dict__.pop('FTn', None)
        pt.FTn = 0
        res = self.XUU(pt, **kwargs)
        # restore old value if needed
        if mem:
            pt.in2polarization = mem
        return 2*pi*res

    def _phiharmonic(self, fun, pt, **kwargs):
        """Return fun evaluated for phi=pt.phi, or proper harmonic."""
        if 'phi' in pt or ('vars' in kwargs
                           and 'phi' in kwargs['vars']):
            return fun(pt, **kwargs)
        elif 'FTn' in pt:
            if pt.FTn < 0:
                res = quadrature.Hquadrature(
                          lambda phi: fun(pt, vars={'phi': phi}, **kwargs)
                          * sin(-pt.FTn*phi), 0, 2*pi)
            elif pt.FTn > 0:
                res = quadrature.Hquadrature(
                          lambda phi: fun(pt, vars={'phi': phi}, **kwargs)
                          * cos(pt.FTn*phi), 0, 2*pi)
            elif pt.FTn == 0:
                res = quadrature.Hquadrature(
                          lambda phi: fun(pt, vars={'phi': phi}, **kwargs), 0, 2*pi)/2
            else:
                raise ValueError('FTn = {} looks wrong!'.format(str(pt.FTn)))
            return res / pi
        else:
            raise ValueError(
              '{} has neither azimuthal angle phi nor harmonic FTn defined!'.format(pt))

# Observables: cross-sections

    def _XUU(self, pt, **kwargs):
        """4-fold beam helicity-independent cross section."""
        R = kwargs.copy()
        R.update({'flip': 'in1polarization'})
        return (self.XS(pt, **kwargs) + self.XS(pt, **R)) / 2

    def XUU(self, pt, **kwargs):
        """4-fold beam helicity-independent cross section."""
        return self._phiharmonic(self._XUU, pt, **kwargs)

    def XUUw(self, pt, **kwargs):
        """Weighted 4-fold beam helicity-independent cross section."""
        kwargs['weighted'] = True
        return self._phiharmonic(self.XUU, pt, **kwargs)

    def _XLU(self, pt, **kwargs):
        """4-fold beam helicity-dependent cross section."""
        R = kwargs.copy()
        R.update({'flip': 'in1polarization'})
        return (self.XS(pt, **kwargs) - self.XS(pt, **R)) / 2

    def XLU(self, pt, **kwargs):
        """4-fold beam helicity-dependent cross section."""
        return self._phiharmonic(self._XLU, pt, **kwargs)

    def XLUw(self, pt, **kwargs):
        """Weighted 4-fold beam helicity-dependent cross section."""
        kwargs['weighted'] = True
        return self._phiharmonic(self.XLU, pt, **kwargs)

    def _XUD(self, pt, **kwargs):
        """Cross-section difference on (long. or trans.) polarized target."""
        pol = kwargs.copy()
        pol.update({'flip': 'in2polarization'})
        o = self.XS(pt, **kwargs)
        f = self.XS(pt, **pol)
        return (o-f)/2

    def XUL(self, pt, **kwargs):
        """Cross-section difference on longitudinally polarized target."""
        assert pt.in2polarizationvector == 'L'
        return self._XUD(pt, **kwargs)

    def XUT(self, pt, **kwargs):
        """Cross-section difference on transversally polarized target."""
        assert pt.in2polarizationvector == 'T'
        return self._XUD(pt, **kwargs)

    def _XCLU(self, pt, **kwargs):
        """4-fold beam charge-spin cross section difference.

        Args:
            pt: instance of DataPoint
            **kwargs: keyword arguments

        Returns:
            D_{CS,U}/2, as measured by COMPASS.

        """
        R = kwargs.copy()
        R.update({'flip': ['in1polarization', 'in1charge']})
        return (self.XS(pt, **kwargs) - self.XS(pt, **R)) / 2

    def XCLU(self, pt, **kwargs):
        """4-fold beam charge-spin cross section difference.

        Args:
            pt: instance of DataPoint
            **kwargs: keyword arguments

        Returns:
            D_{CS,U}/2, as measured by COMPASS. (NOT TESTED.)

        """
        return self._phiharmonic(self._XCLU, pt, **kwargs)

    def _XCUU(self, pt, **kwargs):
        """4-fold beam charge-spin cross section sum.

        Args:
            pt: instance of DataPoint
            **kwargs: keyword arguments

        Returns:
            S_{CS,U}/2, as measured by COMPASS.

        """
        R = kwargs.copy()
        R.update({'flip': ['in1polarization', 'in1charge']})
        return (self.XS(pt, **kwargs) + self.XS(pt, **R)) / 2

    def XCUU(self, pt, **kwargs):
        """4-fold beam charge-spin cross section sum.

        Args:
            pt: instance of DataPoint
            **kwargs: keyword arguments

        Returns:
            S_{CS,U}/2, as measured by COMPASS. (NOT TESTED.)

        """
        return self._phiharmonic(self._XCUU, pt, **kwargs)

# Observables:  asymmetries

    def _AC(self, pt, **kwargs):
        """Calculate beam charge asymmtery."""
        chg = kwargs.copy()
        chg.update({'flip': 'in1charge'})
        o = self.XS(pt, **kwargs)
        c = self.XS(pt, **chg)
        return (o - c) / (o + c)
        # optimized formula (don't calculate parts which cancel anyway)
        # return  self.TINTunp(pt, phi, 0, 1) / (
        #               self.TBH2unp(pt, phi) + self.TDVCS2unp(pt, phi) )

    def AC(self, pt, **kwargs):
        """Calculate beam charge asymmetry."""
        res = self._phiharmonic(self._AC, pt, **kwargs)
        return res

    def _ALU(self, pt, **kwargs):
        """Calculate beam spin asymmetry."""
        return self._XLU(pt, **kwargs) / self._XUU(pt, **kwargs)

    def _ALUapprox(self, pt, **kwargs):
        """Calculate beam spin asymmetry.

        Args:
            pt: instance of DataPoint
            **kwargs: keyword arguments

        Returns:
            Beam spin asymmetry A_LU.
            For harmonics, fast approximate formula is used.

        """
        if 'phi' in pt:
            return self._ALU(pt, **kwargs)
        elif 'FTn' in pt and pt.FTn == -1:
            if 'vars' in kwargs:
                kwargs['vars'].update({'phi': pi/2.})
            else:
                kwargs['vars'] = {'phi': pi/2.}
            return self._ALU(pt, **kwargs)
        else:
            raise ValueError('[%s] has neither azimuthal angle phi\
 nor harmonic FTn = -1 defined!' % pt)

    def _ALUexact(self, pt, **kwargs):
        """Calculate beam spin asymmetry."""
        res = self._phiharmonic(self._ALU, pt, **kwargs)
        return res

    # Make a choice
    ALU = _ALUexact

    def _TSA(self, pt, **kwargs):
        """Target spin asymmetry (transversal or longitudinal)."""
        pol = kwargs.copy()
        pol.update({'flip': 'in2polarization'})
        o = self.XS(pt, **kwargs)
        p = self.XS(pt, **pol)
        return (o-p)/(o+p)

    def TSA(self, pt, **kwargs):
        """Calculate target spin asymmetry (AUL or AUT).

        Args:
            pt: instance of DataPoint
            **kwargs: keyword arguments

        Returns:
            Target spin asymmetry.
            Whether AUL or AUT is calculated is determined by the
            pt.in2polarizationvector being 'L' or 'T'.

        """
        return self._phiharmonic(self._TSA, pt, **kwargs)

    def AUL(self, pt, **kwargs):
        """Calculate longitudinal target spin asymmetry."""
        assert pt.in2polarizationvector == 'L'
        return self.TSA(pt, **kwargs)

    def AUT(self, pt, **kwargs):
        """Calculate transversal target spin asymmetry."""
        assert pt.in2polarizationvector == 'T'
        return self.TSA(pt, **kwargs)

    def _BTSA(self, pt, **kwargs):
        """Calculate beam-target spin asymmetry (BTSA)."""
        # According to 1004.0177 Eq. (1.8)
        bpol = kwargs.copy()
        bpol.update({'flip': 'in1polarization'})
        tpol = kwargs.copy()
        tpol.update({'flip': 'in2polarization'})
        both = kwargs.copy()
        both.update({'flip': ['in1polarization', 'in2polarization']})
        o = self.XS(pt, **kwargs)
        p = self.XS(pt, **bpol)
        t = self.XS(pt, **tpol)
        b = self.XS(pt, **both)
        return ((o+b) - (p+t)) / ((o+b) + (p+t))

    def BTSA(self, pt, **kwargs):
        """Calculate beam-target spin asymmetry (ALL or ALT).

        Args:
            pt: instance of DataPoint
            **kwargs: keyword arguments

        Returns:
            Beam-target spin asymmetry.
            Whether ALL or ALT is calculated is determined by the
            pt.in2polarizationvector being 'L' or 'T'.

        """
        return self._phiharmonic(self._BTSA, pt, **kwargs)

    def ALL(self, pt, **kwargs):
        """Calculate beam and longitudinal target spin asymmetry."""
        assert pt.in2polarizationvector == 'L'
        return self.BTSA(pt, **kwargs)

    def ALT(self, pt, **kwargs):
        """Calculate beam and transversal target spin asymmetry."""
        assert pt.in2polarizationvector == 'T'
        return self.BTSA(pt, **kwargs)

    def _CBTSA(self, pt, chargepar=-1, **kwargs):
        """Calculate double charge and spin asymmetries."""
        # According to 1106.2990 Eq. (18)(chargepar=1) and (19)(chargepar=-1)
        T_args = kwargs.copy()
        T_args.update({'flip': 'in2polarization'})
        B_args = kwargs.copy()
        B_args.update({'flip': 'in1polarization'})
        BT_args = kwargs.copy()
        BT_args.update({'flip': ['in1polarization', 'in2polarization']})
        C_args = kwargs.copy()
        C_args.update({'flip': 'in1charge'})
        CT_args = kwargs.copy()
        CT_args.update({'flip': ['in1charge', 'in2polarization']})
        CB_args = kwargs.copy()
        CB_args.update({'flip': ['in1charge', 'in1polarization']})
        CBT_args = kwargs.copy()
        CBT_args.update({'flip': ['in1charge', 'in1polarization', 'in2polarization']})
        o = self.XS(pt, **kwargs)
        t = self.XS(pt, **T_args)
        b = self.XS(pt, **B_args)
        bt = self.XS(pt, **BT_args)
        c = self.XS(pt, **C_args)
        ct = self.XS(pt, **CT_args)
        cb = self.XS(pt, **CB_args)
        cbt = self.XS(pt, **CBT_args)
        return ((o-t-b+bt) + chargepar*(c-ct-cb+cbt)) / (o+t+b+bt+c+ct+cb+cbt)

    def ALTI(self, pt, **kwargs):
        """Calculate {A_LT,I} as defined by HERMES 1106.2990 Eq. (19)."""
        return self._phiharmonic(self._CBTSA, pt, **kwargs)

    def ALTBHDVCS(self, pt, **kwargs):
        """Calculate {A_LT,BHDVCS} as defined by HERMES 1106.2990 Eq. (18)."""
        return self._phiharmonic(self._CBTSA, pt, chargepar=+1, **kwargs)

    def _ALUI(self, pt, **kwargs):
        """Beam spin asymmetry, interference part."""
        pol = kwargs.copy()
        pol.update({'flip': 'in1polarization'})
        chg = kwargs.copy()
        chg.update({'flip': 'in1charge'})
        both = kwargs.copy()
        both.update({'flip': ['in1polarization', 'in1charge']})
        o = self.XS(pt, **kwargs)
        p = self.XS(pt, **pol)
        c = self.XS(pt, **chg)
        b = self.XS(pt, **both)
        return ((o-p) - (c-b)) / ((o+p) + (c+b))

    def ALUI(self, pt, **kwargs):
        """Calculate beam spin asymmetry, interference part."""
        # As defined by HERMES 0909.3587 Eq. (2.2).
        return self._phiharmonic(self._ALUI, pt, **kwargs)

    def _ALUDVCS(self, pt, **kwargs):
        """Calculate ALU as defined by HERMES 0909.3587 Eq. (2.3)."""
        pol = kwargs.copy()
        pol.update({'flip': 'in1polarization'})
        chg = kwargs.copy()
        chg.update({'flip': 'in1charge'})
        both = kwargs.copy()
        both.update({'flip': ['in1polarization', 'in1charge']})
        o = self.XS(pt, **kwargs)
        p = self.XS(pt, **pol)
        c = self.XS(pt, **chg)
        b = self.XS(pt, **both)
        return ((o-p) + (c-b)) / ((o+p) + (c+b))

    def ALUDVCS(self, pt, **kwargs):
        """Calculate beam spin asymmetry, DVCS part."""
        # As defined by HERMES 0909.3587 Eq. (2.3).
        return self._phiharmonic(self._ALUDVCS, pt, **kwargs)

    def _AUTI(self, pt, **kwargs):
        """Calculate trans. target asymmetry, defined by HERMES 0802.2499 Eq. (15)."""
        pol = kwargs.copy()
        pol.update({'flip': 'in2polarization'})
        chg = kwargs.copy()
        chg.update({'flip': 'in1charge'})
        both = kwargs.copy()
        both.update({'flip': ['in2polarization', 'in1charge']})
        o = self.XS(pt, **kwargs)
        p = self.XS(pt, **pol)
        c = self.XS(pt, **chg)
        b = self.XS(pt, **both)
        return ((o-p) - (c-b)) / ((o+p) + (c+b))

    def AUTI(self, pt, **kwargs):
        """Calculate trans. target asymmetry, defined by HERMES 0802.2499 Eq. (15)."""
        return self._phiharmonic(self._AUTI, pt, **kwargs)

    def _AUTDVCS(self, pt, **kwargs):
        """Calculate trans. target asymmetry, defined by HERMES 0802.2499 Eq. (14)."""
        pol = kwargs.copy()
        pol.update({'flip': 'in2polarization'})
        chg = kwargs.copy()
        chg.update({'flip': 'in1charge'})
        both = kwargs.copy()
        both.update({'flip': ['in2polarization', 'in1charge']})
        o = self.XS(pt, **kwargs)
        p = self.XS(pt, **pol)
        c = self.XS(pt, **chg)
        b = self.XS(pt, **both)
        return ((o-p) + (c-b)) / ((o+p) + (c+b))

    def AUTDVCS(self, pt, **kwargs):
        """Calculate trans. target asymmetry, defined by HERMES 0802.2499 Eq. (14)."""
        return self._phiharmonic(self._AUTDVCS, pt, **kwargs)

# Observables:  ad-hoc, one-off stuff

    def AC0minusr1(self, pt):
        """Return some combination liked by Dieter."""
        return self.ACcos0(pt) - pt.r * self.ACcos1(pt)

    def XLUw2C(self, pt):
        """Im(C^I) as defined by HALL A."""
        return self.ImCCALINTunp(pt)

    def XUUw2C(self, pt):
        """Re(C^I) or Re(C^I + Del C^I) as defined by HALL A.

        Args:
            pt: instance of DataPoint

        Returns:
            Harmonic coefficients.

        Todo:
            Although it is attributed to FTn=0, Re(C^I + Del C^I)
            is only a part of zeroth harmonic.

        """
        if pt.FTn == 0:
            return self.ReCCALINTunp(pt) + self.ReDELCCALINTunp(pt)
        elif pt.FTn == 1:
            return self.ReCCALINTunp(pt)

    def XwA(self, pt):
        """Ratio of first two cos harmonics of w-weighted cross section."""
        # In BMK, not Trento?
        b0 = quadrature.Hquadrature(
                lambda phi: self.XUU(pt, vars={'phi': phi}, weighted=True),
                0, 2.0*pi) / (2.0*pi)
        b1 = quadrature.Hquadrature(
                lambda phi: self.XUU(pt, vars={'phi': phi}, weighted=True) * cos(phi),
                0, 2.0*pi) / pi
        return b1/b0

    def BCSA(self, pt, **kwargs):
        """Beam charge-spin asymmetry as measured by COMPASS."""
        return self.BCSD(pt, **kwargs) / self.BCSS(pt, **kwargs)


# This PEP8 violating end-of-file import serves just to bring BMK into dvcs namespace
from .bmk import BM10, BMK, BM10ex, BM10tw2, hotfixedBMK  # noqa: F401, E402
