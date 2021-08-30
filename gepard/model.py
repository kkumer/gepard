"""Definitions of models.

Models of "soft" hadronic structure functions:
    * Elastic electromagnetic form factors
    * Generalized parton distribution functions (GPDs)
    * Compton form factors (for DVCS)
"""

from math import log, pi
from typing import Dict, List

import numpy as np
from joblib import Parallel, delayed

import gepard as g


class Model(object):
    """Base class for all models.

    Instance of this class specifies structure of relevant hadrons.
    Methods provided are typically GPDs, CFFs, elastic FFs,
    DVMP transition TFFs, DAs, etc.
    Main subclasses are:

     - ParameterModel which depends on real parameters (some of which can be
       provided by minimization routine in the fitting procedure).
     - NeuralModel where structure functions are represented as neural nets
       (not yet implemented)
     - GridModel where values of structure functions are represented as
       grids of numbers, which may be interpolated
       (not yet implemented)

    These are then further subclassed to model actual structure functions.

    """
    def __init__(self, name: str = None, texname: str = None,
                 description: str = None) -> None:
        """Init Model class.

        Args:
            name: short unique model name
            texname: TeX model name for (e.g. for plot annotation)
            description: longer description of the model

        """
        self.name = 'N/A'
        if not texname:
            self.texname = name
        else:
            self.texname = texname
        self.description = 'N/A'


class ParameterModel(Model):
    """Base class for all parametrized models.

    Attributes:
        parameters: dict {'par0': float, ...}. Actual value of parameter.
        parameters_fix: dict {'par0': bool, ...}. Is parameter value
            fixed? Considered False if non-existent.
        parameters_limit: dict {'par0: (float, float), ...}. Range for
            fitting. Considered (-inf, inf) if non-existent.
    """
    parameters: dict = {}
    parameters_fix: dict = {}
    parameters_limit: dict = {}

    def __init__(self) -> None:
        """Init and pre-calculate stuff."""
        super().__init__()

    def release_parameters(self, *pars: str):
        """Release parameters for fitting.

        Args:
            *pars: Names of parameters to be released

        Notes:
            If allowed parameter ranges have to be changed from default, user needs
            to modify parameters dictionary directly.

        Todo:
            Note that relasing and fixing parameters for model instance is sensible
            only before creating Fitter instance! Afterwards, one has to fix and release
            params both for Fitter and for model.

        """
        for par in pars:
            if par not in self.parameters:
                raise ValueError('Parameter {} is not defined in model {}'.format(
                        par, self))
            self.parameters_fix[par] = False

    def fix_parameters(self, *pars: str):
        """Fix parameters so they are not fitting variables.

        Args:
            *pars: Names of parameters to be fixed. If first name is 'ALL'
                then fix all parameters.

        Todo:
            Note that relasing and fixing parameters for model instance is sensible
            only before creating Fitter instance! Afterwards, one has to fix and release
            params both for Fitter and for model.

        """
        if pars[0] == 'ALL':
            # fix 'em all
            for par in self.parameters:
                self.parameters_fix[par] = True
        else:
            for par in pars:
                if par not in self.parameters:
                    raise ValueError('Parameter {} is not defined in model {}'.format(
                            par, self))
                self.parameters_fix[par] = True

    def free_parameters(self) -> List[str]:
        """Return list of names of free fitting parameters."""
        return [p for p in self.parameters if p not in self.parameters_fix
                or not self.parameters_fix[p]]


# --- Models for Elastic Form Factors --- #

class ElasticFormFactors(Model):
    """Dirac and Pauli elastic form factors F_1 and F_2."""


class ElasticDipole(ElasticFormFactors):
    """Dipole approximation from DM's notebook."""

    def F1(self, pt):
        """Dirac elastic proton form factor - dipole parametrization."""
        t = pt.t
        if 'in2particle' in pt and pt.in2particle == 'p':
            return (1.41 * (1.26 - t))/((0.71 - t)**2 * (3.53 - t))
        else:
            print('Neutron dipole elastic FFs are not implemented yet! Use Kelly.')

    def F2(self, pt):
        """Pauli elastic proton form factor - dipole parametrization."""
        t = pt.t
        if 'in2particle' in pt and pt.in2particle == 'p':
            return 3.2 / ((0.71 - t)**2 * (3.53 - t))
        else:
            print('Neutron dipole elastic FFs are not implemented yet! Use Kelly.')


class ElasticKelly(ElasticFormFactors):
    """Kelly's approximation from DM's notebook."""

    def F1(self, pt):
        """Dirac elastic nucleon form factor - Kelly's parametrization."""
        if 'in2particle' in pt and pt.in2particle == 'n':
            return self._nF1(pt.t)
        else:   # proton is default
            return self._pF1(pt.t)

    def F2(self, pt):
        """Dirac elastic nucleon form factor - Kelly's parametrization."""
        if 'in2particle' in pt and pt.in2particle == 'n':
            return self._nF2(pt.t)
        else:   # proton is default
            return self._pF2(pt.t)

    def _pF1(self, t):
        """Dirac elastic proton form factor - Kelly's parametrization."""
        return ((1 + 0.06815437285120148*t)/(1 - 3.118062557942468*t +
                1.0338391956016382*t**2 - 0.5031268669574522*t**3) -
                (0.7931031653189349*(1 - 0.03407718642560074*t)*t) /
                (1 - 3.115222792407001*t + 1.520921000705686*t**2 -
                0.14999913420898098*t**3))/(1 - 0.28397655354667284*t)

    def _pF2(self, t):
        """Pauli elastic proton form factor - Kelly's parametrization."""
        return (-((1 + 0.06815437285120148*t)/(1 - 3.118062557942468*t +
                1.0338391956016382*t**2 - 0.5031268669574522*t**3)) +
                (2.792847351*(1 - 0.03407718642560074*t))/(1 -
                3.115222792407001*t + 1.520921000705686*t**2 -
                0.14999913420898098*t**3)) / (1 - 0.28397655354667284*t)

    def _nF1(self, t):
        """Dirac elastic neutron form factor - Kelly's parametrization."""
        return ((-0.4842637275288574*t)/((1 - 1.4084507042253522*t)**2 *
                (1 - 0.9345440355820054*t)) +
                (0.5417644379086957*(1 - 0.6598447281533554*t)*t) /
                (1 - 4.168632789020339*t + 1.9408278987597791*t**2 -
                1.9100884849907935*t**3))/(1 - 0.2831951622975774*t)

    def _nF2(self, t):
        """Pauli elastic neutron form factor - Kelly's parametrization."""
        return ((0.4842637275288574*t)/((1 - 1.4084507042253522*t)**2 *
                (1 - 0.9345440355820054*t)) - (1.9130427*(1 - 0.6598447281533554*t)) /
                (1 - 4.168632789020339*t + 1.9408278987597791*t**2 -
                1.9100884849907935*t**3))/(1 - 0.2831951622975774*t)


class ElasticZero(ElasticFormFactors):
    """Set F1=F2=0 to get just DVCS^2."""

    def F1(self, pt):
        """Elastic em Dirac form factor F1 set to zero."""
        return 0.

    def F2(self, pt):
        """Elastic em Pauli form factor F2 set to zero."""
        return 0.


# --- Models for GPDs --- #


class ConformalSpaceGPD(ParameterModel):
    """Base class of GPD models built in conformal moment space."""

    def __init__(self, p: int = 0, scheme: str = 'msbar', nf: int = 4,
            q02: float = 4.0, r20: float = 2.5,
            asp: np.array = np.array([0.0606, 0.0518, 0.0488]),
            c: float = 0.35, phi: float = 1.57079632) -> None:
        """Init and pre-calculate stuff.

        Args:
            p: pQCD order (0 = LO, 1 = NLO, 2 = NNLO)
            scheme: pQCD scheme  (msbar or csbar)
            nf: number of active quark flavors
            q02: Initial Q0^2 for pQCD evolution.
            r20: Initial mu0^2 for alpha_strong definition.
            asp: alpha_strong/(2*pi) at scale r20 for (LO, NLO, NNLO)
            c: intersection of Mellin-Barnes curve with real axis
            phi: angle of Mellin-Barnes curve with real axis

        Notes:
            This just takes care of initialization of Mellin-Barnes
            contour points, and evolved Wilson coeffs, and MB Gauss
            integration weights. Actual choice and code for GPDs
            is provided by subclasses.

        """
        self.p = p
        self.scheme = scheme
        self.nf = nf
        self.q02 = q02
        self.r20 = r20
        self.asp = asp
        self.c = c
        self.phi = phi
        npoints, weights = g.quadrature.mellin_barnes(self.c, self.phi)
        self.npts = len(npoints)
        self.npoints = npoints
        self.jpoints = npoints - 1
        self.wg = weights  # Gauss integration weights
        # Initial parameters:
        self.parameters = {'ns': 2./3. - 0.4,
                           'al0s': 1.1,
                           'alps': 0.25,
                           'ms': 1.1,
                           'secs': 0.,
                           'this': 0.,
                           'kaps': 0.,
                           'ng': 0.4,
                           'al0g': 1.2,
                           'alpg': 0.25,
                           'mg': 1.2,
                           'secg': 0.,
                           'thig': 0.,
                           'kapg': 0.}
        super().__init__()

    def pw_strengths(self):
        """Strengths of SO(3) partial waves."""
        # We take maximally three partial waves atm:
        # pw_strengths = (no. pws x no. flavors)
        return np.array([[1., 1., 1, 1],
                         [self.parameters['secs'],
                             self.parameters['secg'], 0, 0],
                         [self.parameters['this'],
                             self.parameters['thig'], 0, 0]])

    def gpd_H(self, eta: float, t: float) -> np.ndarray:
        """Return (npts, 4) array H_j^a for all j-points and 4 flavors."""
        return np.zeros((self.npts, 4), dtype=complex)

    def gpd_E(self, eta: float, t: float) -> np.ndarray:
        """Return (npts, 4) array E_j^a for all j-points and 4 flavors."""
        return np.zeros((self.npts, 4), dtype=complex)


class Test(ConformalSpaceGPD):
    """Simple testing ansatz for GPDs."""

    def __init__(self, **kwargs) -> None:
        """See parent `ConformalSpaceGPD` class for docs."""
        kwargs.setdefault('scheme', 'csbar')
        kwargs.setdefault('nf', 3)
        kwargs.setdefault('q02', 1.0)
        super().__init__(**kwargs)
        self.nf = 3
        self.q02 = 1.0
        self.asp = np.array([0.05, 0.05, 0.05])
        self.r20 = 2.5

    def gpd_H(self, eta: float, t: float) -> np.ndarray:
        """Return (npts, 4) array H_j^a for all j-points and 4 flavors."""
        # For testing purposes, we use here sub-optimal non-numpy algorithm
        h = []
        for j in self.jpoints:
            h.append(g.gpdj.test(j, t, self.parameters))
        return np.array(h)


class FitBP(ConformalSpaceGPD):
    """GPD ansatz from paper hep-ph/0703179."""

    def __init__(self, **kwargs) -> None:
        """See parent `ConformalSpaceGPD` class for docs."""
        kwargs.setdefault('scheme', 'msbar')
        kwargs.setdefault('nf', 4)
        kwargs.setdefault('q02', 2.5)
        kwargs.setdefault('asp', np.array([0.05, 0.05, 0.05]))
        kwargs.setdefault('r20', 2.5)
        kwargs.setdefault('phi', 1.9)
        super().__init__(**kwargs)

    def gpd_H(self, eta: float, t: float) -> np.ndarray:
        """Return (npts, 4) array H_j^a for all j-points and 4 flavors."""
        # For testing purposes, we use here sub-optimal non-numpy algorithm
        h = []
        for j in self.jpoints:
            h.append(g.gpdj.fitbp(j, t, self.parameters))
        return np.array(h)


class Fit(ConformalSpaceGPD):
    """Singlet fitting ansatz for GPDs with three SO(3) partial waves."""

    def __init__(self, **kwargs) -> None:
        """See parent `ConformalSpaceGPD` class for docs."""
        kwargs.setdefault('scheme', 'msbar')
        kwargs.setdefault('nf', 4)
        kwargs.setdefault('q02', 4.0)
        super().__init__(**kwargs)

    def gpd_H(self, eta: float, t: float) -> np.ndarray:
        """Return (npts, 4) array H_j^a for all j-points and 4 flavors."""
        return g.gpdj.fit(self.jpoints, t, self.parameters).transpose()

    def gpd_H_para(self, eta: float, t: float) -> np.ndarray:
        """Return (npts, 4) array H_j^a for all j-points and 4 flavors."""
        # This multiprocessing version is actually 2x slower!
        h = Parallel(n_jobs=20)(delayed(g.gpdj.fit)(j, t, self.parameters)
                                for j in self.jpoints)
        return np.array(h)

    def gpd_E(self, eta: float, t: float) -> np.ndarray:
        """Return (npts, 4) array E_j^a for all j-points and 4 flavors."""
        kappa = np.array([self.parameters['kaps'], self.parameters['kapg'], 0, 0])
        return kappa * g.gpdj.fit(self.jpoints, t, self.parameters).transpose()


# --- Models for Compton Form Factors --- #

class ComptonFormFactors(ParameterModel):
    """Base class for CFFs.

    CFFs are set to be zero here. Actual models are built by subclassing this.
    """

    allCFFs = ['ImH', 'ReH', 'ImE', 'ReE', 'ImHt', 'ReHt', 'ImEt', 'ReEt']
    allCFFsb = ['ImH', 'ReH', 'ImE', 'ReE', 'ImHt', 'ReHt', 'ImEt', 'ReEb']
    allCFFeffs = ['ImHeff', 'ReHeff', 'ImEeff', 'ReEeff',
                  'ImHteff', 'ReHteff', 'ImEteff', 'ReEteff']
    # allGPDs = []

    def print_CFFs(self, pt, format=None):
        """Print values of CFFs at given kinematic point."""
        vals = [getattr(self, cff)(pt) for cff in self.allCFFs]
        if format == 'mma':
            s = "{" + 8*"%s -> %f, "
            s = s[:-2] + "}"
        else:
            s = 8*"%4s = %5.2f\n"
        print(s % g.utils.flatten(tuple(zip(self.allCFFs, vals))))

    # Initial definition of all CFFs. All just return zero.
    for name in allCFFs:
        exec('def %s(self, pt): return 0.' % name)

    for name in allCFFeffs:
        exec('def %s(self, pt): return 0.' % name)

    # Define E-bar as xi*E-tilde
    def ReEb(self, pt):
        return (pt.xi)*self.ReEt(pt)

    def is_within_model_kinematics(self, pt):
        """Check if kinematics is sensible."""
        return ((1.5 <= pt.Q2 <= 5.) and
                (pt.tm < min(1., pt.Q2/4)) and
                (1e-3 < pt.xB < 0.5))


class MellinBarnesModel(ParameterModel):
    """Class of models built by Mellin-Barnes integration."""

    def __init__(self, gpds: ConformalSpaceGPD) -> None:
        """Init MellinBarnes class and pre-calculate stuff.

        Args:
            gpds: provide GPDs as gpds.gpd_H, gpds.gpd_Ht etc.

        """
        self.p = gpds.p
        self.scheme = gpds.scheme
        self.nf = gpds.nf
        self.q02 = gpds.q02
        self.asp = gpds.asp
        self.r20 = gpds.r20
        self.npts = gpds.npts
        self.npoints = gpds.npoints
        self.jpoints = gpds.jpoints
        self.c = gpds.c
        self.phi = gpds.phi
        self.wg = gpds.wg
        self.gpds = gpds
        self.parameters = gpds.parameters
        # Consolidate parameters, both the same and updated from above
        self.gpds.parameters = self.parameters
        # scales
        self.rr2 = 1     # ratio of Q2/renorm. scale squared
        self.rf2 = 1     # ratio of Q2/GPD fact. scale sq.
        self.rdaf2 = 1   # ratio of Q2/DA fact. scale sq. (for DVMP)
        #
        self.tgj = np.tan(pi*self.jpoints/2.)
        # wce[Q2] = wce[spw, j, a] - Wilson coeffs evolved; local to model instance
        self.wce: Dict[float, np.ndarray] = {}  # DVCS
        self.wce_dvmp: Dict[float, np.ndarray] = {}  # DVMP
        self.wce_dis: Dict[float, np.ndarray] = {}  # DIS
        super().__init__()

    def _mellin_barnes_integral(self, xi, wce, h):
        """Return Mellin-Barnes integral.

        Integrates convolution of evolved GPDs with Wilson coefs.

        Todo:
            Only GPD H treated atm.
        """
        eph = np.exp(self.phi*1j)
        cfacj = eph * np.exp((self.jpoints + 1) * log(1/xi))  # eph/xi**(j+1)
        # Temporary singlet part only!:
        cch = np.einsum('j,sa,sja,ja->j', cfacj,
                        self.gpds.pw_strengths()[:, :2], wce, h[:, :2])
        imh = np.dot(self.wg, cch.imag)
        np.multiply(cch, self.tgj, out=cch)
        reh = np.dot(self.wg, cch.imag)
        return reh, imh

    def _dis_mellin_barnes_integral(self, xi, wce, h):
        """Return Mellin-Barnes integral relevant for DIS."""
        eph = np.exp(self.phi*1j)
        cfacj = eph * np.exp((self.jpoints) * log(1/xi))  # eph/xi**j
        # Temporary singlet part only!:
        cch = np.einsum('j,ja,ja->j', cfacj, wce, h[:, :2])
        mb_int = np.dot(self.wg, cch.imag)
        return mb_int

    def cff(self, pt: g.data.DataPoint) -> np.ndarray:
        """Return array(ReH, ImH, ReE, ...) for kinematic point."""
        if self.nf == 3:
            chargefac = 2./9.
        else:  # nf = 4
            chargefac = 5./18.

        try:
            wce_ar = self.wce[pt.Q2]
        except KeyError:
            # calculate it
            wce_ar = g.evolc.calc_wce(self, pt.Q2, 'DVCS')
            # memorize it for future
            self.wce[pt.Q2] = wce_ar
        h = self.gpds.gpd_H(pt.xi, pt.t)
        reh, imh = self._mellin_barnes_integral(pt.xi, wce_ar, h)
        e = self.gpds.gpd_E(pt.xi, pt.t)
        ree, ime = self._mellin_barnes_integral(pt.xi, wce_ar, e)
        return chargefac * np.array([reh, imh, ree, ime, 0, 0, 0, 0])

    def ImH(self, pt: g.data.DataPoint) -> float:
        """Return Im(CFF H) for kinematic point."""
        return self.cff(pt)[1]

    def ReH(self, pt: g.data.DataPoint) -> float:
        """Return Re(CFF H) for kinematic point."""
        return self.cff(pt)[0]

    def ImE(self, pt: g.data.DataPoint) -> float:
        """Return Im(CFF E) for kinematic point."""
        return self.cff(pt)[3]

    def ReE(self, pt: g.data.DataPoint) -> float:
        """Return Re(CFF E) for kinematic point."""
        return self.cff(pt)[2]

    def tff(self, xi: float, t: float, Q2: float) -> np.ndarray:
        """Return array(ReH_rho, ImH_rho, ReE_rho, ...) of DVrhoP transition FFs."""
        assert self.nf == 4

        astrong = 2 * pi * g.qcd.as2pf(self.p, self.nf,  Q2, self.asp[self.p], self.r20)

        try:
            wce_ar_dvmp = self.wce_dvmp[Q2]
        except KeyError:
            # calculate it
            wce_ar_dvmp = g.evolc.calc_wce(self, Q2, 'DVMP')
            # memorize it for future
            self.wce_dvmp[Q2] = wce_ar_dvmp
        # Evaluations depending on model parameters:
        h_prerot = self.gpds.gpd_H(xi, t)
        # Flavor rotation matrix: (sea,G,uv,dv) --> (SIG, G, NS+, NS-)
        # FIXME: should be calculated only once!
        frot_rho_4 = np.array([[1, 0, 1, 1],
                               [0, 1, 0, 0],
                               [3./20., 0, 5./12., 1./12.],
                               [0, 0, 0, 0]]) / np.sqrt(2)
        h = np.einsum('fa,ja->jf', frot_rho_4, h_prerot)
        reh, imh = self._mellin_barnes_integral(xi, wce_ar_dvmp, h)
        return (g.constants.CF * g.constants.F_rho * astrong / g.constants.NC
                / np.sqrt(Q2) * np.array([reh, imh, 0, 0, 0, 0, 0, 0]))

    def ImHrho(self, pt: g.data.DataPoint) -> np.ndarray:
        """Return Im(TFF H) for kinematic point."""
        tffs = self.tff(pt.xi, pt.t, pt.Q2)
        return tffs[1]

    def ReHrho(self, pt: g.data.DataPoint) -> np.ndarray:
        """Return Re(TFF H) for kinematic point."""
        tffs = self.tff(pt.xi, pt.t, pt.Q2)
        return tffs[0]

    def F2(self, pt: g.data.DataPoint) -> float:
        """Return DIS F2 for kinematic point."""
        if self.nf == 3:
            chargefac = 2./9.
        else:  # nf = 4
            chargefac = 5./18.

        try:
            wce_ar_dis = self.wce_dis[pt.Q2]
        except KeyError:
            # calculate it, first PW is the only relevant one
            wce_ar_dis = g.evolc.calc_wce(self, pt.Q2, 'DIS')[0, :, :]
            # memorize it for future
            self.wce_dis[pt.Q2] = wce_ar_dis
        pdf = self.gpds.gpd_H(0, 0)  # forward limit
        # print('pdf = {}'.format(pdf[0, :2]))
        # print('wce = {}'.format(wce_ar_dis[0, :2]))
        mb_int = self._dis_mellin_barnes_integral(pt.xB, wce_ar_dis, pdf)
        return chargefac * mb_int / np.pi


class ComptonDispersionRelations(ComptonFormFactors):
    """Use dispersion relations for real parts of CFFs.

    methods: ReH, ReE, ReHt, ReEt, subtraction
    Subclass should implement ansaetze for ImH, ImE, ImHt, ImEt
    and subtraction. This class implements just dispersion integrals.
    """

    def __init__(self) -> None:
        """Init."""
        super().__init__()

    def dispargV(self, x, fun, pt):
        """Integrand of the dispersion integral (vector case).

        fun -- Im(CFF)
        With variable change x->x^(1/(1-ga))=u in
        order to tame the singularity at x=0.

        """
        ga = 0.9  # Nice value obtained by experimentation in Mathematica
        u = x**(1./(1.-ga))
        res = u**ga * (fun(pt, u) - fun(pt))
        return (2.*u) / (pt.xi**2 - u**2) * res / (1.-ga)

    def dispargA(self, x, fun, pt):
        """Integrand of the dispersion integral (axial-vector case).

        fun -- Im(CFF)
        With variable change x->x^(1/(1-ga))=u in
        order to tame the singularity at x=0.

        """
        ga = 0.9  # Value same as for V-case (TODO: is this the best choice?)
        u = x**(1./(1.-ga))
        res = u**ga * (fun(pt, u) - fun(pt))
        return (2. * pt.xi) / (pt.xi**2 - u**2) * res / (1.-ga)

    def subtraction(self, pt):
        """Subtraction constant."""
        return 0  # default

    def ReH(self, pt):
        """Real part of CFF H.

        Given by dispersion integral over ImH minus subtraction constant.

        """
        res = g.quadrature.PVquadrature(self.dispargV, 0, 1, (self.ImH, pt))
        pvpi = (res + log(pt.xi**2 / (1.-pt.xi**2)) * self.ImH(pt)) / pi
        # P.V./pi - subtraction constant C/(1-t/MC^2)^2
        return pvpi - self.subtraction(pt)

    def ReHt(self, pt):
        """Real part of CFF Ht.

        Given by dispersion integral over ImHt.

        """
        res = g.quadrature.PVquadrature(self.dispargA, 0, 1, (self.ImHt, pt))
        pvpi = (res + log((1.+pt.xi)/(1.-pt.xi)) * self.ImHt(pt))/pi
        return pvpi   # this is P.V./pi

    def ReE(self, pt):
        """Real part of CFF E.

        Given by dispersion integral over ImE plus subtraction constant.

        """
        res = g.quadrature.PVquadrature(self.dispargV, 0, 1, (self.ImE, pt))
        pvpi = (res + log(pt.xi**2 / (1.-pt.xi**2)) * self.ImE(pt)) / pi
        # This is same subtraction constant
        # as for H, but with opposite sign
        return pvpi + self.subtraction(pt)

    def ReEt(self, pt):
        """Real part of CFF Et.

        Given by dispersion integral over ImEt

        """
        res = g.quadrature.PVquadrature(self.dispargA, 0, 1, (self.ImEt, pt))
        pvpi = (res + log((1.+pt.xi)/(1.-pt.xi)) * self.ImEt(pt))/pi
        return pvpi   # this is P.V./pi


class PionPole(object):
    """Various options for pion-pole contribution."""

    def DMfixpole(self, pt):
        """Fixed pion-pole as used by Dieter."""
        pole = (2.2390424 * (1. - (1.7*(0.0196 - pt.t))/(1.
                - pt.t/2.)**2))/((0.0196 - pt.t)*pt.xi)
        if 'in2particle' in pt and pt.in2particle == 'n':
            return -pole  # neutron
        else:
            return pole  # proton

    def DMfreepole(self, pt):
        """Free pion-pole as proposed by Dieter."""
        pole = (self.parameters['rpi'] * 2.16444 / (0.0196 - pt.t) / (1.
                - pt.t/self.parameters['Mpi']**2)**2 / pt.xi)
        if 'in2particle' in pt and pt.in2particle == 'n':
            return -pole  # neutron
        else:
            return pole  # proton


class ComptonModelDR(ComptonDispersionRelations, PionPole):
    """Model for CFFs as in arXiv:0904.0458."""

    def __init__(self, **kwargs):
        """Constructor."""
        # initial values of parameters and limits on their values
        self.parameters = {'Nsea': 1.5, 'alS': 1.13, 'alpS': 0.15,
                           'MS': 0.707, 'rS': 1.0,
                           'bS': 2.0,      'limit_bS': (0.4, 5.0),
                           'Nv': 1.35,
                           'alv': 0.43,
                           'alpv': 0.85,
                           'Mv': 1.0,     'limit_Mv': (0.4, 1.5),
                           'rv': 0.5,     'limit_rv': (0., 8.),
                           'bv': 2.2,     'limit_bv': (0.4, 5.),
                           'C': 7.0,      'limit_C': (-10., 10.),
                           'MC': 1.3,     'limit_MC': (0.4, 2.),
                           'tNv': 0.0,
                           'tal': 0.43,
                           'talp': 0.85,
                           'tMv': 2.7,    'limit_tMv': (0.4, 2.),
                           'trv': 6.0,    'limit_trv': (0., 8.),
                           'tbv': 3.0,    'limit_tbv': (0.4, 5.)}

        # order matters to fit.MinuitFitter, so it is defined by:
#         self.parameter_names = ['Nsea', 'alS', 'alpS', 'MS', 'rS', 'bS',
#                                 'Nv', 'alv', 'alpv', 'Mv', 'rv', 'bv',
#                                 'C', 'MC',
#                                 'tNv', 'tal', 'talp',
#                                 'tMv', 'trv', 'tbv']
        super().__init__()

    def subtraction(self, pt):
        return self.parameters['C']/(1.-pt.t/self.parameters['MC']**2)**2

    def ImH(self, pt, xi=0):
        """Imaginary part of CFF H."""
        p = self.parameters  # just a shortcut
        # FIXME: The following solution is not elegant
        if isinstance(xi, np.ndarray):
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
        if 'in2particle' in pt and pt.in2particle == 'n':
            chgfac = (1.*4./9. + 2./9.)  # neutron
        else:
            chgfac = (2.*4./9. + 1./9.)  # proton
        val = (chgfac * p['Nv'] * p['rv'] * twox**(-p['alv']-p['alpv']*t) *
               onex**p['bv'] / (1. - onex*t/(p['Mv']**2)))
        sea = ((2./9.) * p['Nsea'] * p['rS'] * twox**(-p['alS']-p['alpS']*t) *
               onex**p['bS'] / (1. - onex*t/(p['MS']**2))**2)
        return pi * (val + sea) / (1.+x)

    def ImHt(self, pt, xi=0):
        """Imaginary part of CFF Ht i.e. tilde{H}."""
        p = self.parameters  # just a shortcut
        # FIXME: The following solution is not elegant
        if isinstance(xi, np.ndarray):
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
        if 'in2particle' in pt and pt.in2particle == 'n':
            chgfac = (1.*4./9. + 2./9.)  # neutron
        else:
            chgfac = (2.*4./9. + 1./9.)  # proton
        val = (chgfac * p['tNv'] * p['trv']
               * twox**regge * onex**p['tbv'] / (1. - onex*t/(p['tMv']**2)))
        return pi * val / (1.+x)

    def ImE(self, pt, xi=0):
        """Imaginary part of CFF E."""
        # Just changing function signature w.r.t. ComptonFormFactors
        # to make it compatible for dispersion integral
        return 0

    def ReEt(self, pt):
        """Instead of disp. rel. use pole formula."""
        return self.DMfixpole(pt)

# For compatibility with old models in database:
# ComptonModelDRsea = ComptonModelDR


class ComptonModelDRPP(ComptonModelDR):
    """Model for CFFs as in arXiv:0904.0458. + free pion pole."""

    def __init__(self, **kwargs):
        """Constructor."""
        # First inhert what's needed
        ComptonModelDR.__init__(self, **kwargs)
        # Adding two extra parameters:
        self.parameters.update({
             'rpi': 1.0,    'limit_rpi': (-8, 8.),
             'Mpi': 1.0,    'limit_Mpi': (0.4, 4.)})
        # self.parameter_names.append('rpi')
        # self.parameter_names.append('Mpi')
        # now do whatever else is necessary
        # ComptonFormFactors.__init__(self, **kwargs)

    def ReEt(self, pt):
        """Instead of disp. rel. use pole formula."""
        return self.DMfreepole(pt)


class ComptonHybrid(ComptonFormFactors):
    """This combines MB model for small xB and DR model for valence xB."""

    def __init__(self, instMB: MellinBarnesModel, instDR: ComptonModelDR, **kwargs):
        """Initializes with one instance of MB model and one of DR model."""
        self.MB = instMB
        self.DR = instDR
        self.DR.parameters['Nsea'] = 0.  # sea comes from Gepard part
        self.parameters = {**self.MB.parameters, **self.DR.parameters}
        # Consolidate parameters of all models, base and derived
        self.MB.parameters = self.parameters
        self.DR.parameters = self.parameters
        self.MB.gpds.parameters = self.parameters
        # self.parameter_names = self.DR.parameter_names + self.Gepard.parameter_names

    def is_within_model_kinematics(self, pt):
        """Is kinematics of datapoint ok?"""
        # relaxing xBmin and removing Q2max
        return ((1.5 <= pt.Q2) and
                (pt.tm < min(1., pt.Q2/4)) and
                (1e-5 < pt.xB < 0.65))

    # FIXME: this below looks inconsistent generally for xi != pt.xi !!

    def ImH(self, pt, xi=0):
        return self.MB.ImH(pt) + self.DR.ImH(pt, xi)

    def ReH(self, pt):
        return self.MB.ReH(pt) + self.DR.ReH(pt)

    def ImE(self, pt, xi=0):
        return self.MB.ImE(pt) + self.DR.ImE(pt, xi)

    def ReE(self, pt):
        return self.MB.ReE(pt) + self.DR.ReE(pt)

    # Tildes are not provided by MB model

    def ImHt(self, pt, xi=0):
        return self.DR.ImHt(pt, xi)

    def ReHt(self, pt):
        return self.DR.ReHt(pt)

    def ImEt(self, pt):
        return self.DR.ImEt(pt)

    def ReEt(self, pt):
        return self.DR.ReEt(pt)


#  --- Complete models ---
class ModelDR(ComptonModelDR, ElasticDipole):
    """Complete model as in arXiv:0904.0458.."""


class ModelDRKelly(ComptonModelDR, ElasticKelly):
    """Same, but with Kelly elastic form factors."""


class HybridDipole(ComptonHybrid, ElasticDipole):
    """Complete hybrid model."""


class HybridKelly(ComptonHybrid, ElasticKelly):
    """Complete hybrid model."""

