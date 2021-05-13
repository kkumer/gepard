"""Definitions of models.

Todo:
    * Refactor and transfer everything from old Model.py
"""

from cmath import exp, log, pi
from typing import Dict, List

import numpy as np
from joblib import Parallel, delayed

import gepard as g


class Model(object):
    """Base class for all models.

    Instance of Model Python class specifies structure of relevant
    hadrons so that observables can be calculated.  Methods provided are
    typically GPDs, CFFs, elastic FFs, DVMP transition TFFs, DAs, etc.
    Main subclasses are:

     - ParameterModel which depends on parameters (some of which can be
       provided by minimization routine in the fitting procedure).
     - NeuralModel where structure functions are represented as neural nets
     - NumericModel where values of structure functions are represented as
       grids of numbers, which may be interpolated
    """
    def __init__(self, name: str = None, texname: str = None,
                 description: str = None) -> None:
        """Init Model class.

        Args:
            name: short model name (used as model database)
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


class ConformalSpaceGPD(ParameterModel):
    """Base class of GPD models built in conformal moment space."""

    def __init__(self, p: int = 0, scheme: str = 'MSBAR', nf: int = 4,
                 q02: float = 4.0) -> None:
        """Init and pre-calculate stuff.

        Args:
            p: pQCD order (0 = LO, 1 = NLO, 2 = NNLO)
            scheme: pQCD scheme  (MSBAR or CSBAR)
            nf: number of active quark flavors
            q02: Initial Q0^2 for pQCD evolution.

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
        # alpha_strong/(2*pi) at scale r0^2
        self.asp = np.array([0.0606, 0.0518, 0.0488])
        self.r20 = 2.5
        npoints, weights = g.quadrature.mellin_barnes()
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
                           'ng': 0.4,
                           'al0g': 1.2,
                           'alpg': 0.25,
                           'mg': 1.2,
                           'secg': 0.,
                           'thig': 0.}
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

        # gpd_H: np.ndarray    # trying to appease mypy


class Test(ConformalSpaceGPD):
    """Simple testing ansatz for GPDs."""

    def __init__(self, **kwargs) -> None:
        """See parent `ConformalSpaceGPD` class for docs."""
        kwargs.setdefault('scheme', 'CSBAR')
        kwargs.setdefault('nf', 3)
        kwargs.setdefault('q02', 1.0)
        super().__init__(**kwargs)
        self.nf = 3
        self.q02 = 1.0
        self.asp = np.array([0.05, 0.05, 0.05])
        self.r20 = 2.5

    def gpd_H(self, eta: float, t: float) -> np.ndarray:
        """Return (npts, 4) array H_j^a for all j-points and 4 flavors."""
        h = []
        for j in self.jpoints:
            h.append(g.gpdj.test(j, t, self.parameters))
        return np.array(h)


class Fit(ConformalSpaceGPD):
    """Singlet fitting ansatz for GPDs with three SO(3) partial waves."""

    def __init__(self, **kwargs) -> None:
        """See parent `ConformalSpaceGPD` class for docs."""
        kwargs.setdefault('scheme', 'MSBAR')
        kwargs.setdefault('nf', 4)
        kwargs.setdefault('q02', 4.0)
        super().__init__(**kwargs)

    def gpd_H_single(self, eta: float, t: float) -> np.ndarray:
        """Return (npts, 4) array H_j^a for all j-points and 4 flavors."""
        h = []
        for j in self.jpoints:
            h.append(g.gpdj.fit(j, t, self.parameters))
        return np.array(h)

    def gpd_H_para(self, eta: float, t: float) -> np.ndarray:
        """Return (npts, 4) array H_j^a for all j-points and 4 flavors."""
        h = Parallel(n_jobs=20)(delayed(g.gpdj.fit)(j, t, self.parameters)
                                for j in self.jpoints)
        return np.array(h)

    gpd_H = gpd_H_single  # multiprocessing version is actually slower


class MellinBarnesModel(ParameterModel):
    """Class of models built by Mellin-Barnes integration.

    Todo:
        It should most likely NOT inherit from ConformalSpaceModel!
    """

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
        self.wg = gpds.wg
        self.gpds = gpds
        self.parameters = gpds.parameters
        self.tgj = np.tan(pi*self.jpoints/2.)
        # wce[q2] = wce[spw, j, a] - Wilson coeffs evolved; local to model instance
        self.wce: Dict[float, np.ndarray] = {}  # DVCS
        self.wce_dvmp: Dict[float, np.ndarray] = {}  # DVMP
        super().__init__()

    def _mellin_barnes_integral(self, xi, wce, h):
        """Return Mellin-Barnes integral.

        Integrates convolution of evolved GPDs with Wilson coefs.

        Todo:
            Only GPD H treated atm.
        """
        phij = 1.57079632j
        eph = exp(phij)
        cfacj = eph * np.exp((self.jpoints + 1) * log(1/xi))  # eph/xi**(j+1)
        # Temporary singlet part only!:
        cch = np.einsum('j,sa,sja,ja->j', cfacj,
                        self.gpds.pw_strengths()[:, :2], wce, h[:, :2])
        imh = np.dot(self.wg, cch.imag)
        np.multiply(cch, self.tgj, out=cch)
        reh = np.dot(self.wg, cch.imag)
        return reh, imh

    def cff(self, xi: float, t: float, q2: float) -> np.ndarray:
        """Return array(ReH, ImH, ReE, ...) for kinematic point."""
        if self.nf == 3:
            chargefac = 2./9.
        else:  # nf = 4
            chargefac = 5./18.

        try:
            wce_ar = self.wce[q2]
        except KeyError:
            # calculate it
            wce_ar = g.evolc.calc_wce(self, q2)
            # memorize it for future
            self.wce[q2] = wce_ar
        h = self.gpds.gpd_H(xi, t)
        reh, imh = self._mellin_barnes_integral(xi, wce_ar, h)
        return chargefac * np.array([reh, imh, 0, 0, 0, 0, 0, 0])

    def tff(self, xi: float, t: float, q2: float) -> np.ndarray:
        """Return array(ReH_rho, ImH_rho, ReE_rho, ...) of DVrhoP transition FFs."""
        assert self.nf == 4

        astrong = 2 * pi * g.qcd.as2pf(self.p, self.nf,  q2, self.asp[self.p], self.r20)

        try:
            wce_ar_dvmp = self.wce_dvmp[q2]
        except KeyError:
            # calculate it
            wce_ar_dvmp = g.evolc.calc_wce_dvmp(self, q2)
            # memorize it for future
            self.wce_dvmp[q2] = wce_ar_dvmp
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
                / np.sqrt(q2) * np.array([reh, imh, 0, 0, 0, 0, 0, 0]))
