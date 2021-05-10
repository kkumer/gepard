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

        Warning:
            Presently, only one instance of this class can be reliably used
            in a single Python session, since it relies on Fortran extension
            module common block global state!

        Todo:
            - Make different Q2's possible
            - Move at least evolved Wilson coeffs from Fortran.

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
        """Return (4, npts) array H^a_j for 4 flavors and all j-points."""
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
        """Return (4, npts) array H^a_j for 4 flavors and all j-points."""
        h = []
        for j in self.jpoints:
            h.append(g.gpdj.fit(j, t, self.parameters))
        return np.array(h)

    def gpd_H_para(self, eta: float, t: float) -> np.ndarray:
        """Return (4, npts) array H^a_j for 4 flavors and all j-points."""
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

        Todo:
            This now simply imports everything from old Fortran
            extension. We should consider Python-only implementation.
            (Wilson coeffs and anomalous dimensions can stay Fortran
            for some longer time.)

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
        self.wce: Dict[float, np.ndarray] = {}
        super().__init__()

    def cff(self, xi: float, t: float, q2: float) -> np.ndarray:
        """Return array(ReH, ImH, ReE, ...) for kinematic point."""
        if self.nf == 3:
            chargefac = 2./9.
        else:  # nf = 4
            chargefac = 5./18.

        # Evaluations depending on model parameters:
        h = self.gpds.gpd_H(xi, t)
        pw_strengths = np.array([[1., 1., 1, 1],
                                 [self.parameters['secs'],
                                     self.parameters['secg'], 0, 0],
                                 [self.parameters['this'],
                                     self.parameters['thig'], 0, 0]])
        try:
            wce_ar = self.wce[q2]
        except KeyError:
            # calculate it
            wce_ar = g.evolc.calc_wce(self, q2)
            # memorize it for future
            self.wce[q2] = wce_ar
        phij = 1.57079632j
        eph = exp(phij)
        cfacj = eph * np.exp((self.jpoints + 1) * log(1/xi))  # eph/xi**(j+1)
        # print('pw_strengths[1, 0] = {}'.format(pw_strengths[1, 0]))
        # if t < -0.9:
        #     print('t, q2 = {}, {}'.format(t, q2))
        #     print('wce[0, 0, 0] = {}'.format(wce[0, 0, 0]))
        #     print('qind, qs = {} -> {}'.format(qind, self.qs[5, :4]))
        #     print('id(wce) = {}'.format(id(wce)))
        # print('h[0, 0] = {}'.format(h[0, 0]))
        # cch = np.einsum('j,sa,sja,ja->j', cfacj, pw_strengths, wce, h)
        # Temporary SEC=0 WCE with singlet part only!:
        # cch = np.einsum('j,ja,ja->j', cfacj, wce_ar, h[:, :2])
        # Temporary singlet part only!:
        cch = np.einsum('j,sa,sja,ja->j', cfacj, pw_strengths[:, :2], wce_ar, h[:, :2])
        imh = chargefac * np.dot(self.wg, cch.imag)
        np.multiply(cch, self.tgj, out=cch)
        reh = chargefac * np.dot(self.wg, cch.imag)
        # FIXME: Only CFF H at the moment.
        return np.array([reh, imh, 0, 0, 0, 0, 0, 0])
