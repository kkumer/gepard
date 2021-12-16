"""Definitions of Compton Form Factors.

"""
from __future__ import annotations

from math import log, pi, sqrt
from typing import Dict, List

import numpy as np
import scipy.stats

from . import constants, data, eff, gpd, mellin, model, qcd, quadrature, utils, wilson

# from joblib import Parallel, delayed

class ComptonFormFactors(model.ParameterModel):
    """Base class for CFFs.

    CFFs are set to be zero here. Actual models are built by subclassing this.
    """

    allCFFs = ['ImH', 'ReH', 'ImE', 'ReE', 'ImHt', 'ReHt', 'ImEt', 'ReEt']
    allCFFsb = ['ImH', 'ReH', 'ImE', 'ReE', 'ImHt', 'ReHt', 'ImEt', 'ReEb']
    allCFFeffs = ['ImHeff', 'ReHeff', 'ImEeff', 'ReEeff',
                  'ImHteff', 'ReHteff', 'ImEteff', 'ReEteff']
    # allGPDs = []

    def __init__(self, **kwargs):
        # squared DVCS charge factors
        if self.nf == 3:
            qs = 2/9
            qns = 1/9
        else:  # nf = 4
            qs = 5/18
            qns = 1/6
        self.dvcs_charges = (qs, qs, qns)
        super().__init__(**kwargs)

    def print_CFFs(self, pt, format=None):
        """Print values of CFFs at given kinematic point."""
        vals = [getattr(self, cff)(pt) for cff in self.allCFFs]
        if format == 'mma':
            s = "{" + 8*"%s -> %f, "
            s = s[:-2] + "}"
        else:
            s = 8*"%4s = %5.2f\n"
        print(s % utils.flatten(tuple(zip(self.allCFFs, vals))))

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


class MellinBarnesCFF(ComptonFormFactors, mellin.MellinBarnes):
    """Class of models built by Mellin-Barnes integration."""

    def __init__(self, **kwargs) -> None:
        """Init MellinBarnes class and pre-calculate stuff.

        """
        #
        self.tgj = np.tan(pi*self.jpoints/2.)
        # wce[Q2] = wce[spw, j, a] - Wilson coeffs evolved; local to model instance
        self.wce: Dict[float, np.ndarray] = {}  # DVCS Wilson coefs.
        mellin.MellinBarnes.__init__(self, **kwargs)
        super().__init__(**kwargs)

    def cff(self, pt: data.DataPoint) -> np.ndarray:
        """Return array(ReH, ImH, ReE, ...) for kinematic point."""
        try:
            wce_ar = self.wce[pt.Q2]
        except KeyError:
            # calculate it
            wce_ar = wilson.calc_wce(self, pt.Q2, 'DVCS')
            # memorize it for future
            self.wce[pt.Q2] = wce_ar
        # Evaluations depending on model parameters:
        h_prerot = self.gpd_H(pt.xi, pt.t)
        h = np.einsum('f,fa,ja->jf', self.dvcs_charges, self.frot, h_prerot)
        reh, imh = self._mellin_barnes_integral(pt.xi, wce_ar, h)
        e_prerot = self.gpd_E(pt.xi, pt.t)
        e = np.einsum('f,fa,ja->jf', self.dvcs_charges, self.frot, e_prerot)
        ree, ime = self._mellin_barnes_integral(pt.xi, wce_ar, e)
        return np.array([reh, imh, ree, ime, 0, 0, 0, 0])

    def ReH(self, pt: data.DataPoint) -> float:
        """Return Re(CFF H) for kinematic point."""
        return self.cff(pt)[0]

    def ImH(self, pt: data.DataPoint) -> float:
        """Return Im(CFF H) for kinematic point."""
        return self.cff(pt)[1]

    def ReE(self, pt: data.DataPoint) -> float:
        """Return Re(CFF E) for kinematic point."""
        return self.cff(pt)[2]

    def ImE(self, pt: data.DataPoint) -> float:
        """Return Im(CFF E) for kinematic point."""
        return self.cff(pt)[3]

class ComptonDispersionRelations(ComptonFormFactors):
    """Use dispersion relations for real parts of CFFs.

    methods: ReH, ReE, ReHt, ReEt, subtraction
    Subclass should implement ansaetze for ImH, ImE, ImHt, ImEt
    and subtraction. This class implements just dispersion integrals.
    """

    def __init__(self, **kwargs) -> None:
        """Init."""
        super().__init__(**kwargs)

    def dispargV(self, x, fun, pt):
        """Integrand of the dispersion integral (vector case).

        fun -- Im(CFF)
        With variable change x->x^(1/(1-ga))=u in
        order to tame the singularity at x=0.

        """
        ga = 0.9  # Nice value obtained by experimentation in Mathematica
        u = x**(1./(1.-ga))
        if hasattr(fun, '__self__'):
            res = u**ga * (fun(pt, u) - fun(pt))
        else:  # unbound method
            res = u**ga * (fun(self, pt, u) - fun(self, pt))
        return (2.*u) / (pt.xi**2 - u**2) * res / (1.-ga)

    def dispargA(self, x, fun, pt):
        """Integrand of the dispersion integral (axial-vector case).

        fun -- Im(CFF)
        With variable change x->x^(1/(1-ga))=u in
        order to tame the singularity at x=0.

        """
        ga = 0.9  # Value same as for V-case (TODO: is this the best choice?)
        u = x**(1./(1.-ga))
        if hasattr(fun, '__self__'):
            res = u**ga * (fun(pt, u) - fun(pt))
        else:  # unbound method
            res = u**ga * (fun(self, pt, u) - fun(self, pt))
        return (2. * pt.xi) / (pt.xi**2 - u**2) * res / (1.-ga)

    def subtraction(self, pt):
        """Subtraction constant."""
        return 0  # default

    def _ReV(self, pt, imfun, subsign):
        """Real part of vector CFFs H and E.

        Given by dispersion integral over imfun minus/plus subtraction constant.

        """
        res = quadrature.PVquadrature(self.dispargV, 0, 1, (imfun, pt))
        if hasattr(imfun, '__self__'):
            pvpi = (res + log(pt.xi**2 / (1.-pt.xi**2)) * imfun(pt)) / pi
        else:
            pvpi = (res + log(pt.xi**2 / (1.-pt.xi**2)) * imfun(self, pt)) / pi
        # P.V./pi - subtraction constant C/(1-t/MC^2)^2
        return pvpi + subsign * self.subtraction(pt)

    def _ReA(self, pt, imfun):
        """Real part of axial vector CFFs Ht and Et.

        Given by dispersion integral over imfun.

        """
        res = quadrature.PVquadrature(self.dispargA, 0, 1, (imfun, pt))
        if hasattr(imfun, '__self__'):
            pvpi = (res + log((1.+pt.xi)/(1.-pt.xi)) * imfun(pt))/pi
        else:  # unbound method
            pvpi = (res + log((1.+pt.xi)/(1.-pt.xi)) * imfun(self, pt))/pi
        return pvpi   # this is P.V./pi

    def ReH(self, pt, imfun=None):
        """Real part of CFF H."""
        if imfun is None:
            imfun = self.ImH
        return self._ReV(pt, imfun, -1)

    def ReE(self, pt, imfun=None):
        """Real part of CFF E."""
        if imfun is None:
            imfun = self.ImE
        return self._ReV(pt, imfun, +1)

    def ReHt(self, pt, imfun=None):
        """Real part of CFF Ht."""
        if imfun is None:
            imfun = self.ImHt
        return self._ReA(pt, imfun)

    def ReEt(self, pt, imfun=None):
        """Real part of CFF Et."""
        if imfun is None:
            imfun = self.ImEt
        return self._ReA(pt, imfun)


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
        """Init ComptonModelDR object.

        Args:
            nf: number of active quark flavors

        """
        self.nf = kwargs.setdefault('nf', 4)
        # initial values of parameters and limits on their values
        self.parameters = {'Nsea': 1.5, 'alS': 1.13, 'alpS': 0.15,
                           'MS': 0.707, 'rS': 1.0, 'bS': 2.0,
                           'Nv': 1.35, 'alv': 0.43, 'alpv': 0.85,
                           'Mv': 1.0, 'rv': 0.5, 'bv': 2.2,
                           'C': 7.0, 'MC': 1.3,
                           'tNv': 0.0, 'tal': 0.43, 'talp': 0.85,
                           'tMv': 2.7, 'trv': 6.0, 'tbv': 3.0}

        self.parameters_limits = {'bS': (0.4, 5.0),
                                  'Mv': (0.4, 1.5), 'rv': (0., 8.), 'bv': (0.4, 5.),
                                  'C': (-10., 10.), 'MC': (0.4, 2.),
                                  'tMv': (0.4, 2.), 'trv': (0., 8.), 'tbv': (0.4, 5.)}

        # order matters to fit.MinuitFitter, so it is defined by:
#         self.parameter_names = ['Nsea', 'alS', 'alpS', 'MS', 'rS', 'bS',
#                                 'Nv', 'alv', 'alpv', 'Mv', 'rv', 'bv',
#                                 'C', 'MC',
#                                 'tNv', 'tal', 'talp',
#                                 'tMv', 'trv', 'tbv']
        super().__init__(**kwargs)

    def subtraction(self, pt):
        """Dispersion relations subtraction constant."""
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
        # Adding two extra parameters:
        self.parameters.update({'rpi': 1.0,  'Mpi': 1.0})

        self.parameters_limits.update({
             'rpi': (-8, 8.),
             'Mpi': (0.4, 4.)})
        super().__init__(**kwargs)

    def ReEt(self, pt):
        """Instead of disp. rel. use pole formula."""
        return self.DMfreepole(pt)


class ComptonHybrid(MellinBarnesCFF, ComptonModelDR, ComptonFormFactors):
    """This combines MB model for small xB and DR model for valence xB."""


    def is_within_model_kinematics(self, pt):
        """Is kinematics of datapoint ok?"""
        # relaxing xBmin and removing Q2max
        return ((1.5 <= pt.Q2) and
                (pt.tm < min(1., pt.Q2/4)) and
                (1e-5 < pt.xB < 0.65))

    # FIXME: this below looks inconsistent generally for xi != pt.xi !!

    def ImH(self, pt, xi=0):
        return MellinBarnesCFF.ImH(self, pt) + ComptonModelDR.ImH(self, pt, xi)

    def ReH(self, pt):
        return MellinBarnesCFF.ReH(self, pt) + ComptonModelDR.ReH(self, pt, imfun=ComptonModelDR.ImH)

    def ImE(self, pt, xi=0):
        return MellinBarnesCFF.ImE(self, pt) + ComptonModelDR.ImE(self, pt, xi)

    def ReE(self, pt):
        return MellinBarnesCFF.ReE(self, pt) + ComptonModelDR.ReE(self, pt, imfun=ComptonModelDR.ImE)

    # Tildes are not provided by MB model

    def ImHt(self, pt, xi=0):
        return ComptonModelDR.ImHt(self, pt, xi)

    def ReHt(self, pt):
        return ComptonModelDR.ReHt(self, pt, imfun=ComptonModelDR.ImHt)

    def ImEt(self, pt):
        return ComptonModelDR.ImEt(self, pt)

    def ReEt(self, pt):
        """Instead of disp. rel. use fixed pion pole formula."""
        return self.DMfixpole(pt)


class ComptonHybridPP(ComptonHybrid):
    """This combines MB model for small xB and DRPP model for valence xB."""


    def is_within_model_kinematics(self, pt):
        """Is kinematics of datapoint ok?"""
        # relaxing xBmin and removing Q2max
        return ((1.5 <= pt.Q2) and
                (pt.tm < min(1., pt.Q2/4)) and
                (1e-5 < pt.xB < 0.65))

    def ReEt(self, pt):
        """Instead of disp. rel. use free pion pole formula."""
        return self.DMfreepole(pt)


#  --- Complete models ---
class ModelDR(ComptonModelDR, eff.ElasticDipole):
    """Complete model as in arXiv:0904.0458.."""


class ModelDRKelly(ComptonModelDR, eff.ElasticKelly):
    """Same, but with Kelly elastic form factors."""


class HybridDipole(ComptonHybrid, eff.ElasticDipole):
    """Complete hybrid model."""


class HybridKelly(ComptonHybrid, eff.ElasticKelly):
    """Complete hybrid model."""
