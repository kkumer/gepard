"""Definitions of Compton Form Factors."""

from __future__ import annotations

from math import log, pi
from typing import Dict, Union

import numpy as np

from . import data, gpd, model, quadrature, wilson  # noqa: F401


def _flatten(T):
    """Flatten the tuple."""
    if not isinstance(T, tuple):
        return (T,)
    elif len(T) == 0:
        return ()
    else:
        return _flatten(T[0]) + _flatten(T[1:])


class CFF(model.ParameterModel):
    """Base class for CFFs.

    CFFs are set to be zero here. Actual models are built by subclassing this.
    """

    allCFFs = ['ImH', 'ReH', 'ImE', 'ReE', 'ImHt', 'ReHt', 'ImEt', 'ReEt']
    allCFFsb = ['ImH', 'ReH', 'ImE', 'ReE', 'ImHt', 'ReHt', 'ImEt', 'ReEb']
    allCFFeffs = ['ImHeff', 'ReHeff', 'ImEeff', 'ReEeff',
                  'ImHteff', 'ReHteff', 'ImEteff', 'ReEteff']
    # allGPDs = []

    def __init__(self, **kwargs):
        self.nf = kwargs.setdefault('nf', 4)
        # squared DVCS charge factors
        if self.nf == 3:
            qs = 2/9
            qns = 1/9
        else:  # nf = 4
            qs = 5/18
            qns = 1/6
        self.dvcs_charges = (qs, qs, qns)
        super().__init__(**kwargs)

    def print_CFFs(self, pt: data.DataPoint, format: str = None):
        """Print values of CFFs at given kinematic point.

        Args:
            pt: DataPoint instance
            format: 'mma' to get output pastable into Mathematica

        """
        vals = [getattr(self, cff)(pt) for cff in self.allCFFs]
        if format == 'mma':
            s = "{" + 8*"%s -> %f, "
            s = s[:-2] + "}"
        else:
            s = 8*"%4s = %5.2f\n"
        print(s % _flatten(tuple(zip(self.allCFFs, vals))))

    # Initial definition of all CFFs. All just return zero.
    for name in allCFFs:
        exec('def %s(self, pt): return 0.' % name)

    for name in allCFFeffs:
        exec('def %s(self, pt): return 0.' % name)

    def ReEb(self, pt: data.DataPoint):
        """Return Re(E-bar).

        Args:
            pt: DataPoint instance

        Returns:
            xi*ReEt

        """
        return (pt.xi)*self.ReEt(pt)

    def cff(self, pt: data.DataPoint) -> np.ndarray:
        """Return array(ReH, ImH, ReE, ...) for kinematic point.

        Args:
            pt: DataPoint instance

        Returns:
            Eight CFFs (ReH, ImH, ReE, ..., ImEt).

        This is generic top-class method that only gathers results
        of dedicated ReH(), ImH(), etc. functions.
        The idea is that this method can be implemented more efficiently
        in subclases where components of CFFs can be calculated at the
        same time or in parallel.

        """
        return np.array([self.ReH(pt), self.ImH(pt), self.ReE(pt), self.ImE(pt),
                         self.ReHt(pt), self.ImHt(pt), self.ReEt(pt), self.ImEt(pt)])

    def is_within_model_kinematics(self, pt: data.DataPoint):
        """Check if kinematics is sensible.

        Args:
            pt: DataPoint instance

        Returns:
            True if pt kinematics is inside the model validity range.

        """
        return ((1.5 <= pt.Q2 <= 5.) and
                (pt.tm < min(1., pt.Q2/4)) and
                (1e-3 < pt.xB < 0.5))


class MellinBarnesCFF(CFF):
    """Class of models built by Mellin-Barnes integration."""

    def __init__(self, **kwargs) -> None:
        self.tgj = np.tan(pi*self.jpoints/2.)
        # wce[Q2] = wce[spw, j, a] - Wilson coeffs evolved; local to model instance
        self.wce: Dict[float, np.ndarray] = {}  # DVCS Wilson coefs.
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
        h_prerot = self.H(pt.xi, pt.t)
        h = np.einsum('f,fa,ja->jf', self.dvcs_charges, self.frot, h_prerot)
        e_prerot = self.E(pt.xi, pt.t)
        e = np.einsum('f,fa,ja->jf', self.dvcs_charges, self.frot, e_prerot)
        reh, imh, ree, ime = self._mellin_barnes_integral_HE(pt.xi, wce_ar, h, e)
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


class DispersionCFF(CFF):
    """Use dispersion relations for real parts of CFFs.

    methods: ReH, ReE, ReHt, ReEt, subtraction
    Subclass should implement ansaetze for ImH, ImE, ImHt, ImEt
    and subtraction. This class implements just dispersion integrals.

    """
    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)

    def dispargV(self, x, fun, pt):
        """Integrand of the dispersion integral (vector case).

        Args:
            x: integration variable, longitudinal momentum fraction
            fun: function providing imaginary part of CFF
            pt: DataPoint instance

        Returns:
            Integrand of dispersion integral.
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

        Args:
            x: integration variable, longitudinal momentum fraction
            fun: function providing imaginary part of CFF
            pt: DataPoint instance

        Returns:
            Integrand of dispersion integral.
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

        Args:
            pt: DataPoint instance
            imfun: function providing imaginary part of CFF
            subsign: sign of subtraction constant (-1 for H, +1 for E)

        Returns:
            Dispersion integral over imfun minus/plus subtraction constant.

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

        Args:
            pt: DataPoint instance
            imfun: function providing imaginary part of CFF

        Returns:
            Dispersion integral over imfun.

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
        """Return fixed pion-pole as used by Dieter."""
        pole = (2.2390424 * (1. - (1.7*(0.0196 - pt.t))/(1.
                - pt.t/2.)**2))/((0.0196 - pt.t)*pt.xi)
        if 'in2particle' in pt and pt.in2particle == 'n':
            return -pole  # neutron
        else:
            return pole  # proton

    def DMfreepole(self, pt):
        """Return free pion-pole as proposed by Dieter."""
        pole = (self.parameters['rpi'] * 2.16444 / (0.0196 - pt.t) / (1.
                - pt.t/self.parameters['mpi2'])**2 / pt.xi)
        if 'in2particle' in pt and pt.in2particle == 'n':
            return -pole  # neutron
        else:
            return pole  # proton


class DispersionFixedPoleCFF(DispersionCFF, PionPole):
    """Model for CFFs as in arXiv:0904.0458."""

    def __init__(self, **kwargs):
        self.nf = kwargs.setdefault('nf', 4)
        # initial values of parameters and limits on their values
        self.add_parameters({'Nsea': 1.5, 'alS': 1.13, 'alpS': 0.15,
                             'mS2': 0.499849, 'rS': 1.0, 'bS': 2.0,
                             'Nv': 1.35, 'alv': 0.43, 'alpv': 0.85,
                             'mv2': 1.0, 'rv': 0.5, 'bv': 2.2,
                             'C': 7.0, 'mC2': 1.69,
                             'tNv': 0.0, 'tal': 0.43, 'talp': 0.85,
                             'tmv2': 7.29, 'trv': 6.0, 'tbv': 3.0})

        self.add_parameters_limits({'bS': (0.4, 5.0), 'mv2': (0.16, 2.25),
                                    'rv': (0., 8.), 'bv': (0.4, 5.),
                                    'C': (-10., 10.), 'mC2': (0.16, 4.),
                                    'tmv2': (0.16, 4.), 'trv': (0., 8.),
                                    'tbv': (0.4, 5.)})

        super().__init__(**kwargs)

    def subtraction(self, pt):
        """Dispersion relations subtraction constant."""
        return self.parameters['C']/(1.-pt.t/self.parameters['mC2'])**2

    def ImH(self, pt: data.DataPoint, xi: Union[float, np.ndarray] = 0):
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
               onex**p['bv'] / (1. - onex*t/(p['mv2'])))
        sea = ((2./9.) * p['Nsea'] * p['rS'] * twox**(-p['alS']-p['alpS']*t) *
               onex**p['bS'] / (1. - onex*t/(p['mS2']))**2)
        return pi * (val + sea) / (1.+x)

    def ImHt(self, pt: data.DataPoint, xi: Union[float, np.ndarray] = 0):
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
               * twox**regge * onex**p['tbv'] / (1. - onex*t/(p['tmv2'])))
        return pi * val / (1.+x)

    def ImE(self, pt: data.DataPoint, xi: Union[float, np.ndarray] = 0):
        """Imaginary part of CFF E."""
        # Just changing function signature w.r.t. CFF
        # to make it compatible for dispersion integral
        return 0

    def ReEt(self, pt: data.DataPoint):
        """Instead of disp. rel. use pole formula."""
        return self.DMfixpole(pt)


class DispersionFreePoleCFF(DispersionFixedPoleCFF):
    """Model for CFFs as in arXiv:0904.0458. + free pion pole."""

    def __init__(self, **kwargs):
        # Adding two extra parameters:
        self.add_parameters({'rpi': 1.0,  'mpi2': 1.0})

        self.add_parameters_limits({
             'rpi': (-8, 8.),
             'mpi2': (0.16, 16.)})
        super().__init__(**kwargs)

    def ReEt(self, pt):
        """Instead of disp. rel. use pole formula."""
        return self.DMfreepole(pt)


class HybridCFF(MellinBarnesCFF, DispersionFixedPoleCFF, CFF):
    """Combines MB model for small xB and DR model for valence xB.

    Todo:
        Check general consistency for xi != pt.xi for imag CFFs.

    """

    def is_within_model_kinematics(self, pt):
        """Test kinematics of datapoint."""
        # relaxing xBmin and removing Q2max
        return ((1.5 <= pt.Q2) and
                (pt.tm < min(1., pt.Q2/4)) and
                (1e-5 < pt.xB < 0.65))

    # FIXME: this below looks inconsistent generally for xi != pt.xi !!

    def ImH(self, pt, xi=0):
        """Imaginary part of CFF H."""
        return MellinBarnesCFF.ImH(self, pt) + DispersionFixedPoleCFF.ImH(self, pt, xi)

    def ReH(self, pt):
        """Real part of CFF H."""
        return MellinBarnesCFF.ReH(self, pt) + DispersionFixedPoleCFF.ReH(
                self, pt, imfun=DispersionFixedPoleCFF.ImH)

    def ImE(self, pt, xi=0):
        """Imaginary part of CFF E."""
        return MellinBarnesCFF.ImE(self, pt) + DispersionFixedPoleCFF.ImE(self, pt, xi)

    def ReE(self, pt):
        """Real part of CFF E."""
        return MellinBarnesCFF.ReE(self, pt) + DispersionFixedPoleCFF.ReE(
                self, pt, imfun=DispersionFixedPoleCFF.ImE)

    # Tildes are not yet provided by MellinBarnesCFF model

    def ImHt(self, pt, xi=0):
        """Imaginary part of CFF Ht."""
        return DispersionFixedPoleCFF.ImHt(self, pt, xi)

    def ReHt(self, pt):
        """Real part of CFF Ht."""
        return DispersionFixedPoleCFF.ReHt(self, pt, imfun=DispersionFixedPoleCFF.ImHt)

    def ImEt(self, pt):
        """Imaginary part of CFF Et."""
        return DispersionFixedPoleCFF.ImEt(self, pt)

    def ReEt(self, pt):
        """Real part of CFF Et. Provided by subclasses."""
        return 0


class HybridFixedPoleCFF(HybridCFF):
    """HybridCFF model with fixed pion pole model for ReEt."""

    def ReEt(self, pt):
        """Instead of disp. rel. use fixed pion pole formula."""
        return self.DMfixpole(pt)


class HybridFreePoleCFF(HybridCFF):
    """HybridCFF model with free pion pole model for ReEt."""

    def ReEt(self, pt):
        """Instead of disp. rel. use free pion pole formula."""
        return self.DMfreepole(pt)


# This PEP8 violating end-of-file import serves just to bring GK code into cff namespace
from .gk import GoloskokovKrollCFF  # noqa: F401, E402
