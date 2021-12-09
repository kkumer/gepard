"""Definitions of Compton Form Factors.

"""
from __future__ import annotations

from math import log, pi, sqrt
from typing import Dict, List

import numpy as np
import scipy.stats

from . import constants, data, eff, gpd, model, qcd, quadrature, utils, wilson

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


class MellinBarnesCFF(ComptonFormFactors):
    """Class of models built by Mellin-Barnes integration."""

    def __init__(self, gpds: gpd.ConformalSpaceGPD) -> None:
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
        self.frot = gpds.frot
        self.dvcs_charges = gpds.dvcs_charges
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
        # correction factors for NLO expressions
        # needed to be able to have some tests w.r.t. old wrong notebooks
        # 1. correction introduced below Eq. (20) of 1612.01937. Set 
        #  to zero to get agreement with older results
        self.corr_c1dvmp_one = 1
        # 2. correction to get results from "Towards DVMP" paper.
        #  Set to -1 to get agreement with Dieter's notebook.
        self.corr_c1dvmp_sgn = 1
        super().__init__()

    def _mellin_barnes_integral(self, xi, wce, gpd):
        """Return convolution of evolved Wilson coefs and GPDs."""
        eph = np.exp(self.phi*1j)
        cfacj = eph * np.exp((self.jpoints + 1) * log(1/xi))  # eph/xi**(j+1)
        # Temporary singlet part only!:
        cch = np.einsum('j,sa,sja,ja->j', cfacj,
                        self.gpds.pw_strengths(), wce, gpd)
        imh = np.dot(self.wg, cch.imag)
        np.multiply(cch, self.tgj, out=cch)
        reh = np.dot(self.wg, cch.imag)
        return reh, imh

    def _dis_mellin_barnes_integral(self, xi, wce, pdf):
        """Return convolution of evolved Wilson coefs and PDFs."""
        eph = np.exp(self.phi*1j)
        cfacj = eph * np.exp((self.jpoints) * log(1/xi))  # eph/xi**j
        # Temporary singlet part only!:
        cch = np.einsum('j,ja,ja->j', cfacj, wce, pdf)
        mb_int = np.dot(self.wg, cch.imag)
        return mb_int

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
        h_prerot = self.gpds.gpd_H(pt.xi, pt.t)
        h = np.einsum('f,fa,ja->jf', self.dvcs_charges, self.frot, h_prerot)
        reh, imh = self._mellin_barnes_integral(pt.xi, wce_ar, h)
        e_prerot = self.gpds.gpd_E(pt.xi, pt.t)
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

    def tff(self, xi: float, t: float, Q2: float) -> np.ndarray:
        """Return array(ReH_rho, ImH_rho, ReE_rho, ...) of DVrhoP transition FFs."""
        assert self.nf == 4

        astrong = 2 * pi * qcd.as2pf(self.p, self.nf,  Q2, self.asp[self.p], self.r20)

        try:
            wce_ar_dvmp = self.wce_dvmp[Q2]
        except KeyError:
            # calculate it
            wce_ar_dvmp = wilson.calc_wce(self, Q2, 'DVMP')
            # memorize it for future
            self.wce_dvmp[Q2] = wce_ar_dvmp
        # Evaluations depending on model parameters:
        h_prerot = self.gpds.gpd_H(xi, t)
        # Flavor rotation matrix: (sea,G,uv,dv) --> (SIG, G, NS+, NS-)
        # FIXME: should be constructed only once!
        frot_rho_4 = np.array([[1, 0, 1, 1],
                               [0, 1, 0, 0],
                               [0., 0, 0., 0.]]) / np.sqrt(2)
                               # [3./20., 0, 5./12., 1./12.]]) / np.sqrt(2)
        h = np.einsum('fa,ja->jf', frot_rho_4, h_prerot)
        reh, imh = self._mellin_barnes_integral(xi, wce_ar_dvmp, h)
        return (constants.CF * constants.F_rho * astrong / constants.NC
                / np.sqrt(Q2) * np.array([reh, imh, 0, 0, 0, 0, 0, 0]))

    def ImHrho(self, pt: data.DataPoint) -> np.ndarray:
        """Return Im(TFF H) for kinematic point."""
        tffs = self.tff(pt.xi, pt.t, pt.Q2)
        return tffs[1]

    def ReHrho(self, pt: data.DataPoint) -> np.ndarray:
        """Return Re(TFF H) for kinematic point."""
        tffs = self.tff(pt.xi, pt.t, pt.Q2)
        return tffs[0]

    def F2(self, pt: data.DataPoint) -> float:
        """Return DIS F2 for kinematic point."""
        if self.nf == 3:
            chargefac = 2./9.
        else:  # nf = 4
            chargefac = 5./18.

        try:
            wce_ar_dis = self.wce_dis[pt.Q2]
        except KeyError:
            # calculate it, first PW is the only relevant one
            wce_ar_dis = wilson.calc_wce(self, pt.Q2, 'DIS')[0, :, :]
            # memorize it for future
            self.wce_dis[pt.Q2] = wce_ar_dis
        pdf_prerot = self.gpds.gpd_H(0, 0)  # forward limit
        # Flavor rotation matrix: (sea,G,uv,dv) --> (SIG, G, NS+, NS-)
        # FIXME: should be constructed only once!
        frot_pdf = np.array([[1, 0, 0, 0],
                             [0, 1, 0, 0],
                             [0., 0, 0., 0.]])
        pdf = np.einsum('fa,ja->jf', frot_pdf, pdf_prerot)
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
        res = quadrature.PVquadrature(self.dispargV, 0, 1, (self.ImH, pt))
        pvpi = (res + log(pt.xi**2 / (1.-pt.xi**2)) * self.ImH(pt)) / pi
        # P.V./pi - subtraction constant C/(1-t/MC^2)^2
        return pvpi - self.subtraction(pt)

    def ReHt(self, pt):
        """Real part of CFF Ht.

        Given by dispersion integral over ImHt.

        """
        res = quadrature.PVquadrature(self.dispargA, 0, 1, (self.ImHt, pt))
        pvpi = (res + log((1.+pt.xi)/(1.-pt.xi)) * self.ImHt(pt))/pi
        return pvpi   # this is P.V./pi

    def ReE(self, pt):
        """Real part of CFF E.

        Given by dispersion integral over ImE plus subtraction constant.

        """
        res = quadrature.PVquadrature(self.dispargV, 0, 1, (self.ImE, pt))
        pvpi = (res + log(pt.xi**2 / (1.-pt.xi**2)) * self.ImE(pt)) / pi
        # This is same subtraction constant
        # as for H, but with opposite sign
        return pvpi + self.subtraction(pt)

    def ReEt(self, pt):
        """Real part of CFF Et.

        Given by dispersion integral over ImEt

        """
        res = quadrature.PVquadrature(self.dispargA, 0, 1, (self.ImEt, pt))
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
        super().__init__()

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

    def __init__(self, instMB: MellinBarnesCFF, instDR: ComptonModelDR, **kwargs):
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
class ModelDR(ComptonModelDR, eff.ElasticDipole):
    """Complete model as in arXiv:0904.0458.."""


class ModelDRKelly(ComptonModelDR, eff.ElasticKelly):
    """Same, but with Kelly elastic form factors."""


class HybridDipole(ComptonHybrid, eff.ElasticDipole):
    """Complete hybrid model."""


class HybridKelly(ComptonHybrid, eff.ElasticKelly):
    """Complete hybrid model."""

