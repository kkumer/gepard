"""Deeply Virtual Meson Production (DVMP) observables."""
from __future__ import annotations

from math import pi
from typing import Dict
import numpy as np

from . import constants, data, model, qcd, theory, wilson


class DVMP(theory.Theory):
    """DVMP observables.

    Implements cross-section for electroproduction of meson.
    """

    def _XGAMMA_DVMP_t_Approx(self, pt: data.DataPoint) -> float:
        """Partial (longitudinal) gamma* p -> V p cross section w.r.t. Mandelstam t.

        Args:
            pt: instance of DataPoint

        Returns:
            Cross-section differential in t.
            Approximate formula valid for small xB.

        """
        # 4 * pi**2 * alpha_em * GeV2nb = 112175.5
        res = 112175.5 * pt.xB**2 * (
                self.m.ImH_V(pt)**2 + self.m.ReH_V(pt)**2) / pt.Q2**2
        return res

    _XGAMMA_DVMP_t = _XGAMMA_DVMP_t_Approx


class MellinBarnesTFF(model.ParameterModel):
    """DVMP Transition Form Factors modelled as Mellin-Barnes integral."""

    def __init__(self, **kwargs):
        self.wce_dvmp: Dict[float, np.ndarray] = {}
        self.wce_dvmp2: Dict[float, np.ndarray] = {}
        # correction factors for NLO expressions
        # needed to be able to have some tests w.r.t. old wrong notebooks
        # 1. correction introduced below Eq. (20) of 1612.01937. Set
        #  to zero to get agreement with older results
        self.corr_c1dvmp_one = 1
        # 2. Prefactors to kill parts of NLO hard scattering amplitude
        #  Set to zero to kill corresponding part.
        self.Q1_prefac = 1
        self.PS1_prefac = 1
        self.G1_prefac = 1
        super().__init__(**kwargs)

    def tff(self, xi: float, t: float, Q2: float, meson: str = 'rho0') -> np.ndarray:
        """Return array(ReH_V, ImH_V, ReE_V, ...) of DVMP transition FFs.

        Only vector mesons are implemented, thus _V.
        """
        assert self.nf == 4

        astrong = 2 * pi * qcd.as2pf(self.p, self.nf,  Q2, self.asp[self.p], self.r20)

        try:
            wce_ar_dvmp = self.wce_dvmp[Q2]
            wce_ar_dvmp2 = self.wce_dvmp2[Q2]
        except KeyError:
            # calculate it
            wce_ar_dvmp = wilson.calc_wce(self, Q2, 'DVMP', 0) # added DA moment "kk"
            wce_ar_dvmp2 = wilson.calc_wce(self, Q2, 'DVMP', 2) # added DA moment "kk"

            # memorize it for future
            self.wce_dvmp[Q2] = wce_ar_dvmp
            self.wce_dvmp2[Q2] = wce_ar_dvmp2
        # Evaluations depending on model parameters:
        h_prerot = self.H(xi, t)
        if meson == 'rho0':
            frot = self.frot_rho0_4
            FV = constants.F_rho0
        elif meson == 'phi':
            frot = self.frot_phi_4
            FV = constants.F_phi
        else:
            raise ValueError("{} unknown. Use 'rho0' or 'phi'".format(meson))
        h = np.einsum('fa,ja->jf', frot, h_prerot)
        reh0, imh0 = self._mellin_barnes_integral(xi, wce_ar_dvmp, h)
        reh2, imh2 = self._mellin_barnes_integral(xi, wce_ar_dvmp2, h)
        reh2 *=self.parameters['a2']
        imh2 *=self.parameters['a2']
        reh = reh0 + reh2
        imh = imh0 + imh2
        #print(self.parameters['a2'])
        #print(self._mellin_barnes_integral(xi, wce_ar_dvmp, h))
        #print(reh,imh)
        return (constants.CF * FV * astrong / constants.NC
                / np.sqrt(Q2) * np.array([reh, imh, 0, 0, 0, 0, 0, 0]))
         
    def ImH_V(self, pt: data.DataPoint) -> np.ndarray:
        """Return Im(TFF H) for kinematic point."""
        if pt.process == 'gammastarp2rho0p':
            tffs = self.tff(pt.xi, pt.t, pt.Q2, 'rho0')
        elif pt.process == 'gammastarp2phip':
            tffs = self.tff(pt.xi, pt.t, pt.Q2, 'phi')
        return tffs[1]

    def ReH_V(self, pt: data.DataPoint) -> np.ndarray:
        """Return Re(TFF H) for kinematic point."""
        if pt.process == 'gammastarp2rho0p':
            tffs = self.tff(pt.xi, pt.t, pt.Q2, 'rho0')
        elif pt.process == 'gammastarp2phip':
            tffs = self.tff(pt.xi, pt.t, pt.Q2, 'phi')
        return tffs[0]
