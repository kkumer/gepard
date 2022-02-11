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

    def _XGAMMA_rho0_t_Approx(self, pt: data.DataPoint) -> float:
        """Partial (longitudinal) gamma* p -> rho0 p cross section w.r.t. Mandelstam t.

        Args:
            pt: instance of DataPoint

        Returns:
            Cross-section differential in t.
            Approximate formula valid for small xB.

        """
        # 4 * pi**2 * alpha_em * GeV2nb = 112175.5
        res = 112175.5 * pt.xB**2 * (
                self.m.ImH_rho0(pt)**2 + self.m.ReH_rho0(pt)**2) / pt.Q2**2
        return res

    _XGAMMA_rho0_t = _XGAMMA_rho0_t_Approx


class MellinBarnesTFF(model.ParameterModel):
    """DVMP Transition Form Factors modelled as Mellin-Barnes integral."""

    def __init__(self, **kwargs):
        self.wce_dvmp: Dict[float, np.ndarray] = {}
        # correction factors for NLO expressions
        # needed to be able to have some tests w.r.t. old wrong notebooks
        # 1. correction introduced below Eq. (20) of 1612.01937. Set
        #  to zero to get agreement with older results
        self.corr_c1dvmp_one = 1
        # 2. correction to get results from "Towards DVMP" paper.
        #  Set to -1 to get agreement with Dieter's notebook.
        self.corr_c1dvmp_sgn = 1
        super().__init__(**kwargs)

    def tff(self, xi: float, t: float, Q2: float) -> np.ndarray:
        """Return array(ReH_rho0, ImH_rho0, ReE_rho0, ...) of DVrho0P transition FFs."""
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
        h_prerot = self.H(xi, t)
        h = np.einsum('fa,ja->jf', self.frot_rho0_4, h_prerot)
        reh, imh = self._mellin_barnes_integral(xi, wce_ar_dvmp, h)
        return (constants.CF * constants.F_rho0 * astrong / constants.NC
                / np.sqrt(Q2) * np.array([reh, imh, 0, 0, 0, 0, 0, 0]))

    def ImH_rho0(self, pt: data.DataPoint) -> np.ndarray:
        """Return Im(TFF H) for kinematic point."""
        tffs = self.tff(pt.xi, pt.t, pt.Q2)
        return tffs[1]

    def ReH_rho0(self, pt: data.DataPoint) -> np.ndarray:
        """Return Re(TFF H) for kinematic point."""
        tffs = self.tff(pt.xi, pt.t, pt.Q2)
        return tffs[0]
