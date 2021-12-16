"""Deeply Virtual Meson Production (DVMP) observables."""

from . import theory

class DVMP(theory.Theory):
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
