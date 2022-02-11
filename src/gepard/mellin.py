"""Mellin-Barnes integrations.

Todo:
    Various MB routines should be integrated into one.

"""

from math import log, pi

import numpy as np


class MellinBarnes(object):
    """Base class for models built by Mellin-Barnes integration."""

    def __init__(self, **kwargs) -> None:
        self.tgj = np.tan(pi*self.jpoints/2)
        # print('MellinBarnes init done.')

    def _mellin_barnes_integral(self, xi, wce, gpd):
        """Return convolution of evolved Wilson coefs and GPDs."""
        eph = np.exp(self.phi*1j)
        cfacj = eph * np.exp((self.jpoints + 1) * log(1/xi))  # eph/xi**(j+1)
        cch = np.einsum('j,sa,sja,ja->j', cfacj,
                        self.pw_strengths(), wce, gpd)
        imh = np.dot(self.wg, cch.imag)
        np.multiply(cch, self.tgj, out=cch)
        reh = np.dot(self.wg, cch.imag)
        return reh, imh

    def _mellin_barnes_integral_HE(self, xi, wce, H, E):
        """Return convolution of evolved Wilson coefs and GPDs H and E."""
        eph = np.exp(self.phi*1j)
        cfacj = eph * np.exp((self.jpoints + 1) * log(1/xi))  # eph/xi**(j+1)
        cch_H = np.einsum('j,sa,sja,ja->j', cfacj,
                          self.pw_strengths(), wce, H)
        cch_E = np.einsum('j,sa,sja,ja->j', cfacj,
                          self.pw_strengths_E(), wce, E)
        imh = np.dot(self.wg, cch_H.imag)
        ime = np.dot(self.wg, cch_E.imag)
        np.multiply(cch_H, self.tgj, out=cch_H)
        np.multiply(cch_E, self.tgj, out=cch_E)
        reh = np.dot(self.wg, cch_H.imag)
        ree = np.dot(self.wg, cch_E.imag)
        return reh, imh, ree, ime

    def _dis_mellin_barnes_integral(self, xi, wce, pdf):
        """Return convolution of evolved Wilson coefs and PDFs."""
        eph = np.exp(self.phi*1j)
        cfacj = eph * np.exp((self.jpoints) * log(1/xi))  # eph/xi**j
        cch = np.einsum('j,ja,ja->j', cfacj, wce, pdf)
        mb_int = np.dot(self.wg, cch.imag)
        return mb_int

    def _j2x_mellin_barnes_integral(self, x, eta, wce, gpd):
        """Return convolution of j->x coef, evolution operator and GPD."""
        # difference wrt above integrations is that here we do NOT sum over flavors
        eph = np.exp(self.phi*1j)
        cfacj = eph * np.exp((self.jpoints + 1) * log(1/x))  # eph/x**(j+1)
        if eta < 1e-8:
            # forward limit, PDF-like, so only zero-th PW is taken
            cch = np.einsum('j,ja,ja->ja', cfacj, wce[0, :, :], gpd)
        elif abs(eta-x) < 1e-8:
            # cross-over, border eta=x limit
            cch = np.einsum('j,sa,sja,ja->ja', cfacj,
                            self.pw_strengths(), wce, gpd)
        else:
            raise Exception('eta has to be either 0 or equal to x')
        mb_int_flav = np.dot(self.wg, cch.imag)
        return mb_int_flav

    def _j2x_mellin_barnes_integral_E(self, x, eta, wce, gpd):
        """Return convolution of j->x coef, evolution operator and GPD E."""
        # difference wrt above integrations is that here we do NOT sum over flavors
        eph = np.exp(self.phi*1j)
        cfacj = eph * np.exp((self.jpoints + 1) * log(1/x))  # eph/x**(j+1)
        if eta < 1e-8:
            # forward limit, PDF-like, so only zero-th PW is taken
            cch = np.einsum('j,ja,ja->ja', cfacj, wce[0, :, :], gpd)
        elif abs(eta-x) < 1e-8:
            # cross-over, border eta=x limit
            cch = np.einsum('j,sa,sja,ja->ja', cfacj,
                            self.pw_strengths_E(), wce, gpd)
        else:
            raise Exception('eta has to be either 0 or equal to x')
        mb_int_flav = np.dot(self.wg, cch.imag)
        return mb_int_flav
