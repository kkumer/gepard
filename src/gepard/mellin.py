"""Mellin-Barnes integrations."""

from math import log, pi

import numpy as np

class MellinBarnes(object):
    """Base class for models built by Mellin-Barnes integration."""

    def __init__(self, **kwargs) -> None:
        """Init base MellinBarnes object."""
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

    def _mellin_barnes_integral_HE(self, xi, wce, gpd_H, gpd_E):
        """Return convolution of evolved Wilson coefs and GPDs H and E."""
        eph = np.exp(self.phi*1j)
        cfacj = eph * np.exp((self.jpoints + 1) * log(1/xi))  # eph/xi**(j+1)
        cch_H = np.einsum('j,sa,sja,ja->j', cfacj,
                        self.pw_strengths(), wce, gpd_H)
        cch_E = np.einsum('j,sa,sja,ja->j', cfacj,
                        self.pw_strengths_E(), wce, gpd_E)
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
