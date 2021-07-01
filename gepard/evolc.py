"""Initialization of evolved Wilson coefficients.

Todo:
    Complete docoupling from Fortran

"""
# import pickle, sys
#
# import logzero
# _lg = logzero.logger
#
# from numpy import log, pi, imag, real, sqrt, cos, sin, exp
# from numpy import ndarray, array, sum, loadtxt
# import scipy.stats
# from scipy.special import j0, j1, gamma, beta
# from scipy.interpolate import SmoothBivariateSpline
#
# from quadrature import PVquadrature, bquadrature, rthtquadrature
# from utils import flatten, hubDictNew, stringcolor
# from constants import tolerance2, GeVfm, Mp
# import dispersion as DR

import sys

import numpy as np
from scipy.special import loggamma  # type: ignore

import gepard as g

sys.path.append('.')


def calc_gam(npoints, nf):
    """Calculate (singlet) anomalous dimensions.

    Args:
       npoints: coordinates of MB contour
            nf: number of active quark flavors

    Returns:
         gam[s,k,i,j]: s in range(npwmax), k in range(npts), i,j in [Q,G]

    """
    # LO only
    gam = []
    for pw_shift in [0, 2, 4]:
        gam.append(g.adim.singlet(npoints+pw_shift, 0, nf, 1).transpose(2, 0, 1))
    return np.array(gam)


def _fshu(j: np.ndarray) -> np.ndarray:
    """Shuvaev factor."""
    #  = 2^(J+1) Gamma(5/2+J) / Gamma(3/2) / Gamma(3+J) =
    #  = 2^N Gamma(3/2+N) / Gamma(3/2) / Gamma(2+N)
    return (2**(j+1) * np.exp(loggamma(2.5 + j)
                              - loggamma(3 + j) - loggamma(3/2)))


def calc_wc(m, process):
    """Calculate DVCS or DVMP Wilson coeffs.

    Args:
       m: instance of the model
       process: 'DVCS' or 'DVMP'

    Returns:
         wc[s,k,j]: s in range(npwmax), k in range(npts), j in [Q,G]

    """
    wc = []
    for pw_shift in [0, 2, 4]:
        fshu = _fshu(m.jpoints + pw_shift)
        if process == 'DVCS':
            quark = fshu
            gluon = np.zeros_like(quark)
        elif process == 'DVMP':
            quark = 3 * fshu / m.nf
            gluon = 3 * fshu * 2 / g.constants.CF / (m.jpoints + pw_shift + 3)
        else:
            raise Exception('{} is not DVCS or DVMP!'.format(process))
        wc.append(np.array((quark, gluon)).transpose())
    return np.array(wc)


def calc_wce(m, q2: float, process: str):
    """Calculate evolved Wilson coeffs for given q2.

    Args:
            q2: final evolution scale
             m: instance of the model
       process: 'DVCS' or 'DVMP'

    Returns:
         wce[s,k,j]: s in range(npwmax), k in range(npts), j in [Q,G]

    """
    # Instead of type hint (which leads to circular import for some reason)
    if not isinstance(m, g.model.MellinBarnesModel):
        raise Exception("{} is not of type MellinBarnesModel".format(m))
    evola0 = g.evolution.evolop(m.npoints, m.nf, q2, m.q02, m.asp[m.p], m.r20)
    c0 = calc_wc(m, process)
    return np.einsum('ski,skij->skj', c0, evola0)
