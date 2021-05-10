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


def calc_wc(m):
    """Calculate DVCS Wilson coeffs.

    Args:
       m: instance of the model

    Returns:
         wc[s,k,j]: s in range(npwmax), k in range(npts), j in [Q,G]

    """
    # LO only
    # Shuvaev factor = 2^(J+1) Gamma(5/2+J) / Gamma(3/2) / Gamma(3+J) =
    #   = 2^N Gamma(3/2+N) / Gamma(3/2) / Gamma(2+N)
    wc = []
    for pw_shift in [0, 2, 4]:
        fshu = (2.0**(m.npoints + pw_shift)
                * np.exp(loggamma(1.5 + m.npoints + pw_shift)
                         - loggamma(2 + m.npoints + pw_shift))
                / 0.886226925452758014)
        quark = fshu
        gluon = np.zeros_like(quark)
        wc.append(np.array((quark, gluon)).transpose())
    return np.array(wc)


def calc_wc_dvmp(m):
    """Calculate DVMP Wilson coeffs.

    Args:
       m: instance of the model

    Returns:
         wcdvmp[s,k,j]: s in range(npwmax), k in range(npts), j in [Q,G]

    """
    # LO only
    # Shuvaev factor = 2^(J+1) Gamma(5/2+J) / Gamma(3/2) / Gamma(3+J) =
    #   = 2^N Gamma(3/2+N) / Gamma(3/2) / Gamma(2+N)
    wc = []
    for pw_shift in [0, 2, 4]:
        fshu = (2.0**(m.npoints + pw_shift)
                * np.exp(loggamma(1.5 + m.npoints + pw_shift)
                         - loggamma(2 + m.npoints + pw_shift))
                / 0.886226925452758014)
        sea = 3 * fshu / m.nf
        gluon = 3 * fshu * 2 / g.constants.CF / (m.npoints + pw_shift + 2)
        wc.append(np.array((sea, gluon)).transpose())
    return np.array(wc)


def calc_wce(m, q2: float):
    """Calculate evolved Wilson coeffs for given q2.

    Args:
            q2: final evolution scale
       npoints: coordinates of MB contour

    Returns:
         wce[s,k,j]: s in range(npwmax), k in range(npts), j in [Q,G]

    """
    # Instead of type hint (which leads to circular import for some reason)
    if not isinstance(m, g.model.MellinBarnesModel):
        raise Exception("{} is not of type MellinBarnesModel".format(m))
    evola0 = g.evolution.evolop(m.npoints, m.nf, q2, m.q02, m.asp[m.p], m.r20)
    c0 = calc_wc(m)
    return np.einsum('ski,skij->skj', c0, evola0)


def calc_wce_dvmp(m, q2: float):
    """Calculate evolved DVMP Wilson coeffs for given q2.

    Args:
       q2: final evolution scale
       npoints: coordinates of MB contour

    Returns:
         wce_dvmp[s,k,j]: s in range(npwmax), k in range(npts), j in [Q,G]

    """
    # Instead of type hint (which leads to circular import for some reason)
    if not isinstance(m, g.model.MellinBarnesModel):
        raise Exception("{} is not of type MellinBarnesModel".format(m))
    evola0 = g.evolution.evolop(m.npoints, m.nf, q2, m.q02, m.asp[m.p], m.r20)
    c0 = calc_wc_dvmp(m)
    return np.einsum('ski,skij->skj', c0, evola0)
