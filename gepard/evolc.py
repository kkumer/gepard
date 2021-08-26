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
    """Calculate LO singlet anomalous dimensions matrix on MB contour.

    Args:
       npoints: coordinates of MB contour
            nf: number of active quark flavors

    Returns:
         gam[s,k,i,j]: s in range(npwmax), k in range(npts), i,j in [Q,G]

    """
    # LO only
    gam = []
    for pw_shift in [0, 2, 4]:
        gam.append(g.adim.singlet_LO(npoints+pw_shift, nf, 1).transpose(2, 0, 1))
    return np.array(gam)


def _fshu(j: np.ndarray) -> np.ndarray:
    """Shuvaev factor."""
    #  = 2^(J+1) Gamma(5/2+J) / Gamma(3/2) / Gamma(3+J) =
    #  = 2^N Gamma(3/2+N) / Gamma(3/2) / Gamma(2+N)
    return (2**(j+1) * np.exp(loggamma(2.5 + j)
                              - loggamma(3 + j) - loggamma(3/2)))


def calc_wc(m, j, process: str):
    """Calculate Wilson coeffs for given q2.

    Args:
       m: instance of the model
       j: MB contour point(s) (overrides m.jpoints)
       process: 'DIS, 'DVCS' or 'DVMP'

    Returns:
         wc[k, p, f]: k in range(npts), p in [LO, NLO], f in [Q,G]
    """
    # Instead of type hint (which leads to circular import for some reason)
    if not isinstance(m, g.model.MellinBarnesModel):
        raise Exception("{} is not of type MellinBarnesModel".format(m))
    one = np.ones_like(m.jpoints)
    zero = np.zeros_like(m.jpoints)
    fshu = _fshu(j)
    if process in ['DVCS', 'DIS']:
        if process == 'DVCS':
            quark_norm = fshu
            gluon_norm = fshu
        else:  # DIS
            quark_norm = one
            gluon_norm = one
        q0, g0 = (one, zero)   # LO Q and G
        q1, g1 = (zero, zero)  # NLO if only LO is asked for (m.p=0)
        if m.p == 1:
            # take only singlet part atm:
            q1, g1 = g.c1dvcs.C1(m, j, process)[:, :2].transpose()
    elif process == 'DVMP':
        # Normalizations. Factor 3 is from normalization of DA, so not NC
        # See p. 37, 39 of "Towards DVMP" paper. Eq. (3.62c)
        quark_norm = 3 * fshu
        gluon_norm = 3 * fshu * 2 / g.constants.CF / (j + 3)
        # Normalizations are chosen so that LO is normalized to:
        q0, g0 = (one/m.nf, one)   # LO Q and G
        q1, g1 = (zero, zero)      # NLO if only LO is asked for (m.p=0)
        if m.p == 1:
            qp1, ps1, g1 = g.c1dvmp.c1dvmp(m, 1, j, 0)
            q1 = qp1/m.nf + ps1
    else:
        raise Exception('{} is not DIS, DVCS or DVMP!'.format(process))
    c_quark = quark_norm * np.stack([q0, q1])
    c_gluon = gluon_norm * np.stack([g0, g1])
    return np.stack((c_quark, c_gluon)).transpose()


def calc_wce(m, q2: float, process: str):
    """Calculate evolved Wilson coeffs for given q2, for all PWs.

    Args:
       q2: final evolution scale
       m: instance of the model
       process: 'DIS, 'DVCS' or 'DVMP'

    Returns:
         wce[s,k,f]: s in range(npwmax), k in range(npts), f in [Q,G]
    """
    # Instead of type hint (which leads to circular import for some reason)
    if not isinstance(m, g.model.MellinBarnesModel):
        raise Exception("{} is not of type MellinBarnesModel".format(m))
    wce = []
    for pw_shift in [0, 2, 4]:
        j = m.jpoints + pw_shift
        wc = calc_wc(m, j, process)
        # evolution operator
        evola = g.evolution.evolop(m, j, q2, process)
        # p_mat: matrix that combines (LO, NLO) evolution operator and Wilson coeffs
        # while canceling NNLO term NLO*NLO:
        asmur2 = g.qcd.as2pf(m.p, m.nf, q2/m.rr2, m.asp[m.p], m.r20)
        asmuf2 = g.qcd.as2pf(m.p, m.nf, q2/m.rf2, m.asp[m.p], m.r20)
        p_mat = np.array([[1, asmuf2], [asmur2, 0]])
        # 3. evolved Wilson coeff.
        wce.append(np.einsum('kpi,pq,kqij->kj', wc, p_mat, evola))
    return np.stack(wce, axis=0)  # stack PWs
