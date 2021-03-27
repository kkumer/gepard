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

from numpy import array

import gepard.pygepard as gfor

sys.path.append('.')


def evol_init(p=0, scheme='MSBAR', nf=4, q02=4.0):
    """Initialize fortran evolution code."""
    # Choice for ansatz should be irrelevant since it is provided by Python!
    # this was in Gepard's GEPARD.INI, which is not needed now
    # but look at it for documentation of what parameters below are
    # importlib.reload(gfor)  # hoping for cleanup
    gfor.parint.speed = 1
    gfor.parint.acc = 3
    gfor.parint.p = p
    gfor.parint.nf = nf
    gfor.parint.czero = 1
    gfor.parint.pid = 1   # FIXME: should go away somewhere

    gfor.astrong.mu02 = 2.5
    gfor.astrong.asp = array([0.0606, 0.0518, 0.0488])

    gfor.parflt.q02 = q02
    gfor.parflt.rf2 = 1.0
    gfor.parflt.rr2 = 1.0
    gfor.parflt.rdaf2 = 1.0
    gfor.parflt.rgpdf2 = 1.0

    gfor.mbcont.c = 0.35
    gfor.mbcont.phi = 1.57079632
    gfor.mbcont.cnd = -0.25
    gfor.mbcont.phind = 1.57

    ansatz = 'TEST'  # Irrelevant since we need only evolved Wilsons and weights!
    gfor.parchr.ansatz = array([c for c in ansatz + (6-len(ansatz))*' '])  # array(6)
    gfor.parchr.scheme = array([c for c in scheme + (5-len(scheme))*' '])  # array(5)

    gfor.init()
    # number of points on MB contour
    gfor.npts = 2**gfor.parint.acc * 12 / gfor.parint.speed
    return gfor


# def _GepardFFs(self, pt, FF='cfff'):
#     """Call gepard routine that calculates CFFs or TFFs."""
#     for i in self.parameters_index:
#         gfor.par.par[i-1] = self.parameters[self.parameters_index[i]]
#
#     gfor.kinematics.q2 = pt.Q2
#     gfor.kinematics.xi = pt.xi
#     gfor.kinematics.del1 = pt.t
#
#     gfor.mt.nmts = 1
#     gfor.mt.mtind = 0
#     gfor.mts.mts[0] = - pt.t
#
#     getattr(gfor, FF)()
#     gfor.newcall = 0
#
# def DISF2(self, pt):
#     """Call gepard routine that calculates DIS F2."""
#     for i in self.parameters_index:
#         gfor.par.par[i-1] = self.parameters[self.parameters_index[i]]
#     gfor.parint.pid = 0
#     gfor.kinematics.q2 = pt.Q2
#     # Note a hack in Fortran gepard where for F2 calculation
#     # XI should actually be xB, and not xB/2 !
#     gfor.kinematics.xi = 2*pt.xi/(1.+pt.xi)
#     gfor.kinematics.del2 = 0
#     gfor.f2f()
#     return gfor.f2.f2[gfor.parint.p]
