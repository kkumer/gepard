"""General kinematics."""

from __future__ import annotations

import math
from numpy import cos, pi, sqrt

from . import data
from gepard.constants import Mp, Mp2, alpha


def tmin(Q2: float, xB: float, eps2: float) -> float:
    """Minimal momentum transfer. BMK Eq. (31)."""
    return -Q2 * (2. * (1.-xB)*(1. - sqrt(1.+eps2)) + eps2) / (
            4. * xB * (1.-xB) + eps2)


def tmax(Q2: float, xB: float, eps2: float) -> float:
    """Maximal momentum transfer."""
    return -Q2 * (2. * (1.-xB)*(1. + sqrt(1.+eps2)) + eps2) / (
            4. * xB * (1.-xB) + eps2)


def xBmin(s: float, Q2: float) -> float:
    """Constrained by xB=Q2/(s-Mp^2)/yMax, with 1-yMax+yMax^2*eps2/4=0."""
    yMax = 1 + Mp2*Q2/(s-Mp2)**2
    return Q2 / (s-Mp2) / yMax


def K2(Q2: float, xB: float, t: float, y: float, eps2: float) -> float:
    """BMK Eq. (30)."""
    tm = tmin(Q2, xB, eps2)
    brace = sqrt(1.+eps2) + (4. * xB * (1.-xB) + eps2) / (
            4. * (1.-xB)) * (t - tm) / Q2
    return -(t/Q2) * (1.-xB) * (1.-y-y*y*eps2/4.) * (
            1. - tm / t) * brace


def J(Q2: float, xB: float, t: float, y: float, eps2: float) -> float:
    """BMK below Eq. (32)."""
    return (1.-y-y*eps2/2.) * (1. + t/Q2) - (1.-xB)*(2.-y)*t/Q2


def is_within_phase_space(pt: data.DataPoint) -> bool:
    """Check if pt kinematics within allowed phase space."""
    return (pt.xB > xBmin(pt.s, pt.Q2) and pt.t < tmin(pt.Q2, pt.xB, pt.eps2))


def r(Q2: float, xB: float, t: float, y: float, eps2: float) -> float:
    """Dieter's fitting notes, below Eq. (13)."""
    K = sqrt(K2(Q2, xB, t, y, eps2))
    brace = (2.-y)**2 * K / (1.-y) + (1./K)*(t/Q2)*(1.-y)*(2.-xB)
    return - (2.-y) / (2.-2.*y+y**2) * brace


def P1P2(pt: data.DataPoint) -> float:
    """Product of Bethe-Heitler propagators, BMK Eq. (32)."""
    P1 = - (J(pt.Q2, pt.xB, pt.t, pt.y, pt.eps2) + 2. *
            sqrt(K2(pt.Q2, pt.xB, pt.t, pt.y, pt.eps2)) * cos(pt.phi)) / (
                    pt.y * (1. + pt.eps2))
    P2 = 1. + pt.t / pt.Q2 - P1
    return P1 * P2


def anintP1P2(pt: data.DataPoint) -> float:
    """Analitical integral of BH propagators."""
    xB, Q2, t, y, eps2, K2 = pt.xB, pt.Q2, pt.t, pt.y, pt.eps2, pt.K2
    brace = ((1 - y - (1+eps2/2.) * y**2 * eps2/2.) * (1. + t/Q2)**2 +
             2.*K2 - (1.-xB)*(2.-y)**2 * (1. + xB*t/Q2) * t/Q2)
    return -2. * pi * brace / (1+eps2)**2 / y**2


def weight_BH(pt: data.DataPoint) -> float:
    """Weight factor removing BH propagators from INT and BH amplitudes."""
    # It is normalized to int_{0}^{2pi} w  2pi as in BMK
    return 2.*pi*pt.P1P2 / pt.intP1P2


def prepare(pt: data.DataPoint) -> None:
    """Pre-calculate GPD-independent kinamatical constants and functions."""
    if not hasattr(pt, "s"):
        # This is for variable beam energy;  code duplication
        if pt.process in ['ep2epgamma', 'en2engamma']:
            if pt.exptype == 'fixed target':
                pt.s = 2 * Mp * pt.in1energy + Mp2
            elif pt.exptype == 'collider':
                pt.s = 2 * pt.in1energy * (pt.in2energy + math.sqrt(
                    pt.in2energy**2 - Mp2)) + Mp2
            else:
                pass  # FIXME: raise error
        else:
            pass  # FIXME: should I raise error here?
    pt.y = (pt.W**2 + pt.Q2 - Mp2) / (pt.s - Mp2)
    pt.eps = 2. * pt.xB * Mp / sqrt(pt.Q2)
    pt.eps2 = pt.eps**2
    if 't' in pt:
        pt.J = J(pt.Q2, pt.xB, pt.t, pt.y, pt.eps2)
        pt.K2 = K2(pt.Q2, pt.xB, pt.t, pt.y, pt.eps2)
        pt.K = sqrt(pt.K2)
        pt.tK2 = pt.K2*pt.Q2/(1-pt.y-pt.eps2*pt.y**2/4.)
        pt.tK = sqrt(pt.tK2)
        pt.r = r(pt.Q2, pt.xB, pt.t, pt.y, pt.eps2)
        # Needed for BMP higher-twist stuff:
        pt.chi0 = sqrt(2.*pt.Q2)*pt.tK/sqrt(1+pt.eps2)/(pt.Q2+pt.t)
        pt.chi = (pt.Q2-pt.t+2.*pt.xB*pt.t)/sqrt(1+pt.eps2)/(pt.Q2+pt.t) - 1.
        # First option is numerical, second is analytical and faster
        # pt.intP1P2 = g.quadrature.Hquadrature(lambda phi: P1P2(pt, phi), 0, 2.0*pi)
        pt.intP1P2 = anintP1P2(pt)
    if 'phi' in pt:
        pt.P1P2 = P1P2(pt)


def long2trans(pt: data.DataPoint) -> float:
    """Ratio of longitudinal to transverse photon flux 1304.0077 Eq. (2.9)."""
    return (1.-pt.y-pt.eps2*pt.y**2/4.)/(1-pt.y+pt.y**2/2+pt.eps2*pt.y**2/4.)


def HandFlux(pt: data.DataPoint) -> float:
    """Virtual photon flux (Hand convention) 1304.0077 Eq. (2.9)."""
    return (alpha/2./pi)*(pt.y**2/(1.-long2trans(pt)))*(1-pt.xB)/pt.xB/pt.Q2
