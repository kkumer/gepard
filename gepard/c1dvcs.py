"""NLO DVCS Wilson coefficients for vector case."""

import math

import numpy as np

from gepard.adim import block
from gepard.constants import CF
from gepard.special import S1, S2


def c1_F2(n: complex, nf: int) -> np.ndarray:
    """c_F2 NLO non-singlet Wilson coefficient."""
    NSP = CF*(-9+2/n**2+3/n+4/(1+n)+(3-2/(n*(1+n)))*S1(n)-2*S2(n)+2*S1(n)**2)/2
    NSM = NSP
    Q = NSP
    G = nf*(n**(-2)+4/(1+n)-4/(2+n)-((1+S1(n))*(2+n+n**2))/(n*(1+n)*(2+n)))
    return np.array((Q, G, NSP, NSM))


def c1_FL(n: complex, nf: int) -> np.ndarray:
    """c_FL NLO non-singlet Wilson coefficient."""
    NSP = (2*CF)/(1+n)
    NSM = NSP
    Q = NSP
    G = (4*nf)/((1+n)*(2+n))
    return np.array((Q, G, NSP, NSM))


def c1_F1(n: complex, nf: int) -> np.ndarray:
    """c_F1 NLO non-singlet Wilson coefficient."""
    return c1_F2(n, nf) - c1_FL(n, nf)


def shift1(m, j: np.ndarray, process: str) -> np.ndarray:
    """Calculate NLO shift coeff s_1.

    Args:
       m: instance of g.model.MellinBarnesModel
       j: MB contour points
       process: 'DVCS' or 'DIS'

    """
    LRF2 = math.log(m.rf2)
    if process == 'DIS':
        s1 = - LRF2*np.ones_like(j)
    elif process == 'DVCS':    # Eq. (88a)
        s1 = S1(j+3/1) - S1(j+2) + 2*math.log(2) - LRF2
    else:
        raise Exception('Process {} is neither DVCS nor DIS!'.format(process))
    return s1


# def shift2(m, j: np.ndarray, process: str) -> np.ndarray:
#     """Calculate NNLO shift coeff s_2.
#
#     Args:
#        m: instance of g.model.MellinBarnesModel
#        j: MB contour points
#        process: 'DVCS' or 'DIS'
#
#     """
#     if process == 'DIS':
#         s2 = shift1(m, j, process)**2
#     elif process == 'DVCS':    # Eq. (88b)
#         s2 = shift1(m, j, process)**2 - S2(j+3/2) + S2(j+2)
#     else:
#         raise Exception('Process {} is neither DVCS nor DIS!'.format(process))
#     return s2


def C1(m, j: np.ndarray, process: str) -> np.ndarray:
    """Calculate NLO Wilson coeff C_1 for DVCS or DIS.

    Args:
       m: instance of g.model.MellinBarnesModel
       j: MB contour points
       process: 'DVCS' or 'DIS'

    "Big C" from eqs. (91) and (101) from "Towards ... DVCS." paper
    """
    c0 = np.array([1, 0, 1, 1])  # LO, valid for DIS and DVCS
    if process == 'DIS':
        c1 = c1_F2(j+1, m.nf)
    elif process == 'DVCS':
        c1 = c1_F1(j+1, m.nf)
    else:
        raise Exception('Process {} is neither DVCS nor DIS!'.format(process))
    return c1 + shift1(m, j, process)*np.einsum('i,ij->j', c0,
                                                block(j+1, m.nf)[0, :, :])/2
