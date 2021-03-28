"""Components of the evolution operator mathcal{E}."""

from cmath import sqrt
from typing import Tuple

import numpy as np

import gepard as g
import gepard.pygepard as gfor


def lambdaf(sec: int, k: int) -> Tuple[complex, complex]:
    """Eigenvalues of the LO singlet anomalous dimensions matrix.

    Args:
        sec: index of SO(3) partial wave
          k: index of point on MB contour

    Returns:
        two eigenvalues (lam_+, lam_-)

    """
    # To avoid crossing of the square root cut on the
    # negative real axis we use trick by Dieter Mueller

    # gam(sec, k, p, flav1, flav2)
    gam = gfor.gam.gam   # pre-calculated anomalous dimensions

    aux = ((gam[sec, k, 0, 0, 0] - gam[sec, k, 0, 1, 1]) *
           sqrt(1. + 4.0 * gam[sec, k, 0, 0, 1] * gam[sec, k, 0, 1, 0] /
                (gam[sec, k, 0, 0, 0] - gam[sec, k, 0, 1, 1])**2))

    lam1 = 0.5 * (gam[sec, k, 0, 0, 0] + gam[sec, k, 0, 1, 1] - aux)
    lam2 = lam1 + aux

    return lam1, lam2


def rnnlof(nf: int, sec: int, k: int) -> Tuple[np.ndarray, np.ndarray]:
    """Projectors on evolution quark-gluon singlet eigenaxes.

    Also, projected NLO mu-independent R(bet*gam1 - gam0)R is returned

    Args:
         nf: number of active quark flavors
        sec: index of SO(3) partial wave
          k: index of point on MB contour

    Returns:
          pr: Projector pr[a,i,j]
      r1proj: r1proj[a,b] = sum_ij pr[a,i,j] R1[i,j] pr[b,i,j]
                 a,b in {+,-};  i,j in {Q, G}

    Todo:
        * NLO part is not checked at all

    """
    # cf. my DIS notes p. 61

    lamp, lamm = lambdaf(sec, k)
    den = 1. / (lamp - lamm)

    # P+ and P-
    gam = gfor.gam.gam
    prp = den * (gam[sec, k, 0, :, :] - lamm * np.identity(2))
    prm = np.identity(2) - prp
    pr = np.stack([prp, prm])   # full projector 3D matrix

    inv = 1.0 / g.qcd.beta(0, nf)

    r1 = inv * (gam[sec, k, 1, :, :] -
                0.5 * inv * g.qcd.beta(1, nf) * gam[sec, k, 0, :, :])

    r1proj = np.einsum('afg,fg,bfg->ab', pr, r1, pr)

    return pr, r1proj


def erfunc(p: int, nf: int,  q2: float, sec: int, k: int) -> Tuple[complex, np.ndarray]:
    """Parts of evolution operator dependent on q2.

    Args:
          p: pQCD order, 0=LO, 1=NLO, ...
         nf: number of active quark flavors
         q2: final evolution momentum squared
        sec: index of SO(3) partial wave
          k: index of point on MB contour

    Returns:
          (R, erfunc) where R is ratio of alpha-strongs, and erfunc is
          array from Eq. (126) of Towards NPB paper.
          erfunc[a, b], with a,b in [+, -]

    """
    lamp, lamm = lambdaf(sec, k)
    b0 = g.qcd.beta(0, nf)

    levi_civita = np.array([[0, 1], [-1, 0]])
    bll = b0 * np.ones((2, 2)) + (lamp-lamm)*levi_civita

    r20 = gfor.astrong.mu02  # 2.5
    as0 = gfor.astrong.asp[0]  # 0.0606
    q02 = gfor.parflt.q02

    asmuf2 = g.qcd.as2pf(p, nf, q2, as0, r20)
    asq02 = g.qcd.as2pf(p, nf, q02, as0, r20)
    R = asmuf2/asq02

    er1 = (1. - (1./R)**(bll/b0)) / bll

    return R, er1


def evola(p: int, nf: int,  q2: float, sec: int, k: int) -> np.ndarray:
    """GPD evolution operator terms.

    Args:
          p: pQCD order, 0=LO, 1=NLO, ...
         nf: number of active quark flavors
         q2: final evolution momentum squared
        sec: index of SO(3) partial wave
          k: index of point on MB contour

    Returns:
          Array which is term from Eq. (121) of Towards NPB paper.
          evola[i, j], with i,j in [Q, G]

    """
    pr, r1proj = rnnlof(nf, sec, k)
    R, er1 = erfunc(p, nf,  q2, sec, k)
    lamp, lamm = lambdaf(sec, k)
    lam = np.stack([lamp, lamm])   # Should work with this all the time!
    b0 = g.qcd.beta(0, nf)

    assert p == 0   # Only LO implemented ATM

    Rfact = R**(-lam/b0)
    evola0ab = np.einsum('aij,ab->abij', pr,  np.identity(2))
    evola0 = np.einsum('abij,b->ij', evola0ab, Rfact)

    return evola0
