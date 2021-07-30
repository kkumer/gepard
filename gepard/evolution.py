"""GPD evolution operator mathcal{E}.

Returns:
   Numpy array evol[p, k, j]  where flavor index j takes values
   (1) -- singlet quark
   (2) -- gluon
   (3) -- NS(+)
   (4) -- NS(-)  (not tested!)

Notes:
   Careful: GPD models may have u_v and d_v at slots 3 and 4.

Todo:
    * Implement non-singlet
    * Implement NLO and NNLO
    * Implement MSBAR scheme

"""

from typing import Tuple

import numpy as np

import gepard as g


def lambdaf(m, j) -> np.ndarray:
    """Eigenvalues of the LO singlet anomalous dimensions matrix.

    Args:
          m: instance of the model
          j: MB contour points (overrides m.jpoints)

    Returns:
        lam[a, k]
        a in [+, -] and k is MB contour point index

    """
    # To avoid crossing of the square root cut on the
    # negative real axis we use trick by Dieter Mueller

    gam0 = g.adim.singlet_LO(j+1, m.nf)   # LO singlet anomalous dimensions

    aux = ((gam0[0, 0, :] - gam0[1, 1, :]) *
           np.sqrt(1. + 4.0 * gam0[0, 1, :] * gam0[1, 0, :] /
                   (gam0[0, 0, :] - gam0[1, 1, :])**2))
    lam1 = 0.5 * (gam0[0, 0, :] + gam0[1, 1, :] - aux)
    lam2 = lam1 + aux
    return np.stack([lam1, lam2])


def rnnlof(m, j) -> Tuple[np.ndarray, np.ndarray]:
    """Projectors on evolution quark-gluon singlet eigenaxes.

    Also, projected NLO mu-independent R(bet*gam1 - gam0)R is implemented

    Args:
          m: instance of the model
          j: MB contour points (overrides m.jpoints)

    Returns:
          pr: Projector pr[k, a, i, j]  # Eq. (122)
      r1proj: r1proj[a,b] = sum_ij pr[a,i,j] R1[i,j] pr[b,i,j]  # Eq. (124)
                 a,b in {+,-};  i,j in {Q, G}
    """
    # cf. my DIS notes p. 61

    lam = lambdaf(m, j)
    den = 1. / (lam[0, :] - lam[1, :])
    gam0 = g.adim.singlet_LO(j+1, m.nf).transpose((2, 0, 1))  # LO singlet anomalous dimensions
    gam1 = g.adim.singlet_NLO(j+1, m.nf).transpose((2, 0, 1))  # NLO singlet anomalous dimensions

    # P+ and P-
    ssm = gam0 - np.einsum('k,ij->kij', lam[1, :], np.identity(2))
    ssp = gam0 - np.einsum('k,ij->kij', lam[0, :], np.identity(2))
    prp = np.einsum('k,kij->kij', den, ssm)
    prm = np.einsum('k,kij->kij', -den, ssp)
    pr = np.stack([prp, prm], axis=1)


    # NLO
    inv = 1.0 / g.qcd.beta(0, m.nf)
    r1 = inv * (gam1 -
                0.5 * inv * g.qcd.beta(1, m.nf) * gam0)
    r1proj = np.einsum('kaim,kmn,kbnj->kabij', pr, r1, pr)

    return pr, r1proj


def evolop(m, j, q2: float) -> np.ndarray:
    """GPD evolution operator.

    Args:
         m: instance of the model
         j: MB contour points (overrides m.jpoints)
         q2: final evolution momentum squared

    Returns:
         Array corresponding Eq. (121) of Towards NPB paper.
         evolop[k, p, i, j]
         -  k is index of point on MB contour,
         -  p is pQCD order (0=LO, 1=NLO)
         -  i, j in [Q, G]
    """
    # 1. Alpha-strong ratio.
    # When m.p=1 (NLO), LO part of the evolution operator
    # will still be multiplied by ratio of alpha_strongs
    # evaluated at NLO, as it should.
    asmuf2 = g.qcd.as2pf(m.p, m.nf, q2, m.asp[m.p], m.r20)
    asq02 = g.qcd.as2pf(m.p, m.nf, m.q02, m.asp[m.p], m.r20)
    R = asmuf2/asq02

    # 2. LO errfunc
    lam = lambdaf(m, j)
    b0 = g.qcd.beta(0, m.nf)
    levi_civita = np.array([[0, 1], [-1, 0]])
    bll = np.einsum('k,ij->kij', lam[0, :] - lam[1, :], levi_civita)
    bll = b0 * np.ones_like(bll) + bll

    er1 = (np.ones_like(bll) - (1./R)**(bll/b0)) / bll  # Eq. (126)

    # 3. projectors
    pr, r1proj = rnnlof(m, j)

    Rfact = R**(-lam/b0)  # LO evolution (alpha(mu)/alpha(mu0))^(-gamma/beta0)
    evola0ab = np.einsum('kaij,ab->kabij', pr,  np.identity(2))
    evola0 = np.einsum('kabij,bk->kij', evola0ab, Rfact)

    if m.p == 1:
        evola1ab = np.einsum('kij,kabij->kabij', er1, r1proj)  # Eq. (124)
        evola1 = np.einsum('kabij,bk->kij', evola1ab, Rfact)
    else:
        evola1 = np.zeros_like(evola0)

    evola = np.stack((evola0, evola1), axis=1)

    return evola
