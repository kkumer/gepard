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

import numpy as np

import gepard as g


def lambdaf(npoints: np.ndarray, nf: int) -> np.ndarray:
    """Eigenvalues of the LO singlet anomalous dimensions matrix.

    Returns:
        lam[a, s, k]
        a in [+, -], s is SO(3) PW index, and k is MB contour point index

    """
    # To avoid crossing of the square root cut on the
    # negative real axis we use trick by Dieter Mueller

    gam = g.evolc.calc_gam(npoints, nf)   # anomalous dimensions

    aux = ((gam[:, :, 0, 0] - gam[:, :, 1, 1]) *
           np.sqrt(1. + 4.0 * gam[:, :, 0, 1] * gam[:, :, 1, 0] /
           (gam[:, :, 0, 0] - gam[:, :, 1, 1])**2))
    lam1 = 0.5 * (gam[:, :, 0, 0] + gam[:, :, 1, 1] - aux)
    lam2 = lam1 + aux
    return np.stack([lam1, lam2])


def rnnlof(npoints, nf) -> np.ndarray:
    """Projectors on evolution quark-gluon singlet eigenaxes.

    Also, projected NLO mu-independent R(bet*gam1 - gam0)R is implemented

    Args:
         nf: number of active quark flavors
        sec: index of SO(3) partial wave
          k: index of point on MB contour

    Returns:
          pr: Projector pr[sec, k, a, i, j]
      r1proj (not yet): r1proj[a,b] = sum_ij pr[a,i,j] R1[i,j] pr[b,i,j]
                 a,b in {+,-};  i,j in {Q, G}

    Todo:
        * NLO part is not checked at all

    """
    # cf. my DIS notes p. 61

    lam = lambdaf(npoints, nf)
    den = 1. / (lam[0, :, :] - lam[1, :, :])

    # P+ and P-
    gam = g.evolc.calc_gam(npoints, nf)   # anomalous dimensions
    ssm = gam - np.einsum('sk,ij->skij', lam[1, :, :], np.identity(2))
    ssp = gam - np.einsum('sk,ij->skij', lam[0, :, :], np.identity(2))
    prp = np.einsum('sk,skij->skij', den, ssm)
    prm = np.einsum('sk,skij->skij', -den, ssp)
    pr = np.stack([prp, prm], axis=2)

    # Needed at NLO
    # inv = 1.0 / g.qcd.beta(0, nf)
    # r1 = inv * (gam[:, :, 1, :, :] -
    #             0.5 * inv * g.qcd.beta(1, nf) * gam[:, :, 0, :, :])
    # r1proj = np.einsum('afg,skfg,bfg->skab', pr, r1, pr)

    return pr   # , r1proj


def evolop(npoints: np.ndarray, nf: int,  q2: float, q02: float,
           as0: float, r20: float) -> np.ndarray:
    """GPD evolution operator.

    Args:
         nf: number of active quark flavors
         q2: final evolution momentum squared

    Returns:
         Array corresponding Eq. (121) of Towards NPB paper.
         evolop[s, k, p, i, j]
         -  s is index of SO(3) partial wave,
         -  k is index of point on MB contour,
         -  p is pQCD order (0=LO, 1=NLO)
         -  i, j in [Q, G]

    """
    # FIXME: p=0 hardcoded, all p's should be done simultaneously
    asmuf2 = g.qcd.as2pf(0, nf, q2, as0, r20)
    asq02 = g.qcd.as2pf(0, nf, q02, as0, r20)
    R = asmuf2/asq02

    # 2. errfunc
    lam = lambdaf(npoints, nf)
    b0 = g.qcd.beta(0, nf)
    levi_civita = np.array([[0, 1], [-1, 0]])
    bll = np.einsum('sk,ij->skij', lam[0, :, :] - lam[1, :, :], levi_civita)
    bll = b0 * np.ones_like(bll) + bll
    # er1 = (np.ones_like(bll) - (1./R)**(bll/b0)) / bll  # NLO

    # 3. LO evolution op
    pr = rnnlof(npoints, nf)

    Rfact = R**(-lam/b0)
    evola0ab = np.einsum('skaij,ab->skabij', pr,  np.identity(2))
    evola0 = np.einsum('skabij,bsk->skij', evola0ab, Rfact)

    evola1 = np.zeros_like(evola0)   # temp NLO operator is zero

    evola = np.stack((evola0, evola1), axis=2)

    return evola
