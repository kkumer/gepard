"""GPD evolution operator mathcal{E}.

Returns:
   Numpy array evol[p, k, j]  where flavor index j takes values
   in the evolution basis:
   (1) -- singlet quark
   (2) -- gluon
   (3) -- NS(+)
   (4) -- NS(-)  (not tested!)

Notes:
   GPD models may be defined in different basis and should
   provide appropriate transformation matrix

Todo:
    * Array indices ordering is a bit of a mess

"""

from typing import Tuple

import numpy as np

from . import adim, constants, evolution, qcd, quadrature, special


def lambdaf(gam0) -> np.ndarray:
    """Eigenvalues of the LO singlet anomalous dimensions matrix.

    Args:
          gam0: matrix of LO anomalous dimensions

    Returns:
        lam[a, k]
        a in [+, -] and k is MB contour point index

    """
    # To avoid crossing of the square root cut on the
    # negative real axis we use trick by Dieter Mueller
    aux = ((gam0[..., 0, 0] - gam0[..., 1, 1]) *
           np.sqrt(1. + 4.0 * gam0[..., 0, 1] * gam0[..., 1, 0] /
                   (gam0[..., 0, 0] - gam0[..., 1, 1])**2))
    lam1 = 0.5 * (gam0[..., 0, 0] + gam0[..., 1, 1] - aux)
    lam2 = lam1 + aux
    return np.stack([lam1, lam2])


def projectors(gam0) -> Tuple[np.ndarray, np.ndarray]:
    """Projectors on evolution quark-gluon singlet eigenaxes.

    Args:
          gam0: LO anomalous dimension

    Returns:
         lam: eigenvalues of LO an. dimm matrix lam[a, k]  # Eq. (123)
          pr: Projector pr[k, a, i, j]  # Eq. (122)
               k is MB contour point index
               a in [+, -]
               i,j in {Q, G}

    """
    lam = lambdaf(gam0)
    den = 1. / (lam[0, ...] - lam[1, ...])

    # P+ and P-
    ssm = gam0 - np.einsum('...,ij->...ij', lam[1, ...], np.identity(2))
    ssp = gam0 - np.einsum('...,ij->...ij', lam[0, ...], np.identity(2))
    prp = np.einsum('...,...ij->...ij', den, ssm)
    prm = np.einsum('...,...ij->...ij', -den, ssp)
    # We insert a-axis before i,j-axes, i.e. on -3rd place
    pr = np.stack([prp, prm], axis=-3)
    return lam, pr


def rnlof(m, j) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return projected NLO mu-independent P(bet*gam0 - gam1)P.

    Args:
          m: instance of the model
          j: MB contour points (overrides m.jpoints)

    Returns:
         Tuple of three arrays, namely:
         lam: eigenvalues of LO an. dimm matrix lam[a, k]  # Eq. (123)
         pr: Projector pr[k, a, i, j]  # Eq. (122)
         r1proj: r1proj[a,b] = sum_ij pr[a,i,j] R1[i,j] pr[b,i,j]  # Eq. (124)
         a,b in {+,-};  i,j in {Q, G}

    """
    # cf. my DIS notes p. 61
    gam0 = adim.singlet_LO(j+1, m.nf).transpose((2, 0, 1))  # LO singlet an. dim.
    gam1 = adim.singlet_NLO(j+1, m.nf).transpose((2, 0, 1))  # NLO singlet an. dim.
    lam, pr = projectors(gam0)
    inv = 1.0 / qcd.beta(0, m.nf)
    # Cf. Eq. (124), but defined with opposite sign
    # and with additional 1/b0, as in gepard-fortran
    r1 = inv * (gam1 - 0.5 * inv * qcd.beta(1, m.nf) * gam0)
    r1proj = np.einsum('kaim,kmn,kbnj->kabij', pr, r1, pr)

    return lam, pr, r1proj


def rnlonsf(m, j, prty) -> np.ndarray:
    """Return NLO mu-independent part of evolution operator.

    Args:
          m: instance of the model
          j: MB contour points (overrides m.jpoints)
          prty: 1 for NS^{+}, -1 for NS^{-}

    Returns:
         r1: beta/gama ratio from Eq. (117)

    """
    # cf. my DIS notes p. 61
    gam0 = adim.non_singlet_LO(j+1, m.nf, prty)   # LO non-singlet an. dim.
    gam1 = adim.non_singlet_NLO(j+1, m.nf, prty)  # NLO non-singlet an. dim.
    inv = 1.0 / qcd.beta(0, m.nf)
    # Cf. Eq. (117), but defined with opposite sign as in gepard-fortran
    r1 = inv * (gam1 - 0.5 * inv * qcd.beta(1, m.nf) * gam0)
    return gam0, r1


def erfunc(m, lamj, lamk, R) -> np.ndarray:
    """Mu-dep. part of NLO evolution operator. Eq. (126)."""
    b0 = qcd.beta(0, m.nf)
    levi_civita = np.array([[0, 1], [-1, 0]])
    bll = np.einsum('...,ij->...ij', lamj[0, ...] - lamk[1, ...], levi_civita)
    bll = b0 * np.ones_like(bll) + bll

    er1 = (np.ones_like(bll) - (1./R)**(bll/b0)) / bll  # Eq. (126)
    er1 = b0 * er1   # as defined in gepard-fortran
    return er1


def erfunc_nd(m, lamj, lamk, R) -> np.ndarray:
    """Mu-dep. part of NLO evolution operator. Eq. (126)."""
    assert lamk.shape == (2,)  # FIX THIS
    b0 = qcd.beta(0, m.nf)
    bll = np.subtract.outer(lamj.transpose(), lamk)
    bll = b0 * np.ones_like(bll) + bll

    er1 = (np.ones_like(bll) - (1./R)**(bll/b0)) / bll  # Eq. (126)
    er1 = b0 * er1   # as defined in gepard-fortran
    return er1


def cb1(m, Q2, zn, zk, NS: bool = False):
    """Non-diagonal part of NLO evol op.

    Args:
          m: instance of the model
          Q2: evolution point
          zn: non-diagonal evolution Mellin-Barnes integration point (array)
          zk: COPE Mellin-Barnes integration point (not array! - FIXME)
          NS: do we want non-singlet?

    Returns:
         B_jk: non-diagonal part of evol. op. from Eq. (140)

    Note:
         It's multiplied by GAMMA(3/2) GAMMA(K+3) / (2^(K+1) GAMMA(K+5/2))
         so it's ready to be combined with diagonally evolved C_K = 1 + ...
         where in the end everything will be multiplied by
        (2^(K+1) GAMMA(K+5/2))  / ( GAMMA(3/2) GAMMA(K+3) )

    """
    asmuf2 = qcd.as2pf(m.p, m.nf, Q2, m.asp[m.p], m.r20)
    asQ02 = qcd.as2pf(m.p, m.nf, m.Q02, m.asp[m.p], m.r20)
    R = asmuf2/asQ02
    b0 = qcd.beta(0, m.nf)
    AAA = (special.S1((zn+zk+2)/2) -
           special.S1((zn-zk-2)/2) +
           2*special.S1(zn-zk-1) -
           special.S1(zn+1))
    # GOD = "G Over D" = g_jk / d_jk
    GOD_11 = 2 * constants.CF * (
            2*AAA + (AAA - special.S1(zn+1))*(zn-zk)*(
                zn+zk + 3)/(zk+1)/(zk+2))
    nzero = np.zeros_like(GOD_11)
    GOD_12 = nzero
    GOD_21 = 2*constants.CF*(zn-zk)*(zn+zk+3)/zn/(zk+1)/(zk+2)
    GOD_22 = 2 * constants.CA * (2*AAA + (AAA - special.S1(
        zn + 1)) * (special.poch(zn, 4) / special.poch(zk, 4) -
                    1) + 2 * (zn-zk) * (
            zn+zk + 3) / special.poch(zk, 4)) * zk / zn
    god = np.array([[GOD_11, GOD_12], [GOD_21, GOD_22]])
    dm_22 = zk/zn
    dm_11 = np.ones_like(dm_22)
    dm_12 = np.zeros_like(dm_22)
    dm = np.array([[dm_11, dm_12], [dm_12, dm_22]])
    fac = (zk+1)*(zk+2)*(2*zn+3)/(zn+1)/(zn+2)/(zn-zk)/(zn+zk+3)
    if NS:
        gamn = adim.non_singlet_LO(zn+1, m.nf)
        gamk = adim.non_singlet_LO(zk+1, m.nf)
        r1 = (1 - (1/R)**((b0 + gamn - gamk)/b0)) / (b0+gamn-gamk)
        cb1 = r1 * (gamn-gamk) * (b0 - gamk + GOD_11) * R**(-gamk/b0)
        cb1 = fac * cb1
    else:
        gamn = adim.singlet_LO(zn+1, m.nf).transpose((2, 0, 1))
        gamk = adim.singlet_LO(zk+1, m.nf)
        lamn, pn = evolution.projectors(gamn)
        lamk, pk = evolution.projectors(gamk)
        proj_DM = np.einsum('naif,fgn,bgj->nabij', pn, dm, pk)
        proj_GOD = np.einsum('naif,fgn,bgj->nabij', pn, god, pk)
        er1 = evolution.erfunc_nd(m, lamn, lamk, R)
        bet_proj_DM = np.einsum('b,nabij->nabij', b0-lamk, proj_DM)
        cb1 = np.einsum('n,nab,nabij,b->nij', fac,
                        er1*np.subtract.outer(lamn.transpose(), lamk),
                        bet_proj_DM + proj_GOD,
                        R**(-lamk/b0)) / b0
    return cb1


def evolop(m, j, Q2: float, process_class: str) -> np.ndarray:
    """GPD evolution operator.

    Args:
         m: instance of the model
         j: MB contour points (overrides m.jpoints)
         Q2: final evolution momentum squared
         process_class: DIS, DVCS or DVMP

    Returns:
         Array corresponding Eq. (121) of Towards DVCS paper.
         evolop[k, p, i, j]
         -  k is index of point on MB contour,
         -  p is pQCD order (0=LO, 1=NLO)
         -  i, j in [Q, G]

    Todo:
        Argument should not be a process class but GPD vs PDF, or we
        should avoid it altogether somehow. This serves here only
        to get the correct choice of evolution scheme (msbar vs csbar).

    """
    # 1. Alpha-strong ratio.
    # When m.p=1 (NLO), LO part of the evolution operator
    # will still be multiplied by ratio of alpha_strongs
    # evaluated at NLO, as it should.
    asmuf2 = qcd.as2pf(m.p, m.nf, Q2, m.asp[m.p], m.r20)
    asQ02 = qcd.as2pf(m.p, m.nf, m.Q02, m.asp[m.p], m.r20)
    R = asmuf2/asQ02

    # 2. egeinvalues, projectors, projected mu-indep. part
    lam, pr, r1proj = rnlof(m, j)

    # 3. LO errfunc
    b0 = qcd.beta(0, m.nf)
    er1 = erfunc(m, lam, lam, R)

    Rfact = R**(-lam/b0)  # LO evolution (alpha(mu)/alpha(mu0))^(-gamma/beta0)
    evola0ab = np.einsum('kaij,ab->kabij', pr,  np.identity(2))
    evola0 = np.einsum('kabij,bk->kij', evola0ab, Rfact)

    if m.p == 1:
        # Cf. eq. (124), but with opposite sign, as in gepard-fortran
        evola1ab = - np.einsum('kab,kabij->kabij', er1, r1proj)
        evola1 = np.einsum('kabij,bk->kij', evola1ab, Rfact)
        # adding non-diagonal evolution when needed or asked for
        # if ((process_class == 'DVMP') or (
        #        process_class == 'DVCS' and m.scheme == 'msbar')):
        if ((process_class != 'DIS') and (m.scheme == 'msbar')):
            # FIXME: find a way to do it array-wise i.e. get rid of j_single
            nd = []
            for j_single in j:
                znd, wgnd = quadrature.nd_mellin_barnes()
                ndphij = 1.57j
                ephnd = np.exp(ndphij)
                tginv = ephnd/np.tan(np.pi*znd/2)
                tginvc = ephnd.conjugate()/np.tan(np.pi*znd.conjugate()/2)
                cb1f = cb1(m, Q2, j_single + znd + 2, j_single)
                cb1fc = cb1(m, Q2, j_single + znd.conjugate() + 2, j_single)
                ndint = np.einsum('n,nij,n->ij', wgnd, cb1f, tginv)
                ndint -= np.einsum('n,nij,n->ij', wgnd, cb1fc, tginvc)
                ndint = ndint * 0.25j
                nd.append(ndint)
            evola1 += np.array(nd)
    else:
        evola1 = np.zeros_like(evola0)

    evola = np.stack((evola0, evola1), axis=1)

    return evola


def evolopns(m, j, Q2: float, process_class: str) -> np.ndarray:
    """GPD evolution operator (NSP case only atm).

    Args:
         m: instance of the model
         j: MB contour points (overrides m.jpoints)
         Q2: final evolution momentum squared
         process_class: DIS, DVCS or DVMP

    Returns:
         Array corresponding Eq. (116) of Towards DVCS paper.
         evolopns[k, p]
         -  k is index of point on MB contour,
         -  p is pQCD order (0=LO, 1=NLO)

    Notes:
        Code duplication, should be merged with evolop function

    """
    # 1. Alpha-strong ratio.
    # When m.p=1 (NLO), LO part of the evolution operator
    # will still be multiplied by ratio of alpha_strongs
    # evaluated at NLO, as it should.
    asmuf2 = qcd.as2pf(m.p, m.nf, Q2, m.asp[m.p], m.r20)
    asQ02 = qcd.as2pf(m.p, m.nf, m.Q02, m.asp[m.p], m.r20)
    R = asmuf2/asQ02

    # 2. mu-indep. part
    gam0, r1 = rnlonsf(m, j, 1)   # prty=1 fixed

    # 2. LO errfunc
    b0 = qcd.beta(0, m.nf)
    aux1 = - (1 - 1 / R) * r1  # Cf. eq. (117), but with opposite sign
    aux0 = np.ones_like(aux1)

    evola0 = aux0 * R**(-gam0/b0)  # LO evolution (alpha(mu)/alpha(mu0))^(-gamma/beta0)

    if m.p == 1:
        evola1 = aux1 * R**(-gam0/b0)
        # adding non-diagonal evolution when needed or asked for
        # if ((process_class == 'DVMP') or (
        #    process_class == 'DVCS' and m.scheme == 'msbar')):
        if ((process_class != 'DIS') and (m.scheme == 'msbar')):
            # FIXME: find a way to do it array-wise i.e. get rid of j_single
            nd = []
            for j_single in j:
                znd, wgnd = quadrature.nd_mellin_barnes()
                ndphij = 1.57j
                ephnd = np.exp(ndphij)
                tginv = ephnd/np.tan(np.pi*znd/2)
                tginvc = ephnd.conjugate()/np.tan(np.pi*znd.conjugate()/2)
                cb1f = cb1(m, Q2, j_single + znd + 2, j_single, NS=True)
                cb1fc = cb1(m, Q2, j_single + znd.conjugate() + 2, j_single, NS=True)
                ndint = np.sum(wgnd * cb1f * tginv)
                ndint -= np.sum(wgnd * cb1fc * tginvc)
                ndint = ndint * 0.25j
                nd.append(ndint)
            evola1 += np.array(nd)
    else:
        evola1 = np.zeros_like(evola0)

    evola = np.stack((evola0, evola1), axis=1)

    return evola
