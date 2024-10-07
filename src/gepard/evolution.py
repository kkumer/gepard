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
    """Return SI-NS block LO anomalous dimensions matrix.

    Args:
          gam0: matrix of LO anomalous dimensions

    Returns:
        lam[a, k]
        a in [+, -, NS+, NS-] and k is MB contour point index
           + and - are eigenvalues in singlet quark - gluon space

    """
    # To avoid crossing of the square root cut on the
    # negative real axis we use trick by Dieter Mueller
    aux = ((gam0[..., 0, 0] - gam0[..., 1, 1]) *
           np.sqrt(1. + 4.0 * gam0[..., 0, 1] * gam0[..., 1, 0] /
                   (gam0[..., 0, 0] - gam0[..., 1, 1])**2))
    lam1 = 0.5 * (gam0[..., 0, 0] + gam0[..., 1, 1] - aux)
    lam2 = lam1 + aux
    return np.stack([lam1, lam2, gam0[..., 0, 0], gam0[..., 0, 0]])


def projectors(gam0):
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
    ssm = gam0 - np.einsum('...,ij->...ij', lam[1, ...], np.identity(4))
    ssp = gam0 - np.einsum('...,ij->...ij', lam[0, ...], np.identity(4))
    prp = np.einsum('...,...ij->...ij', den, ssm)
    prm = np.einsum('...,...ij->...ij', -den, ssp)
    # We insert a-axis before i,j-axes, i.e. on -3rd place
    pr = np.stack([prp, prm], axis=-3)

    # Extend this to NS+, NS- extra two dimensions
    # Initialize the new array with shape (..., 4, 4, 4)
    new_pr = np.zeros((pr.shape[:-3] + (4, 4, 4)), dtype=pr.dtype)
    
    # Copy the existing projectors into the new array
    new_pr[..., 0, :2, :2] = pr[..., 0, :2, :2]  
    new_pr[..., 1, :2, :2] = pr[..., 1, :2, :2]
    
    # Set the trivial non-singlet projectors
    new_pr[..., 2, :, :] = np.array([[0, 0, 0, 0],
                                   [0, 0, 0, 0],
                                   [0, 0, 1, 0], 
                                   [0, 0, 0, 0]])
    new_pr[..., 3, :, :] = np.array([[0, 0, 0, 0],
                                   [0, 0, 0, 0],
                                   [0, 0, 0, 0], 
                                   [0, 0, 0, 1]])

    return lam, new_pr


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
    gam = adim.block(j+1, m.nf)  #  gam[k, p, i, j]
    gam[np.isnan(gam)] = 1.  # FIXME: cludge to regulate for integer j (DA evolution)
    lam, pr = projectors(gam[:, 0, :, :])
    inv = 1.0 / qcd.beta(0, m.nf)
    # Cf. Eq. (124), but defined with opposite sign
    # and with additional 1/b0, as in gepard-fortran
    r1 = inv * (gam[:, 1, :, :] - 0.5 * inv * qcd.beta(1, m.nf) * gam[:, 0, :, :])
    r1proj = np.einsum('kaim,kmn,kbnj->kabij', pr, r1, pr)

    return lam, pr, r1proj


def erfunc(m, lamj, lamk, R) -> np.ndarray:
    """Mu-dep. part of NLO evolution operator. Eq. (126)."""
    b0 = qcd.beta(0, m.nf)
    levi_civita = np.array([[0, 1], [-1, 0]])
    bll = np.einsum('...,ij->...ij', lamj[0, ...] - lamk[1, ...], levi_civita)
    bll = b0 * np.ones_like(bll) + bll

    er1 = (np.ones_like(bll) - (1./R)**(bll/b0)) / bll  # Eq. (126)
    er1 = b0 * er1   # as defined in gepard-fortran
    return er1


def erfunc_di(m, lamj, lamk, R) -> np.ndarray:
    """Mu-dep. part of NLO evolution operator. Eq. (126)."""
    b0 = qcd.beta(0, m.nf)
    bll = np.subtract.outer(lamj.transpose(), lamk)
    bll = b0 * np.ones_like(bll) + bll
    er1 = (np.ones_like(bll) - (1./R)**(bll/b0)) / bll  # Eq. (126)
    er1 = b0 * er1   # as defined in gepard-fortran
    return er1


def erfunc_nd(m, lamj, lamk, R) -> np.ndarray:
    """Mu-dep. part of NLO evolution operator. Eq. (126)."""
    b0 = qcd.beta(0, m.nf)
    lamj_ex = lamj[:, np.newaxis, :, :]
    lamk_ex = lamk[:, :, 0][np.newaxis, :, :, np.newaxis]
    bllpre = lamj_ex - lamk_ex
    bllpre = np.transpose(bllpre, (3, 0, 1, 2))
    bll = b0 * np.ones_like(bllpre) + bllpre
    #
    er1 = (np.ones_like(bll) - (1./R)**(bll/b0)) / bll  # Eq. (126)
    er1 = b0 * er1   # as defined in gepard-fortran
    return er1


def cb1(m, R, zn, zk, NS: bool = False):
    """Non-diagonal part of NLO evol op.

    Args:
          m: instance of the model
          zn: non-diagonal evolution Mellin-Barnes integration point (array)
          zk: COPE Mellin-Barnes integration point (not array! - FIXME)
          NS: do we want non-singlet?

    Returns:
         B_jk: non-diagonal part of evol. op. from Eq. (140).
         For NS=True we get array cb1[zn=j] and dependent on zk=k.
         For NS=False we get array cb1[zn=j, i, j] and dependent on zk=k,
         where i,j are in [Q, G].

    Note:
         It's multiplied by GAMMA(3/2) GAMMA(K+3) / (2^(K+1) GAMMA(K+5/2))
         so it's ready to be combined with diagonally evolved C_K = 1 + ...
         where in the end everything will be multiplied by
        (2^(K+1) GAMMA(K+5/2))  / ( GAMMA(3/2) GAMMA(K+3) )

    """
    b0 = qcd.beta(0, m.nf)
    AAA = (special.S1((zn+zk+2)/2) -
           special.S1((zn-zk-2)/2) +
           2*special.S1(zn-zk-1) -
           special.S1(zn+1))
    # GOD = "G Over D" = g_jk / d_jk
    GOD_11 = 2 * constants.CF * (
            2*AAA + (AAA - special.S1(zn+1))*(zn-zk)*(
                zn+zk + 3)/(zk+1)/(zk+2))
    fac = (zk+1)*(zk+2)*(2*zn+3)/(zn+1)/(zn+2)/(zn-zk)/(zn+zk+3)
    nzero = np.zeros_like(GOD_11)
    GOD_12 = nzero
    GOD_21 = 2*constants.CF*(zn-zk)*(zn+zk+3)/zn/(zk+1)/(zk+2)
    GOD_22 = 2 * constants.CA * (2*AAA + (AAA - special.S1(
        zn + 1)) * (special.poch(zn, 4) / special.poch(zk, 4) -
                    1) + 2 * (zn-zk) * (
            zn+zk + 3) / special.poch(zk, 4)) * zk / zn
    god = np.array([[GOD_11, GOD_12, nzero, nzero], 
                    [GOD_21, GOD_22, nzero, nzero],
                    [nzero, nzero, GOD_11, nzero],
                    [nzero, nzero, nzero, GOD_11]])
    dm_22 = zk/zn
    dmzero = np.zeros_like(dm_22)
    dm_11 = np.ones_like(dm_22)
    dm_12 = dmzero
    dm = np.array([[dm_11, dm_12, dmzero, dmzero],
                   [dm_12, dm_22, dmzero, dmzero],
                   [dmzero, dmzero, dm_11, dmzero],
                   [dmzero, dmzero, dmzero, dm_11]])
    # gamn = adim.singlet_LO(zn+1, m.nf).transpose((2, 0, 1))
    # gamk = adim.singlet_LO(zk+1, m.nf)
    gamn = adim.block(zn+1, m.nf)
    gamk = adim.block(zk+1, m.nf)
    lamn, pn = projectors(gamn[..., 0, :, :])
    lamk, pk = projectors(gamk[..., 0, :, :])
    #  np.einsum('naif,fgn,bgj->nabij', pn, dm, pk)
    proj_DM = np.einsum('knaif,fgkn,kbgj->knabij', pn, dm, pk[:, 0,...])
    proj_GOD = np.einsum('knaif,fgkn,kbgj->knabij', pn, god, pk[:, 0,...])
    er1 = erfunc_nd(m, lamn, lamk, R)
    bet_proj_DM = np.einsum('bk,knabij->knabij', b0-lamk[...,0], proj_DM)
    lamn_ex = lamn[:, np.newaxis, :, :]
    lamk_ex = lamk[:, :, 0][np.newaxis, :, :, np.newaxis]
    bllpre = lamn_ex - lamk_ex
    bllpre = np.transpose(bllpre, (3, 0, 1, 2))
    cb1 = np.einsum('kn,nabk,knabij,bk->nkij', fac,
                    er1*bllpre,
                    bet_proj_DM + proj_GOD,
                    R**(-lamk[...,0]/b0)) / b0
    return cb1


def evolop(m, j, R: float, process_class: str, DA=False) -> np.ndarray:
    """GPD evolution operator.

    Args:
         m: instance of the model
         j: MB contour points (overrides m.jpoints)
         R: ratio of astrong(mu_fact)/astrong(mu_input)
        Q2: default endpoint of evolution (can be changed with m.rf2)
         process_class: DIS, DVCS or DVMP
        DA: If DA is evolved (NS only)

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
    # eigenvalues, projectors, projected mu-indep. part
    kronecker = np.identity(4)
    lam, pr, r1proj = rnlof(m, j)
    if DA:
        kronecker = np.identity(2)
        lam = lam[2:, :]
        pr = pr[..., 2:, 2:, 2:]
        r1proj = r1proj[..., 2:, 2:, 2:, 2:]
    # LO errfunc
    b0 = qcd.beta(0, m.nf)
    er1 = erfunc_di(m, lam, lam, R)
    # Next line is needed for integer j (DA evolution)
    # FIXME: Is it too expensive? Should not be needed n>=2, isnt it?
    # lam[np.isnan(lam)] = 1.
    # r1proj[np.isnan(r1proj)] = 1.

    Rfact = R**(-lam/b0)  # LO evolution (alpha(mu)/alpha(mu0))^(-gamma/beta0)
    evola0ab = np.einsum('kaij,ab->kabij', pr,  kronecker)
    evola0 = np.einsum('kabij,bk->kij', evola0ab, Rfact)
    if m.p == 1:
        # Cf. eq. (124), but with opposite sign, as in gepard-fortran
        evola1ab = - np.einsum('kabk,kabij->kabij', er1, r1proj)
        evola1 = np.einsum('kabij,bk->kij', evola1ab, Rfact)
        # adding non-diagonal evolution when needed or asked for
        # if ((process_class == 'DVMP') or (
        #        process_class == 'DVCS' and m.scheme == 'msbar')):
        if ((process_class != 'DIS') and (m.scheme == 'msbar') and not DA):
            zj = j[:, np.newaxis]
            znd, wgnd = quadrature.nd_mellin_barnes()
            znd = znd[np.newaxis, :]
            ndphij = 1.57j
            ephnd = np.exp(ndphij)
            tginv = ephnd/np.tan(np.pi*znd/2)
            tginvc = ephnd.conjugate()/np.tan(np.pi*znd.conjugate()/2)
            cb1f = cb1(m, R, zj + znd + 2, zj)
            cb1fc = cb1(m, R, zj + znd.conjugate() + 2, zj)
            ndint = np.einsum('n,nkij,kn->kij', wgnd, cb1f, tginv)
            ndint -= np.einsum('n,nkij,kn->kij', wgnd, cb1fc, tginvc)
            ndint = ndint * 0.25j
            evola1 += ndint
    else:
        evola1 = np.zeros_like(evola0)

    evola = np.stack((evola0, evola1), axis=1)

    return evola


def evolop_old(m, j, R: float, process_class: str) -> np.ndarray:
    """GPD evolution operator.

    Args:
         m: instance of the model
         j: MB contour points (overrides m.jpoints)
         R: ratio of astrong(mu_fact)/astrong(mu_input)
         Q2: default endpoint of evolution (can be changed with m.rf2)
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
    # eigenvalues, projectors, projected mu-indep. part
    lam, pr, r1proj = rnlof(m, j)
    # LO errfunc
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
                cb1f = cb1(m, R, j_single + znd + 2, j_single)
                cb1fc = cb1(m, R, j_single + znd.conjugate() + 2, j_single)
                if process_class == 'DVMP':
                    Vdvmp = np.array([[np.ones_like(znd), np.zeros_like(znd)],
                                        [np.zeros_like(znd),
                                           (j_single+3)/(j_single+znd+5)]])
                    cb1f = np.einsum('ikn,nkj->nij', Vdvmp, cb1f)
                    Vdvmpc = np.array([[np.ones_like(znd), np.zeros_like(znd)],
                                        [np.zeros_like(znd),
                                           (j_single+3)/(j_single+znd.conjugate()+5)]])
                    cb1fc = np.einsum('ikn,nkj->nij', Vdvmpc, cb1fc)
                ndint = np.einsum('n,nij,n->ij', wgnd, cb1f, tginv)
                ndint -= np.einsum('n,nij,n->ij', wgnd, cb1fc, tginvc)
                ndint = ndint * 0.25j
                nd.append(ndint)
            evola1 += np.array(nd)
    else:
        evola1 = np.zeros_like(evola0)

    evola = np.stack((evola0, evola1), axis=1)

    return evola
