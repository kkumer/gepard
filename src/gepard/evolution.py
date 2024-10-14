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

import copy
import numpy as np

from . import adim, constants, evolution, qcd, quadrature, special


def lambdaf(gam0: np.ndarray) -> np.ndarray:
    """Return eigenvalues of singlet LO anomalous dimensions matrix.

    Args:
          gam0: matrix of LO anomalous dimensions of shape [..., 2, 2]

    Returns:
        lam[..., a],  a in [+, -] where
          + and - are eigenvalues in singlet quark - gluon space

    """
    # To avoid crossing of the square root cut on the
    # negative real axis we use trick by Dieter Mueller
    aux = ((gam0[..., 0, 0] - gam0[..., 1, 1]) *
           np.sqrt(1. + 4.0 * gam0[..., 0, 1] * gam0[..., 1, 0] /
                   (gam0[..., 0, 0] - gam0[..., 1, 1])**2))
    lam1 = 0.5 * (gam0[..., 0, 0] + gam0[..., 1, 1] - aux)
    lam2 = lam1 + aux
    return np.stack([lam1, lam2], axis=-1)


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
    den = 1. / (lam[..., 0] - lam[..., 1])

    # P+ and P-
    ssm = gam0 - np.einsum('...,ij->...ij', lam[..., 1], np.identity(2))
    ssp = gam0 - np.einsum('...,ij->...ij', lam[..., 0], np.identity(2))
    prp = np.einsum('...,...ij->...ij', den, ssm)
    prm = np.einsum('...,...ij->...ij', -den, ssp)
    # We insert a-axis before i,j-axes, i.e. on -3rd place
    pr = np.stack([prp, prm], axis=-3)

    return lam, pr


def rnlof(m) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return projected singlet NLO mu-independent P(bet*gam0 - gam1)P.

    Args:
          m: instance of the model

    Returns:
         Tuple of three arrays, namely:
         lam: eigenvalues of LO an. dimm matrix lam[a, k]  # Eq. (123)
         pr: Projector pr[k, a, i, j]  # Eq. (122)
         r1proj: r1proj[a,b] = sum_ij pr[a,i,j] R1[i,j] pr[b,i,j]  # Eq. (124)
         a,b in {+,-};  i,j in {Q, G}

    """
    gam0 = adim.singlet_LO(m.jpoints_pws + 1, m.nf)
    gam1 = adim.singlet_NLO(m.jpoints_pws + 1, m.nf)
    lam, pr = projectors(gam0)
    inv = 1.0 / qcd.beta(0, m.nf)
    # Cf. Eq. (124), but defined with opposite sign
    # and with additional 1/b0, as in gepard-fortran
    r1 = inv * (gam1 - 0.5 * inv * qcd.beta(1, m.nf) * gam0)
    r1proj = np.einsum('skaim,skmn,skbnj->skabij', pr, r1, pr)

    return lam, pr, r1proj


def rnlonsf(m, prty: int) -> np.ndarray:
    """Return NLO mu-independent part of evolution operator.

    Args:
          m: instance of the model
          j: MB contour points (overrides m.jpoints)
          prty: 1 for NS^{+}, -1 for NS^{-}
                                                                                                          
    Returns:                                                                                              
         r1: second mu-independent factor from Eq. (117)
                                                                                                          
    """                                                                                                   
    gam0 = adim.non_singlet_LO(m.jpoints_pws + 1, m.nf)
    gam1 = adim.non_singlet_NLO(m.jpoints_pws + 1, m.nf, prty)
    inv = 1.0 / qcd.beta(0, m.nf)
    # Cf. Eq. (117), but defined with opposite sign as in gepard-fortran
    r1 = inv * (gam1 - 0.5 * inv * qcd.beta(1, m.nf) * gam0)
    return gam0, r1

def erfunc(m, lamj, lamk, R) -> np.ndarray:
    """Erfunc from Eq. (126)."""
    b0 = qcd.beta(0, m.nf)
    lamj_ex = lamj[..., :, np.newaxis]
    lamk_ex = lamk[..., np.newaxis, :]
    bllpre = lamj_ex - lamk_ex
    bll = b0 * np.ones_like(bllpre) + bllpre
    #
    er1 = (np.ones_like(bll) - (1./R)**(bll/b0)) / bll  # Eq. (126)
    er1 = b0 * er1   # as defined in gepard-fortran
    return er1

def cb1(m, R, zn, zk):
    """Non-diagonal part of singlet NLO evol op.

    Args:
          m: instance of the model
          zn: non-diagonal evolution Mellin-Barnes integration point (array)
          zk: COPE Mellin-Barnes integration point (not array! - FIXME)

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
    god = np.array([[GOD_11, GOD_12], 
                    [GOD_21, GOD_22]])
    dm_22 = zk/zn
    dmzero = np.zeros_like(dm_22)
    dm_11 = np.ones_like(dm_22)
    dm_12 = dmzero
    dm = np.array([[dm_11, dm_12],
                   [dm_12, dm_22]])
    gamn = adim.singlet_LO(zn+1, m.nf)
    gamk = adim.singlet_LO(zk+1, m.nf)
    lamn, pn = projectors(gamn)
    lamk, pk = projectors(gamk)

    proj_DM = np.einsum('sknaif,fgskn,skbgj->sknabij', pn, dm, pk[:, :, 0,...])
    proj_GOD = np.einsum('sknaif,fgskn,skbgj->sknabij', pn, god, pk[:, :, 0,...])

    er1 = erfunc(m, lamn, lamk, R)
    bet_proj_DM = np.einsum('skb,sknabij->sknabij', b0-lamk[...,0, :], proj_DM)
    lamn_ex = lamn[..., np.newaxis]
    lamk_ex = lamk[:, :, np.newaxis, 0, np.newaxis, :]
    bllpre = lamn_ex - lamk_ex
    # bllpre = np.transpose(bllpre, (3, 0, 1, 2))
    cb1 = np.einsum('skn,sknab,sknabij,skb->sknij', fac,
                    er1*bllpre,
                    bet_proj_DM + proj_GOD,
                    R**(-lamk[:,:,0,:]/b0)) / b0
    return cb1


def cb1ns(m, R, zn, zk):
    """Non-diagonal part of non-singlet NLO evol op.

    Args:
          m: instance of the model
          zn: non-diagonal evolution Mellin-Barnes integration point (array)
          zk: COPE Mellin-Barnes integration point (not array! - FIXME)

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

    gamn = adim.non_singlet_LO(zn+1, m.nf)                                                            
    gamk = adim.non_singlet_LO(zk+1, m.nf)                                                            

    r1 = (1 - (1/R)**((b0 + gamn - gamk)/b0)) / (b0+gamn-gamk)                                        
    cb1 = r1 * (gamn-gamk) * (b0 - gamk + GOD_11) * R**(-gamk/b0)                                     
    cb1 = fac * cb1          

    return cb1



def evolop(m, R: float, process_class: str) -> np.ndarray:
    """GPD evolution operator.

    Args:
         m: instance of the model
         R: ratio of astrong(mu_fact)/astrong(mu_input)

    Returns:
         Array corresponding Eq. (121) of Towards DVCS paper.
         evolop[s, k, p, i, j]
         -  s is index of SO(3) partial wave
         -  k is index of point on MB contour,
         -  p is pQCD order (0=LO, 1=NLO)
         -  i, j in [Q, G, NSP, ...] i.e. in m.evolution_basis

    Todo:
        Argument should not be a process class but GPD vs PDF, or we
        should avoid it altogether somehow. This serves here only
        to get the correct choice of evolution scheme (msbar vs csbar).

    """

    ndim = len(m.evolution_basis)
    evola = np.zeros((m.npws, m.npts, m.p + 1, ndim, ndim), dtype=complex)
    b0 = qcd.beta(0, m.nf)

    rest_ebas = 0  # remaining items in evolution basis
    if (ndim > 1) and (m.evolution_basis[1] == 'G'):
        # we need singlet flavor-mixing part first
        # eigenvalues, projectors, projected mu-indep. part
        lam, pr, r1proj = rnlof(m)
        er1 = erfunc(m, lam, lam, R)

        Rfact = R**(-lam/b0)  # LO evolution (alpha(mu)/alpha(mu0))^(-gamma/beta0)
        evola0ab = np.einsum('skaij,ab->skabij', pr,  np.identity(2))
        evola[:, :, 0, :2, :2] = np.einsum('skabij,skb->skij', evola0ab, Rfact)
        if m.p == 1:
            # Cf. eq. (124), but with opposite sign, as in gepard-fortran
            evola1ab = - np.einsum('skab,skabij->skabij', er1, r1proj)
            evola[:, :, 1, :2, :2] = np.einsum('skabij,skb->skij', evola1ab, Rfact)
            # adding non-diagonal evolution when needed or asked for
            # if ((process_class == 'DVMP') or (
            #        process_class == 'DVCS' and m.scheme == 'msbar')):
            if ((process_class != 'DIS') and (m.scheme == 'msbar')):
                zj = m.jpoints_pws[:, :, np.newaxis]
                znd, wgnd = quadrature.nd_mellin_barnes()
                znd = znd[np.newaxis, :]
                ndphij = 1.57j
                ephnd = np.exp(ndphij)
                tginv = ephnd/np.tan(np.pi*znd/2)
                tginvc = ephnd.conjugate()/np.tan(np.pi*znd.conjugate()/2)
                cb1f = cb1(m, R, zj + znd + 2, zj)
                cb1fc = cb1(m, R, zj + znd.conjugate() + 2, zj)
                if process_class == 'DVMP':
                    Vdvmp = np.zeros_like(cb1f)
                    Vdvmp[..., 0, 0] = np.ones_like(zj+znd)
                    Vdvmp[..., 1, 1] =  (zj+3)/(zj+znd+5)
                    cb1f = np.einsum('sknab,sknbc->sknac', Vdvmp, cb1f)
                    # Vdvmpc = copy.deepcopy(Vdvmp)  # maybe faster?
                    # Vdvmpc[..., 1, 1] =  (zj+3)/(zj+znd.conjugate()+5)
                    Vdvmpc = np.zeros_like(cb1fc)
                    Vdvmpc[..., 0, 0] = np.ones_like(zj+znd)
                    Vdvmpc[..., 1, 1] =  (zj+3)/(zj+znd.conjugate()+5)
                    cb1fc = np.einsum('sknab,sknbc->sknac', Vdvmpc, cb1fc)
                ndint = np.einsum('n,sknij,kn->skij', wgnd, cb1f, tginv)
                ndint -= np.einsum('n,sknij,kn->skij', wgnd, cb1fc, tginvc)
                ndint = ndint * 0.25j
                evola[:, :, 1, :2, :2] += ndint
        rest_ebas = 2  # we processed 2 items of evolution basis
    rest_of_basis = m.evolution_basis[rest_ebas:]  # may be empty if only SI neded
    for count, item in enumerate(rest_of_basis):
        prty = {'NSP':1, 'NSM':-1}[item]
        gam0, r1 = rnlonsf(m, prty)
        aux0 = np.ones_like(r1)
        evola[:, :, 0, rest_ebas+count, rest_ebas+count] = aux0 * R**(-gam0/b0)  
        if m.p == 1:
            aux1 = - (1 - 1 / R) * r1  # Cf. eq. (117), but with opposite sign
            evola[:, :, 1, rest_ebas+count, rest_ebas+count] = aux1 * R**(-gam0/b0)
            if ((process_class != 'DIS') and (m.scheme == 'msbar')):
                zj = m.jpoints_pws[:, :, np.newaxis]
                znd, wgnd = quadrature.nd_mellin_barnes()
                znd = znd[np.newaxis, :]
                ndphij = 1.57j
                ephnd = np.exp(ndphij)
                tginv = ephnd/np.tan(np.pi*znd/2)
                tginvc = ephnd.conjugate()/np.tan(np.pi*znd.conjugate()/2)
                cb1f = cb1ns(m, R, zj + znd + 2, zj)
                cb1fc = cb1ns(m, R, zj + znd.conjugate() + 2, zj)
                ndint = np.einsum('n,skn,kn->sk', wgnd, cb1f, tginv)
                ndint -= np.einsum('n,skn,kn->sk', wgnd, cb1fc, tginvc)
                ndint = ndint * 0.25j
                evola[:, :, 1, rest_ebas+count, rest_ebas+count] += ndint

    return evola


def evolop_da(m, R: float, prty=+1) -> np.ndarray:
    """DA evolution operator. Non-singlet only

    Args:
             m: instance of the model
       gpoints: array of Gegenbauer indices
             R: ratio of astrong(mu_fact_da)/astrong(mu_input)
          prty: C parity of DA

    Returns:
         evola_da[g, p]
         -  g is DA Gegenbauer index in m.gpoints
         -  p is pQCD order (0=LO, 1=NLO)

    Todo:
        Merge this with evolop for GPDs. Ugly code duplication of the most
        complex part of the evolution code.

    """

    evola = np.zeros((m.ngegens, m.p + 1), dtype=complex)

    b0 = qcd.beta(0, m.nf)
    gam0 = adim.non_singlet_LO(m.gpoints + 1, m.nf)
    inv = 1.0 / qcd.beta(0, m.nf)
    
    evola[:, 0] = R**(-gam0/b0)  
    if m.p == 1:
        gam1 = adim.non_singlet_NLO(m.gpoints + 1, m.nf, prty)
        r1 = inv * (gam1 - 0.5 * inv * qcd.beta(1, m.nf) * gam0)
        aux1 = - (1 - 1 / R) * r1  # Cf. eq. (117), but with opposite sign
        evola[:, 1] = aux1 * R**(-gam0/b0)
        if m.scheme == 'msbar':
            zj = m.gpoints[:, np.newaxis]
            znd, wgnd = quadrature.nd_mellin_barnes()
            znd = znd[np.newaxis, :]
            ndphij = 1.57j
            ephnd = np.exp(ndphij)
            tginv = ephnd/np.tan(np.pi*znd/2)
            tginvc = ephnd.conjugate()/np.tan(np.pi*znd.conjugate()/2)
            cb1f = cb1ns(m, R, zj + znd + 2, zj)
            cb1fc = cb1ns(m, R, zj + znd.conjugate() + 2, zj)
            ndint = np.einsum('n,kn,kn->k', wgnd, cb1f, tginv)
            ndint -= np.einsum('n,kn,kn->k', wgnd, cb1fc, tginvc)
            ndint = ndint * 0.25j
            evola[:, 1] += ndint

    return evola
