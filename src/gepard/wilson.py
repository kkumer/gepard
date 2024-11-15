"""Wilson coefficients and evolved Wilson coefficients."""

import numpy as np
from scipy.special import loggamma  # type: ignore

from . import adim, c1dvcs, c1dvmp, constants, evolution, qcd, theory, gegenbauer


def _fshu(j: np.ndarray) -> np.ndarray:
    """Shuvaev factor."""
    #  = 2^(J+1) Gamma(5/2+J) / Gamma(3/2) / Gamma(3+J) =
    #  = 2^N Gamma(3/2+N) / Gamma(3/2) / Gamma(2+N)
    return (2**(j+1) * np.exp(loggamma(2.5 + j)
                              - loggamma(3 + j) - loggamma(3/2)))

def calc_wc(m: theory.Theory, process_class: str):
    """Calculate Wilson coeffs.

    Args:
        m: instance of the Theory class
        process_class: 'DIS', 'DVCS' or 'DVMP'

    Returns:
        wc[s, k, p, a]: s in range(npws), k in range(npts), p in [LO, NLO],
                        a in [Q, G, NSP, NSM]

    """
    # Initialize array of correct shape with complex zeros
    nebas = len(m.evolution_basis)
    wc = np.zeros((m.npws, m.npts, m.p + 1, nebas), dtype=complex)
    # LO hard scattering coefs are equal to 1 for all flavors ...
    wc[:, :, 0, :] = np.ones_like(wc[:, :, 0, :])
    # ... apart from gluons which are now set to zero:
    if 'G' in m.evolution_basis:
        wc[:, :, 0, m.evolution_basis.index('G')] = np.zeros_like(wc[:, :, 0, 0])
    # Calculate NLO coefs if needed:
    if m.p == 1:
        # FIXME: C1 should calculate what's needed and we should not cut here
        wc[:, :, 1, :] = c1dvcs.C1(m, process_class)[..., :nebas]
    if process_class == 'DIS':
        prefac = np.ones_like(m.jpoints_pws)
    elif process_class == 'DVCS':
        prefac = _fshu(m.jpoints_pws)
    else:
        raise ValueError('process_class {} is not "DIS" or "DVCS"!'.format(process_class))
    wc = np.einsum('sk,skpa->skpa', prefac, wc)
    return wc


def calc_wc_dvmp(m: theory.Theory, g: int):
    """Calculate Wilson coeffs for DVMP.

    Args:
        m: instance of the Theory class
        g: index of DA Gegenbauer moment (often k in published papers)

    Returns:
        wc[s, k, p, f]: s in range (npws), k in range(npts), p in [LO, NLO], f in [Q,G,NSP]

    """
    nebas = len(m.evolution_basis)  # GPD, not DA
    wc = np.zeros((m.npws, m.npts, m.p + 1, nebas), dtype=complex)
    j = m.jpoints_pws
    fshu = _fshu(j)
    # Normalizations. Factor 3 is from normalization of DA, so not NC
    # See p. 37, 39 of "Towards DVMP" paper. Eq. (3.62c)
    # FIXME: NSP normalizations for DVMP not checked yet
    quark_norm = 3 * fshu
    gluon_norm = 3 * fshu * 2 / constants.CF / (j + 3)
    # Normalizations are chosen so that LO is normalized to:
    wc[:, :, 0, 0] = quark_norm * np.ones_like(j) / m.nf   # Q
    wc[:, :, 0, 1] = gluon_norm * np.ones_like(j)          # G
    wc[:, :, 0, 2] = quark_norm * np.ones_like(j)          # NSP ?!
    if m.p == 1:
        qp1, ps1, g1 = c1dvmp.c1dvmp(m, 1, j, g)
        q1 = qp1/m.nf + ps1
        nsp1 = qp1
        wc[:, :, 1, 0] = quark_norm * q1  # Q
        wc[:, :, 1, 1] = gluon_norm * g1  # Q
        wc[:, :, 1, 2] = quark_norm * nsp1 # NSP ?!
    return wc


def calc_pj(m: theory.Theory, x: float, eta: float):
    """Calculate gegenbauer p_j to transform GPD from j to x space.

    Args:
        m: instance of the Theory class
        x, eta: GPD kinematical variables

    Returns:
        p_j[s, k, a]: s in range(npws), k in range(npts),
                     a in [Q, G, ...] i.e. m.evolution_basis

    Todo:
        Check conventions for various quark combinations.

    """
    # FIXME: vectorize p_j more universally so this looping is not needed
    pj = []
    for j in m.jpoints_pws:
        pjf = []
        rest_ebas = 0  # index of next item in evolution basis
        pjQp = gegenbauer.p_j(j, x, eta, 3/2)
        if 'G' in m.evolution_basis:
            pjQm = gegenbauer.p_j(j, -x, eta, 3/2)
            pjGp =  gegenbauer.p_j(j-1, x, eta, 5/2)
            pjGm =  gegenbauer.p_j(j-1, -x, eta, 5/2)
            pjG = -(pjGp+pjGm)  # MS06 (B.14)
            pjf.append(pjQp-pjQm)
            pjf.append(pjG)
            rest_ebas = 2  # we processed first two items
        rest_of_basis = m.evolution_basis[rest_ebas:]  # may be empty if only SI neded
        for k, item in enumerate(rest_of_basis):
            pjf.append(pjQp)
        pj.append(np.stack(pjf, axis=1))
    return np.stack(pj, axis=0)


def calc_wce(m: theory.Theory, Q2: float, process_class: str):
    """Calculate evolved Wilson coeffs for given Q2, for all PWs.

    Args:
        m: instance of the Theory
        Q2: process scale (so renormalization scale is
            Q2/m.rr2 and factorization scale is Q2/m.rf2)
        process_class: 'DIS' or 'DVCS'

    Returns:
        wce[s,k,b]: s in range(npwmax), k in range(npts), b in [Q,G,NSP]

    """
    # alpha_strong/2pi at renormalization, factorization and input scales:
    asmur2 = qcd.as2pf(m.p, m.nf, Q2/m.rr2, m.asp[m.p], m.r20)
    asmuf2 = qcd.as2pf(m.p, m.nf, Q2/m.rf2, m.asp[m.p], m.r20)
    asQ02 = qcd.as2pf(m.p, m.nf, m.Q02, m.asp[m.p], m.r20)
    R = asmuf2/asQ02
    wc = calc_wc(m, process_class)
    evola = evolution.evolop(m, R, process_class)
    p_mat = np.array([[1,],])   # LO version of p_mat, see below
    if m.p == 1:
        # p_mat: matrix that combines (LO, NLO) evolution operator and Wilson coeffs
        # while canceling NNLO term NLO*NLO:
        p_mat = np.array([[1, asmuf2], [asmur2, 0]])
    # 3. evolved Wilson coeff.
    wce = np.einsum('skpa,pq,skqab->skb', wc, p_mat, evola)
    return wce[..., :3]


def calc_wce_dvmp(m: theory.Theory, Q2: float):
    """Calculate evolved DVMP Wilson coeffs for given Q2, for all PWs.

    Args:
        m: instance of the Theory
        Q2: default renormalization and factorization scale
            (can be changed with m.rdvmpr2, m.rf2 and m.rdaf2)

    Returns:
        wce[s,k,g,j]: s in range(npwmax), k in range(npts),
                       g in gegens, j in [Q,G,NSP]

    Todo:
        Code duplication, merge with calc_wce()

    """
    process_class = 'DVMP'
    # alpha_strong/2pi at various scales:
    asmur2 = qcd.as2pf(m.p, m.nf, Q2/m.rdvmpr2, m.asp[m.p], m.r20)
    # GPD
    asmuf2 = qcd.as2pf(m.p, m.nf, Q2/m.rf2, m.asp[m.p], m.r20)
    asQ02 = qcd.as2pf(m.p, m.nf, m.Q02, m.asp[m.p], m.r20)
    R = asmuf2/asQ02
    # DA
    asdamuf2 = qcd.as2pf(m.p, m.nf, Q2/m.rdaf2, m.asp[m.p], m.r20)
    asdaQ02 = qcd.as2pf(m.p, m.nf, m.daQ02, m.asp[m.p], m.r20)
    RDA = asdamuf2/asdaQ02
    # unfortunately c1dvmp doesn't accept numpy array of Gegenbauers
    # so we have to build wc like this:
    wc = []  
    for g in m.gpoints:
        wc.append(calc_wc_dvmp(m, g))
    wc = np.stack(wc, axis=0)  # gegenbauer index is first
    # FIXME: next line is NS^+ to get agreement with tests but rho is C=-1
    evola_da = evolution.evolop_da(m, RDA, prty=+1) # +1->-1 !!
    evola = evolution.evolop(m, R, process_class)
    # p_mat: matrix that combines (LO, NLO) evolution operator and Wilson coeffs
    # while canceling NNLO term NLO*NLO:
    p_mat = np.array([1]).reshape((1, 1, 1))   # LO
    if m.p == 1:
        p_mat = np.zeros((2, 2, 2))
        p_mat[0, 0, 0] = 1.
        p_mat[1, 0, 0] = asmur2
        p_mat[0, 1, 0] = asmuf2
        p_mat[0, 0, 1] = asdamuf2
    # 3. evolved Wilson coeff. (evola for GPD, evola_ns for DA)
    wce = np.einsum('gskpi,pqr,skqij,gr->skgj', wc, p_mat, evola, evola_da)
    return wce

def calc_j2x(m: theory.Theory, x: float, eta: float, Q2: float):
    """Calculate j2x coeffs, combined with evolution operator.

    Args:
        m: instance of the Theory
        x: long. momentum fraction argument of GPD
        eta: skewness
        Q2: final evolution scale

    Returns:
        wce[s,k,i,j]: s in range(npwmax); k in range(npts); i,j in [Q,G,NSP]

    Todo:
        Implement general (eta != x) j to x transform.

    """
    asmuf2 = qcd.as2pf(m.p, m.nf, Q2/m.rf2, m.asp[m.p], m.r20)
    asQ02 = qcd.as2pf(m.p, m.nf, m.Q02, m.asp[m.p], m.r20)
    R = asmuf2/asQ02
    nebas = len(m.evolution_basis)
    wc = np.zeros((m.npws, m.npts, m.p + 1, nebas), dtype=complex)
    j = m.jpoints_pws
    if eta < 1e-8:
        # forward limit, PDF-like
        one = np.ones_like(wc[:, :, 0, 0]) * x**(-j-1)/np.pi
        rest_ebas = 0
        if 'G' in m.evolution_basis:
            wc[:, :, 0, 0] =  one     # Q
            wc[:, :, 0, 1] = x * one  # G
            rest_ebas = 2
        rest_of_basis = m.evolution_basis[rest_ebas:]
        for k, item in enumerate(rest_of_basis):
            wc[:, :, 0, rest_ebas + k] = one     # NS
        evola = evolution.evolop(m, R, 'DIS')
    elif abs(eta-x) < 1e-8 and x>0:
        # cross-over, border eta=x limit
        fac = _fshu(j) * x**(-j-1)/np.pi
        facnu = np.einsum('ij,i->ij', fac, eta**np.arange(0, 2*m.npws, 2))
        rest_ebas = 0
        if 'G' in m.evolution_basis:
            wc[:, :, 0, 0] = facnu  # Q
            wc[:, :, 0, 1] = facnu * 2*x/(3+j) # G
            rest_ebas = 2
        rest_of_basis = m.evolution_basis[rest_ebas:]
        for k, item in enumerate(rest_of_basis):
            wc[:, :, 0, rest_ebas + k] = facnu     # NS # WRONG!
        evola = evolution.evolop(m, R, 'DVCS')
    else:
        # FIXME: a bit of a kludge follows
        pj = np.squeeze(calc_pj(m, x, eta))
        eta2nu = eta**np.arange(0, 2*m.npws, 2)
        if pj.ndim == 2:
            wc[:, :, 0, 0] = - np.einsum('s,sj,sj->sj', eta2nu, pj,
                                       1 / np.sin(np.pi * j))
        else:
            wc[:, :, 0, :] = - np.einsum('s,sja,sj->sja', eta2nu, pj,
                                       1 / np.sin(np.pi * j))
        evola = evolution.evolop(m, R, 'DVCS')
    
    p_mat = np.array([[1,],])  # LO
    if m.p == 1:
        p_mat = np.array([[1, asmuf2], [0, 0]])
    # Return evolved Wilson coeff.
    # Note the difference w.r.t. usual wce, where we have inner product of
    # wc and evola, while here we have element-wise product, so that we
    # are in gpd.py able to get separately quark and gluon GPDs in x-space
    wce = np.einsum('skpa,pq,skqab->skab', wc, p_mat, evola)
    return wce
