"""GPD models in conformal moment j-space.

Args:
      j: complex-space conformal moment
      t: momentum transfer squared
    par: dictionary {'parmeter1_name': value1, ...}

Returns:
   conformal moments of GPD as npts x 4 array, where
   npts is number of points on Mellin-Barnes contour and
   4 flavors are:
   (1) -- singlet quark
   (2) -- gluon
   (3) -- u_valence
   (4) -- d_valence

   Valence here means "valence-like GPD". See hep-ph/0703179
   for description.

"""

from cmath import exp
from typing import Tuple

import numpy as np
from scipy.special import loggamma  # type: ignore

from . import constants, model, quadrature, special

#  ---- Building block - j-dependence ----


def qj(j: np.ndarray, t: float, poch: int, norm: float, al0: float,
        alp: float, alpf: float = 0, val: int = 0) -> np.ndarray:
    r"""GPD building block Q_j with reggeized t-dependence.

    Args:
        j: complex conformal moment
        t: momentum transfer
        poch: pochhammer (beta+1)
        al0: Regge intercept
        alp: Regge slope (leading pole)
        alpf: Regge slope (full)
        val: 0 for sea, 1 for valence partons

    Returns:
        conformal moment of GPD with Reggeized t-dependence,

    Examples:
        >>> qj(0.5+0.3j, 0.,  4, 3, 0.5, -0.8)  #doctest: +ELLIPSIS
        (5.668107685...-4.002820...j)

    Notes:
        There are two implementations of Regge t dependence, one uses
        `alpf` :math:`\alpha'` parameter (`alp` should be set to zero)

        .. math::

            Q_j = N \frac{B(1-\alpha(t)+j, \beta+1)}{B(2-\alpha_0, \beta+1)}

        This is then equal to Eq. (41) of arXiv:hep-ph/0605237
        without residual :math:`F(\Delta^2)`.

        The default uses `alp` parameter, and expression is
        from Eq. (40) or (41) of arXiv:0904.0458, which takes into account
        only the leading pole:

        .. math::

            Q_j = N \frac{B(1-\alpha_0+j, \beta+1)}{B(2-\alpha_0, \beta+1)}
                  \frac{1+j-\alpha_0}{1+j-\alpha(t)}

        and is again defined without residual beta(t) from Eq. (19).


    """
    assert alp*alpf == 0   # only one version of alpha' can be used
    alpt = al0 + alp * t
    qj = (
        norm
        * special.pochhammer(2 - val - al0 - alpf * t, poch)
        / special.pochhammer(1 - al0 + j, poch)
        * (1 + j - al0)
        / (1 + j - alpt)
    )
    return qj


#  ---- Building block - residual t-dependence ----

def betadip(j: np.ndarray, t: float, m02: float, delm2: float, pp: int) -> np.ndarray:
    r"""GPD residual dipole t-dependence function beta from Eq. (19) of NPB.

    Args:
        j: conformal moment
        t: momentum transfer
        m02: mass param
        delm2: mass param - interplay with j
        pp: exponent (=2 for dipole)

    Returns:
        Dipole residual t-dependence
    """
    return 1. / (1. - t / (m02 + delm2*j))**pp


#  ---- Ansaetze for GPD shapes  ----

def toy(j: complex, *args) -> Tuple[complex, complex, complex, complex]:
    """Return 'toy' (no params) singlet GPD ansatz."""
    singlet = 454760.7514415856 * exp(loggamma(0.5 + j)) / exp(loggamma(10.6 + j))
    gluon = 17.837861981813603 * exp(loggamma(-0.1 + j)) / exp(loggamma(4.7 + j))
    return (singlet, gluon, 0+0j, 0+0j)


def test(j: complex, t: float, par: dict) -> Tuple[complex, complex, complex, complex]:
    """Return simple testing singlet GPD ansatz."""
    singlet = par['ns'] / (1 - t/par['ms2'])**3 / special.pochhammer(
          1.0 - par['al0s'] - par['alps']*t + j, 8) * special.pochhammer(
          2.0 - par['al0s'], 8)
    gluon = par['ng'] / (1 - t/par['mg2'])**3 / special.pochhammer(
          1.0 - par['al0g'] - par['alpg']*t + j, 6) * special.pochhammer(
          2.0 - par['al0g'], 6)
    return (singlet, gluon, 0+0j, 0+0j)


def singlet_ng_constrained(j: np.ndarray, t: float, par: dict) -> np.ndarray:
    """Singlet-only GPD ansatz, with ng param constrained by sum-rule.

    Notes:
        This ansatz is used for all published KM fits, for sea parton part.
    """
    par['ng'] = 0.6 - par['ns']  # first sum-rule constraint
    singlet = (qj(j, t, 9, par['ns'], par['al0s'], par['alps']) *
               betadip(j, t, par['ms2'], 0., 2))
    gluon = (qj(j, t, 7, par['ng'], par['al0g'], par['alpg']) *
             betadip(j, t, par['mg2'], 0., 2))
    return np.array((singlet, gluon, np.zeros_like(gluon), np.zeros_like(gluon)))


def ansatz07(j: np.ndarray, t: float, par: dict) -> np.ndarray:
    """GPD ansatz from paper hep-ph/0703179."""
    uv = (qj(j, t, 4, par['nu'], par['al0u'], par['alpu'], val=1) *
          betadip(j, t, par['mu2'], par['delmu2'], par['powu']))
    dv = (qj(j, t, 4, par['nd'], par['al0d'], par['alpd'], val=1) *
          betadip(j, t, par['md2'], par['delmd2'], par['powd']))
    sea = (qj(j, t, 8, par['ns'], par['al0s'], par['alps']) *
           betadip(j, t, par['ms2'], par['delms2'], par['pows']))
    gluon = (qj(j, t, 6, par['ng'], par['al0g'], par['alpg']) *
             betadip(j, t, par['mg2'], par['delmg2'], par['powg']))
    return np.array((sea, gluon, uv, dv))


def ansatz07_fixed(j: np.ndarray, t: float, type: str) -> np.ndarray:
    """GPD ansatz from hep-ph/0703179 with fixed parameters.

    Args:
        type: 'soft', 'hard', 'softNS', 'hardNS'

    Notes:
        This is the same as ansatz07, only instead of passing
        parameter dict, user passes type string choosing
        particular fixed parameter choices from the paper above.
    """
    # a.k.a. 'FITBP' ansatz from Fortran Gepard
    par = {'al0s': 1.1, 'alps': 0.15, 'alpg': 0.15,
           'nu': 2, 'nd': 1, 'al0v': 0.5, 'alpv': 1}
    if type[:4] == 'hard':
        par['ng'] = 0.4
        par['al0g'] = par['al0s'] + 0.05
        if type[-2:] == 'NS':
            par['nsea'] = 4/15
        else:
            par['nsea'] = 2/3 - par['ng']
    elif type[:4] == 'soft':
        par['ng'] = 0.3
        par['al0g'] = par['al0s'] - 0.2
        if type[-2:] == 'NS':
            par['nsea'] = 0
        else:
            par['nsea'] = 2/3 - par['ng']
    pochs = 8
    pochg = 6
    pochv = 4
    mjt = 1 - t / (constants.Mp2*(4+j))
    uv = qj(j, t, pochv, par['nu'], par['al0v'],
                      alpf=0, alp=par['alpv'], val=1)
    uv = uv / mjt
    dv = qj(j, t, pochv, par['nd'], par['al0v'],
                      alpf=0, alp=par['alpv'], val=1)
    dv = dv / mjt
    sea = qj(j, t, pochs, par['nsea'], par['al0s'],
                       alpf=0, alp=par['alps'])
    sea = sea / mjt**3
    gluon = qj(j, t, pochg, par['ng'], par['al0g'],
                         alpf=0, alp=par['alpg'])
    gluon = gluon / mjt**2
    return np.array((sea, gluon, uv, dv))


# ---- Full GPD models  (for all j-points) ----


class ConformalSpaceGPD(model.ParameterModel):
    """Base class of GPD models built in conformal moment space."""

    def __init__(self, p: int = 0, scheme: str = 'msbar', nf: int = 4,
            q02: float = 4.0, r20: float = 2.5,
            asp: np.array = np.array([0.0606, 0.0518, 0.0488]),
            c: float = 0.35, phi: float = 1.57079632) -> None:
        """Init and pre-calculate stuff.

        Args:
            p: pQCD order (0 = LO, 1 = NLO, 2 = NNLO)
            scheme: pQCD scheme  (msbar or csbar)
            nf: number of active quark flavors
            q02: Initial Q0^2 for pQCD evolution.
            r20: Initial mu0^2 for alpha_strong definition.
            asp: alpha_strong/(2*pi) at scale r20 for (LO, NLO, NNLO)
            c: intersection of Mellin-Barnes curve with real axis
            phi: angle of Mellin-Barnes curve with real axis

        Notes:
            This just takes care of initialization of Mellin-Barnes
            contour points, and evolved Wilson coeffs, and MB Gauss
            integration weights. Actual choice and code for GPDs
            is provided by subclasses.

        """
        self.p = p
        self.scheme = scheme
        self.nf = nf
        self.q02 = q02
        self.r20 = r20
        self.asp = asp
        self.c = c
        self.phi = phi
        npoints, weights = quadrature.mellin_barnes(self.c, self.phi)
        self.npts = len(npoints)
        self.npoints = npoints
        self.jpoints = npoints - 1
        self.wg = weights  # Gauss integration weights
        # Initial parameters:
        self.parameters = {'ns': 2./3. - 0.4,
                           'al0s': 1.1,
                           'alps': 0.25,
                           'ms2': 1.1,
                           'secs': 0.,
                           'this': 0.,
                           'kaps': 0.,
                           'ng': 0.4,
                           'al0g': 1.2,
                           'alpg': 0.25,
                           'mg2': 1.2,
                           'secg': 0.,
                           'thig': 0.,
                           'kapg': 0.}
        # Flavor rotation matrix.
        # ----------------------
        # It transforms GPDs from flavor basis to evolution basis.
        # By default, evolution basis is 3-dim (SIG, G, NS+)
        # (NS- is not completely implemented yet)
        # while default flavor basis is (sea,G,uv,dv).
        # User is free to use more complicated flavor structure of model
        #
        # Default matrix that follows is appropriate for low-x DVCS,
        # with singlet-only contribution.
        # For definitions of sea-like and valence-like GPDs, sea, uv, dv
        # see hep-ph/0703179
        self.frot = np.array([[1, 0, 1, 1],
                              [0, 1, 0, 0],
                              [0, 0, 0, 0]])
        # squared DVCS charge factors
        # This might belong to CFF code
        if self.nf == 3:
            qs = 2/9
            qns = 1/9
        else:  # nf = 4
            qs = 5/18
            qns = 1/6
        self.dvcs_charges = (qs, qs, qns)
        super().__init__()

    def pw_strengths(self):
        """Strengths of SO(3) partial waves."""
        # We take maximally three partial waves atm:
        # pw_strengths = (no. pws x no. flavors)
        # flavors are (Q, G, NSP)
        return np.array([[1., 1., 1],
                         [self.parameters['secs'],
                             self.parameters['secg'], 0],
                         [self.parameters['this'],
                             self.parameters['thig'], 0]])

    def gpd_H(self, eta: float, t: float) -> np.ndarray:
        """Return (npts, 4) array H_j^a for all j-points and 4 flavors."""
        return np.zeros((self.npts, 4), dtype=complex)

    def gpd_E(self, eta: float, t: float) -> np.ndarray:
        """Return (npts, 4) array E_j^a for all j-points and 4 flavors."""
        return np.zeros((self.npts, 4), dtype=complex)


class TestGPD(ConformalSpaceGPD):
    """Simple testing ansatz for GPDs."""

    def __init__(self, **kwargs) -> None:
        """See parent `ConformalSpaceGPD` class for docs."""
        kwargs.setdefault('scheme', 'csbar')
        kwargs.setdefault('nf', 3)
        kwargs.setdefault('q02', 1.0)
        super().__init__(**kwargs)
        self.nf = 3
        self.q02 = 1.0
        self.asp = np.array([0.05, 0.05, 0.05])
        self.r20 = 2.5

    def gpd_H(self, eta: float, t: float) -> np.ndarray:
        """Return (npts, 4) array H_j^a for all j-points and 4 flavors."""
        # For testing purposes, we use here sub-optimal non-numpy algorithm
        h = []
        for j in self.jpoints:
            h.append(test(j, t, self.parameters))
        return np.array(h)


class PWNormGPD(ConformalSpaceGPD):
    """Singlet-only model for GPDs with three SO(3) partial waves.

    Notes:
        Subleading PWs are proportional to the leading one, and
        only their norms are fitting parameters. Norms of second
        PWs is given by parameers 'secs' (quarks) and 'secg' gluons,
        and norm of third PWs is given by 'this' and 'thig'.
    """

    def __init__(self, **kwargs) -> None:
        """See parent `ConformalSpaceGPD` class for docs."""
        kwargs.setdefault('scheme', 'msbar')
        kwargs.setdefault('nf', 4)
        kwargs.setdefault('q02', 4.0)
        super().__init__(**kwargs)

    def gpd_H(self, eta: float, t: float) -> np.ndarray:
        """Return (npts, 4) array H_j^a for all j-points and 4 flavors."""
        return singlet_ng_constrained(self.jpoints,
                                             t, self.parameters).transpose()

#     def gpd_H_para(self, eta: float, t: float) -> np.ndarray:
#         """Return (npts, 4) array H_j^a for all j-points and 4 flavors."""
#         # This multiprocessing version is actually 2x slower!
#         h = Parallel(n_jobs=20)(delayed(singlet_ng_constrained)(j, t,
#                                                                        self.parameters)
#                                 for j in self.jpoints)
#         return np.array(h)

    def gpd_E(self, eta: float, t: float) -> np.ndarray:
        """Return (npts, 4) array E_j^a for all j-points and 4 flavors."""
        kappa = np.array([self.parameters['kaps'], self.parameters['kapg'], 0, 0])
        return kappa * singlet_ng_constrained(self.jpoints, t,
                                                     self.parameters).transpose()

