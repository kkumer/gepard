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

from . import constants, data, mellin, model, quadrature, special, wilson

#  ---- Building block - j-dependence ----


def qj(j: np.ndarray, t: float, poch: int, norm: float, al0: float,
        alp: float, alpf: float = 0, val: int = 0) -> np.ndarray:
    r"""GPD building block Q_j with reggeized t-dependence.

    Args:
        j: complex conformal moment
        t: momentum transfer
        poch: pochhammer (beta+1)
        norm: overall normalization
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
    r"""GPD residual dipole t-dependence function.

    Args:
        j: conformal moment
        t: momentum transfer
        m02: mass param
        delm2: mass param - interplay with j
        pp: exponent (=2 for dipole)

    Returns:
        Dipole residual t-dependence (beta from Eq. (19) of 0904.0458)

    """
    return 1. / (1. - t / (m02 + delm2*j))**pp


def betaexp(t: float, m02: float) -> float:
    r"""GPD residual exponential t-dependence function.

    Args:
        t: momentum transfer
        m02: mass parameter

    Returns:
        Exponential residual t-dependence (beta from Eq. (19) of 0904.0458)

    """
    return exp(t / (2 * m02))


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


def singlet_ng_constrained(j: np.ndarray, t: float, par: dict,
                           residualt: str = 'dipole') -> np.ndarray:
    r"""Singlet-only GPD ansatz, with ng parameter constrained by sum-rule.

    Args:
        j: conformal moment
        t: momentum transfer
        par: parameters dict
        residualt: residual t-dependence type ('dipole' or 'exp')

    Returns:
        GPD ansatz, as numpy array

    Notes:
        This ansatz is used for all published KM fits, for the sea parton part.

    """
    par['ng'] = 0.6 - par['ns']  # first sum-rule constraint
    if residualt == 'dipole':
        tdep_s = betadip(j, t, par['ms2'], 0., 2)
        tdep_g = betadip(j, t, par['mg2'], 0., 2)
    elif residualt == 'exp':
        tdep_s = betaexp(t, par['ms2'])
        tdep_g = betaexp(t, par['mg2'])
    else:
        raise ValueError("{} unknown. Use 'dipole' or 'exp'".format(residualt))

    singlet = (qj(j, t, 9, par['ns'], par['al0s'], par['alps']) * tdep_s)
    gluon = (qj(j, t, 7, par['ng'], par['al0g'], par['alpg']) * tdep_g)
    return np.array((singlet, gluon, np.zeros_like(gluon), np.zeros_like(gluon)))


def singlet_ng_constrained_E(j: np.ndarray, t: float, par: dict,
                             residualt: str = 'dipole') -> np.ndarray:
    r"""Singlet-only GPD ansatz, with ng parameter constrained by sum-rule.

    Args:
        j: conformal moment
        t: momentum transfer
        par: parameters dict
        residualt: residual t-dependence type ('dipole' or 'exp')

    Returns:
        GPD ansatz, as numpy array

    Notes:
        This ansatz is used for all published KM fits, for the sea parton part.

    """
    par['Eng'] = 0.6 - par['Ens']  # first sum-rule constraint
    if residualt == 'dipole':
        tdep_s = betadip(j, t, par['Ems2'], 0., 2)
        tdep_g = betadip(j, t, par['Emg2'], 0., 2)
    elif residualt == 'exp':
        tdep_s = betaexp(t, par['Ems2'])
        tdep_g = betaexp(t, par['Emg2'])
    else:
        raise ValueError("{} unknown. Use 'dipole' or 'exp'".format(residualt))

    singlet = (qj(j, t, 9, par['Ens'], par['Eal0s'], par['Ealps']) * tdep_s)
    gluon = (qj(j, t, 7, par['Eng'], par['Eal0g'], par['Ealpg']) * tdep_g)
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
        j: conformal moment
        t: momentum transfer
        type: 'soft', 'hard', 'softNS', 'hardNS'

    Returns:
        GPD ansatz, as numpy array

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


class GPD(model.ParameterModel):
    """Base class of all GPD models.

    Args:
        p: pQCD order (0 = LO, 1 = NLO, 2 = NNLO)
        scheme: pQCD scheme  ('msbar' or 'csbar')
        nf: number of active quark flavors
        Q02: Initial Q0^2 for GPD evolution.
        r20: Initial mu0^2 for alpha_strong definition.
        asp: alpha_strong/(2*pi) at scale r20 for (LO, NLO, NNLO)
        residualt: residual t dependence ('dipole' or 'exp')

    """
    def __init__(self, **kwargs) -> None:
        self.p = kwargs.setdefault('p', 0)
        self.scheme = kwargs.setdefault('scheme', 'msbar')
        self.nf = kwargs.setdefault('nf', 4)
        self.Q02 = kwargs.setdefault('Q02', 4.0)
        self.r20 = kwargs.setdefault('r20', 2.5)
        self.asp = kwargs.setdefault('asp', np.array([0.0606, 0.0518, 0.0488]))
        self.residualt = kwargs.setdefault('residualt', 'dipole')
        # scales
        self.rr2 = 1     # ratio of Q2/renorm. scale squared
        self.rf2 = 1     # ratio of Q2/GPD fact. scale sq.
        self.rdaf2 = 1   # ratio of Q2/DA fact. scale sq. (for DVMP)
        #
        # Model parameters
        all_pars = [
            "ns", "al0s", "alps", "ms2", "delms2", "pows", "secs", "this",
            "ng", "al0g", "alpg", "mg2", "delmg2", "powg", "secg", "thig",
            "kaps", "kapg",
            "Ens", "Eal0s", "Ealps", "Ems2", "Edelms2", "Epows", "Esecs", "Ethis",
            "Eng", "Eal0g", "Ealpg", "Emg2", "Edelmg2", "Epowg", "Esecg", "Ethig"]
        # subclasses should actually set values
        for par in all_pars:
            self.parameters[par] = self.parameters.setdefault(par, 0)
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
        # For DIS PDFs:
        self.frot_pdf = np.array([[1, 0, 0, 0],
                                  [0, 1, 0, 0],
                                  [0, 0, 0, 0]])
        # For DVMP
        self.frot_rho0_4 = np.array([[1, 0, 1, 1],
                                    [0, 1, 0, 0],
                                    [0., 0, 0., 0.]]) / np.sqrt(2)
        #                           # [3./20., 0, 5./12., 1./12.]]) / np.sqrt(2)
        # For j2x
        self.frot_j2x = self.frot_pdf
        super().__init__(**kwargs)


class ConformalSpaceGPD(GPD, mellin.MellinBarnes):
    """Base class of GPD models built in conformal moment space.

    Args:
        c: intersection of Mellin-Barnes curve with real axis
        phi: angle of Mellin-Barnes curve with real axis

    Notes:
        This just takes care of initialization of Mellin-Barnes
        contour points, and evolved Wilson coeffs, and MB Gauss
        integration weights. Actual choice and code for GPDs
        is provided by subclasses.

    """
    def __init__(self, **kwargs) -> None:
        self.c = kwargs.setdefault('c', 0.35)
        self.phi = kwargs.setdefault('phi', 1.57079632)
        npoints, weights = quadrature.mellin_barnes(self.c, self.phi)
        self.npts = len(npoints)
        self.npoints = npoints
        self.jpoints = npoints - 1
        self.wg = weights  # Gauss integration weights
        # Initial parameters:
        self.add_parameters({'ns': 2./3. - 0.4,
                             'al0s': 1.1,  'Eal0s': 1.1,
                             'alps': 0.25, 'Ealps': 0.25,
                             'ms2': 1.1, 'Ems2': 1.1,
                             'secs': 0., 'Esecs': 0,
                             'this': 0., 'Ethis': 0,
                             'kaps': 0.,
                             'ng': 0.4, 'Eng': 0.4,
                             'al0g': 1.2, 'Eal0g': 1.2,
                             'alpg': 0.25, 'Ealpg': 0.25,
                             'mg2': 1.2, 'Emg2': 1.2,
                             'secg': 0., 'Esecg': 0,
                             'thig': 0., 'Ethig': 0,
                             'kapg': 0.})
        mellin.MellinBarnes.__init__(self, **kwargs)
        super().__init__(**kwargs)

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

    def pw_strengths_E(self):
        """Strengths of SO(3) partial waves for gpd E."""
        return np.array([[1., 1., 1],
                         [self.parameters['Esecs'],
                             self.parameters['Esecg'], 0],
                         [self.parameters['Ethis'],
                             self.parameters['Ethig'], 0]])

    def H(self, eta: float, t: float) -> np.ndarray:
        """Return (npts, 4) array H_j^a for all j-points and 4 flavors."""
        return np.zeros((self.npts, 4), dtype=complex)

    def E(self, eta: float, t: float) -> np.ndarray:
        """Return (npts, 4) array E_j^a for all j-points and 4 flavors."""
        return np.zeros((self.npts, 4), dtype=complex)

    def Hx(self, pt: data.DataPoint) -> np.ndarray:
        """Return x-space GPD H.

        Args:
            pt: datapoint with kinematics info

        Returns:
            x-space GPD. 3-dim vector (singlet quark, gluon, non-singlet quark)
            is returned by transforming original conformal moment space (j-space) model.

        Todo:
            Non-singlet component is set to zero. We need to check
            normalization/symmetrization first.

        """
        # get "Wilson" coef., first PW is the only relevant one
        wce_j2x = wilson.calc_j2x(self, pt.x, pt.eta, pt.Q2)
        gpd_prerot = self.H(pt.eta, pt.t)
        gpd = np.einsum('fa,ja->jf', self.frot_j2x, gpd_prerot)
        mb_int_flav = self._j2x_mellin_barnes_integral(pt.x, pt.eta, wce_j2x, gpd)
        return mb_int_flav / np.pi

    def Ex(self, pt: data.DataPoint) -> np.ndarray:
        """Return x-space GPD E.

        Args:
            pt: datapoint with kinematics info

        Returns:
            x-space GPD. 3-dim vector (singlet quark, gluon, non-singlet quark)
            is returned by transforming original conformal moment space (j-space) model.

        Todo:
            Non-singlet component is set to zero. We need to check
            normalization/symmetrization first.

        """
        # get "Wilson" coef., first PW is the only relevant one
        wce_j2x = wilson.calc_j2x(self, pt.x, pt.eta, pt.Q2)
        gpd_prerot = self.E(pt.eta, pt.t)
        gpd = np.einsum('fa,ja->jf', self.frot_j2x, gpd_prerot)
        mb_int_flav = self._j2x_mellin_barnes_integral_E(pt.x, pt.eta, wce_j2x, gpd)
        return mb_int_flav / np.pi


class TestGPD(ConformalSpaceGPD):
    """Simple testing ansatz for GPDs."""

    def __init__(self, **kwargs) -> None:
        kwargs.setdefault('scheme', 'csbar')
        kwargs.setdefault('nf', 3)
        kwargs.setdefault('Q02', 1.0)
        kwargs.setdefault('r20', 2.5)
        kwargs.setdefault('asp', np.array([0.05, 0.05, 0.05]))
        super().__init__(**kwargs)

    def H(self, eta: float, t: float) -> np.ndarray:
        """Return (npts, 4) array H_j^a for all j-points and 4 flavors."""
        # For testing purposes, we use here sub-optimal non-numpy algorithm
        h = []
        for j in self.jpoints:
            h.append(test(j, t, self.parameters))
        return np.array(h)


class PWNormGPD(ConformalSpaceGPD):
    """Singlet-only model for GPDs with three SO(3) partial waves.

    Args:
        p: pQCD order (0 = LO, 1 = NLO, 2 = NNLO)
        scheme: pQCD scheme  ('msbar' or 'csbar')
        nf: number of active quark flavors
        Q02: Initial Q0^2 for GPD evolution.
        r20: Initial mu0^2 for alpha_strong definition.
        asp: alpha_strong/(2*pi) at scale r20 for (LO, NLO, NNLO)
        residualt: residual t dependence ('dipole' or 'exp')
        c: intersection of Mellin-Barnes curve with real axis
        phi: angle of Mellin-Barnes curve with real axis

    Notes:
        Subleading PWs are proportional to the leading one, and
        only their norms are fitting parameters. Norms of second
        PWs is given by parameers 'secs' (quarks) and 'secg' gluons,
        and norm of third PWs is given by 'this' and 'thig'.

        This is used for modelling sea partons in KM10-KM20 models.

    """
    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)

    def H(self, eta: float, t: float) -> np.ndarray:
        """Return (npts, 4) array H_j^a for all j-points and 4 flavors."""
        return singlet_ng_constrained(self.jpoints,
                                      t, self.parameters, self.residualt).transpose()

#     def H_para(self, eta: float, t: float) -> np.ndarray:
#         """Return (npts, 4) array H_j^a for all j-points and 4 flavors."""
#         # This multiprocessing version is actually 2x slower!
#         h = Parallel(n_jobs=20)(delayed(singlet_ng_constrained)(j, t,
#                                                                        self.parameters)
#                                 for j in self.jpoints)
#         return np.array(h)

    def E(self, eta: float, t: float) -> np.ndarray:
        """Return (npts, 4) array E_j^a for all j-points and 4 flavors."""
        # Implement BS+BG=0 sum rule that fixes 'kapg'
        self.parameters['kapg'] = - self.parameters['kaps'] * self.parameters['ns'] / (
                0.6 - self.parameters['ns'])
        kappa = np.array([self.parameters['kaps'], self.parameters['kapg'], 0, 0])
        return kappa * singlet_ng_constrained_E(
                self.jpoints, t, self.parameters, self.residualt).transpose()
