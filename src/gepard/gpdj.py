"""GPD models in conformal moment j-space.

Args:
      j: complex-space conformal moment
      t: momentum transfer squared
    par: dictionary {'parmeter1_name': value1, ...}

Returns:
   conformal moment of GPD as 4-dim flavor array:
   (1) -- singlet quark
   (2) -- gluon
   (3) -- u_valence
   (4) -- d_valence

Notes:
   Valence means "valence-like GPD". See hep-ph/0703179
   for description.

Todo:
    * Add implementation of parameters
    * Add all used models
    * class/inheritance approach?

"""

from cmath import exp
from typing import Tuple

import numpy as np
from scipy.special import loggamma  # type: ignore

from . import constants, qj, special


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
    singlet = (qj.qj(j, t, 9, par['ns'], par['al0s'], par['alps']) *
               qj.betadip(j, t, par['ms2'], 0., 2))
    gluon = (qj.qj(j, t, 7, par['ng'], par['al0g'], par['alpg']) *
             qj.betadip(j, t, par['mg2'], 0., 2))
    return np.array((singlet, gluon, np.zeros_like(gluon), np.zeros_like(gluon)))


def ansatz07(j: np.ndarray, t: float, par: dict) -> np.ndarray:
    """GPD ansatz from paper hep-ph/0703179."""
    uv = (qj.qj(j, t, 4, par['nu'], par['al0u'], par['alpu'], val=1) *
          qj.betadip(j, t, par['mu2'], par['delmu2'], par['powu']))
    dv = (qj.qj(j, t, 4, par['nd'], par['al0d'], par['alpd'], val=1) *
          qj.betadip(j, t, par['md2'], par['delmd2'], par['powd']))
    sea = (qj.qj(j, t, 8, par['ns'], par['al0s'], par['alps']) *
           qj.betadip(j, t, par['ms2'], par['delms2'], par['pows']))
    gluon = (qj.qj(j, t, 6, par['ng'], par['al0g'], par['alpg']) *
             qj.betadip(j, t, par['mg2'], par['delmg2'], par['powg']))
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
    uv = qj.qj(j, t, pochv, par['nu'], par['al0v'],
                      alpf=0, alp=par['alpv'], val=1)
    uv = uv / mjt
    dv = qj.qj(j, t, pochv, par['nd'], par['al0v'],
                      alpf=0, alp=par['alpv'], val=1)
    dv = dv / mjt
    sea = qj.qj(j, t, pochs, par['nsea'], par['al0s'],
                       alpf=0, alp=par['alps'])
    sea = sea / mjt**3
    gluon = qj.qj(j, t, pochg, par['ng'], par['al0g'],
                         alpf=0, alp=par['alpg'])
    gluon = gluon / mjt**2
    return np.array((sea, gluon, uv, dv))
