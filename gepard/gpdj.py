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
   Careful: Old Fortran Gepard used slots 3 and 4 for
   NS+ and NS- evolution types.

Todo:
    * Add implementation of parameters
    * Add all used models
    * class/inheritance approach?

"""

from cmath import exp
from typing import Tuple

from scipy.special import loggamma

import gepard.qj
import gepard.special


def toy(j: complex, *args) -> Tuple[complex, complex, complex, complex]:
    """Return 'toy' (no params) singlet GPD ansatz."""
    singlet = 454760.7514415856 * exp(loggamma(0.5 + j)) / exp(loggamma(10.6 + j))
    gluon = 17.837861981813603 * exp(loggamma(-0.1 + j)) / exp(loggamma(4.7 + j))
    return (singlet, gluon, 0+0j, 0+0j)


def test(j: complex, t: float, par: dict) -> Tuple[complex, complex, complex, complex]:
    """Return simple testing singlet GPD ansatz."""
    singlet = par['ns'] / (1 - t/par['ms']**2)**3 / gepard.special.pochhammer(
          1.0 - par['al0s'] - par['alps']*t + j, 8) * gepard.special.pochhammer(
          2.0 - par['al0s'], 8)
    gluon = par['ng'] / (1 - t/par['mg']**2)**3 / gepard.special.pochhammer(
          1.0 - par['al0g'] - par['alpg']*t + j, 6) * gepard.special.pochhammer(
          2.0 - par['al0g'], 6)
    return (singlet, gluon, 0+0j, 0+0j)


def fit(j: complex, t: float, par: dict) -> Tuple[complex, complex, complex, complex]:
    """Return default fitting singlet GPD ansatz."""
    par['ng'] = 0.6 - par['ns']  # first sum-rule constraint
    singlet = (gepard.qj.qj(j, t, 9, par['ns'], par['al0s'], par['alps']) *
               gepard.qj.betadip(j, t, par['ms'], 0., 2))
    gluon = (gepard.qj.qj(j, t, 7, par['ng'], par['al0g'], par['alpg']) *
             gepard.qj.betadip(j, t, par['mg'], 0., 2))
    return (singlet, gluon, 0+0j, 0+0j)
