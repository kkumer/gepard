"""Gepard --- Python package for analysis of generalized parton distributions."""

# Import what's needed into gepard's namespace

import gepard.adim  # noqa: F401
import gepard.c1dvcs  # noqa: F401
import gepard.c1dvmp  # noqa: F401
import gepard.evolution  # noqa: F398
import gepard.evolc  # noqa: F398
import gepard.fitter  # noqa: F401
import gepard.gpdj  # noqa: F401
import gepard.kinematics  # noqa: F401
import gepard.plots  # noqa: F401
import gepard.qcd  # noqa: F401
import gepard.qj  # noqa: F401
import gepard.quadrature  # noqa: F401
import gepard.special  # noqa: F401
import gepard.theory  # noqa: F401
import gepard.utils  # noqa: F401

from .constants import Mp, Mp2  # noqa: F401
from .data import DataPoint, DataSet  # noqa: F401
from .model import ConformalSpaceGPD, MellinBarnesModel  # noqa: F401
