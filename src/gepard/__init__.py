"""Gepard --- Python package for analysis of generalized parton distributions."""

from importlib.metadata import version
__version__ = version(__name__)

# Some things are brought to gepard top namespace
# for user's convenience, and for easier preservation of
# backward compatibility

from .cff import (CFF, MellinBarnesCFF, DispersionCFF, PionPole,  # noqa: F401
                  DispersionFixedPoleCFF, DispersionFreePoleCFF,
                  HybridFixedPoleCFF, HybridFreePoleCFF, GoloskokovKrollCFF)
from .data import (DataPoint, DataSet, dset, loaddata,  # noqa: F401
                   select, list_data, describe_data)
from .dis import DIS  # noqa: F401
from .dvcs import DVCS, BM10, BMK, BM10ex, BM10tw2, hotfixedBMK  # noqa: F401
from .dvmp import DVMP, MellinBarnesTFF  # noqa: F401
from .eff import ZeroEFF, DipoleEFF, KellyEFF  # noqa: F401
from .fitter import MinuitFitter  # noqa: F401
from .gpd import GPD, ConformalSpaceGPD, TestGPD, PWNormGPD  # noqa: F401
from .kinematics import tmin, tmax, weight_BH, prepare  # noqa: F401
from .qcd import beta, as2pf  # noqa: F401
from .theory import Theory  # noqa: F401
