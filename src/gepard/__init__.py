"""Gepard --- Python package for analysis of generalized parton distributions."""

# Some things are brought to gepard top namespace
# for user's convenience, and easier preservation of
# backward compatibility

from .cff import (CFF, MellinBarnesCFF, DispersionCFF, PionPole,
        DispersionFixedPoleCFF, DispersionFreePoleCFF,
        HybridFixedPoleCFF, HybridFreePoleCFF)
from .data import DataPoint, DataSet, dset  # noqa: F401
from .dis import DIS
from .dvcs import DVCS, BM10, BMK, BM10ex, BM10tw2, hotfixedBMK
from .dvmp import DVMP, MellinBarnesTFF
from .eff import ZeroEFF, DipoleEFF, KellyEFF
from .fitter import FitterMinuit
from .gpd import GPD, ConformalSpaceGPD, TestGPD, PWNormGPD
from .kinematics import tmin, tmax, weight_BH, prepare
from .qcd import beta, as2pf
from .theory import Theory
from .utils import fill_kinematics, select, list_data, listchis, describe_data
