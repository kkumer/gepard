"""Testing code for GPD models."""

import gepard as g
import numpy as np
from pytest import approx, mark

j_test = 0.5+0.3j
t_test = -0.1
tb_test = -0.25
xi_test = 0.001
Q2_test = 4.

pt0 = g.DataPoint({'xi': xi_test, 't': t_test, 'Q2': Q2_test})
pt0.xB = 2.*pt0.xi/(1.+pt0.xi)
ptzero = g.DataPoint({'x': xi_test, 'eta': 0, 't': t_test, 'Q2': Q2_test})
pttraj = g.DataPoint({'x': xi_test, 'eta': xi_test, 't': t_test, 'Q2': Q2_test})

par_test = {'ns': 2./3. - 0.4, 'al0s': 1.1, 'alps': 0.25, 'ms2': 1.1**2,
            'ng': 0.4, 'al0g': 1.2, 'alpg': 0.25, 'mg2': 1.2**2}

par_fit = {'ns':  0.152039, 'al0s': 1.15751, 'alps': 0.15, 'ms2': 0.478391,
           'secs': -0.15152, 'this': 0.,  # 'ng': 0.4,  # provided by ns
           'al0g': 1.24732, 'alpg': 0.15, 'mg2': 0.7, 'secg': -0.81217, 'thig': 0.}

MP = 0.938272  # proton mass
par_bp = {'ns': 0, 'al0s': 1.1, 'alps': 0.15,
          'ms2': (2*MP)**2, 'delms2': MP**2, 'pows': 3,
          'ng': 0.5, 'al0g': 1.0, 'alpg': 0.15,
          'mg2': (2*MP)**2, 'delmg2': MP**2, 'powg': 2,
          'nu': 2.0, 'al0u': 0.5, 'alpu': 1.0,
          'mu2': (2*MP)**2, 'delmu2': MP**2, 'powu': 1,
          'nd': 1.0, 'al0d': 0.5, 'alpd': 1.0,
          'md2': (2*MP)**2, 'delmd2': MP**2, 'powd': 1}

#  'hard' ansatz:
par_bp_hard = {'ng': 0.4, 'al0g': 1.1 + 0.05, 'ns': 2./3. - 0.4}

# one of Dieter's fits
par_DM12 = {"ns": 0.152, "al0s": 1.1575, "alps": 0.1028, "ms2": 0.1814, "delms2": 0.,
            "pows": 2., "secs": 0.5142, "this": -0.2106, # "SKEWS": 0.,
            "al0g": 1.2473, "alpg": -0.023, "mg2": 0.2272, "delmg2": 0.,
            "powg": 2., "secg": -4.587, "thig": 1.7607, # "SKEWG": 0.,
            "kaps": 1.5265,  # kapg fixed by B_G sum rule
            "Ens": 0.152,  # for fitting should be pegged to ns
            "Eal0s": 1.1754, "Ealps": 0.0239, "Ems2": 0.1808, "Edelms2": 0.,
            "Epows": 2., "Esecs": 0.5235, "Ethis": -0.2194, # "ESKEWS": 0.,
            "Eal0g": 1.8486, "Ealpg": 0.05, "Emg2": 0.1487, "Edelmg2": 0.,
            "Epowg": 2., "Esecg": -4.6366, "Ethig": 1.8638} #, "ESKEWG": 0.}


#  -- Testing actual ansatze i.e. shapes of functions at single j-point


def test_gpd_toy():
    """Test simplest singlet GPD (no params)."""
    assert g.gpd.toy(j_test) == approx(((0.05957785753421334-0.07079531134708178j),
                                        (0.4061603004928009-0.8446849472313395j),
                                        0j, 0j))


def test_gpd_test():
    """Test testing singlet GPD."""
    assert g.gpd.test(j_test, t_test, par_test) == approx((
                             (0.2816440151875255-0.827096278856685j),
                             (0.32840732368726266-1.240107972227522j),
                             0j, 0j))


def test_gpd_fit():
    """Test fitting nl-SO(3) singlet GPD."""
    assert g.gpd.singlet_ng_constrained(j_test, t_test, par_fit) == approx((
                            (0.0952245874481982-0.5039416636123855j),
                            (0.14952730665086755-1.5915106823348937j),
                            0j, 0j))


def test_gpd_ansatz07():
    """Test singlet GPD from hep-ph/0703179 - hard gluons."""
    par_bp.update(par_bp_hard)
    sea, G, uv, dv = g.gpd.ansatz07(j_test, tb_test, par_bp)
    assert (sea+uv+dv, G) == approx((
                             (0.8252585606460219-1.20745720388647j),
                             (0.48405328160880784-1.2858277230367638j)))

# -- Testing complete GPDs on MB contour


def test_ConformalMoment_gpdH_test():
    """Test ConformalMoment GPD class with testing GPD H."""
    test_gpd = g.gpd.TestGPD()
    test_gpd.parameters.update(par_test)
    assert test_gpd.H(0.1, -0.2)[:1, :2] == approx(
            np.array([[1.65552601-0.00182673j, 3.1300559-0.00436946j]]))


def test_ConformalMoment_gpdH_fit():
    """Test ConformalMoment GPD class with testing GPD H."""
    fit_gpd = g.gpd.PWNormGPD()
    fit_gpd.parameters.update(par_fit)
    assert fit_gpd.H(0.1, -0.2)[:1, :2] == approx(
            np.array([[1.1665696086-0.00161121675988j, 5.59105109-0.0109293227j]]))


# @mark.slow
# def test_ConformalMoment_gpd_parallel():
#     """Test parallel evaluation of GPDs on the Mellin-Barnes contour."""
#     fit_gpd = g.gpd.PWNormGPD()
#     fit_gpd.parameters.update(par_fit)
#     assert fit_gpd.H_para(0.1, -0.2)[:1, :2] == approx(
#             np.array([[1.1665696086-0.00161121675988j, 5.59105109-0.0109293227j]]))

# class MyTheory(g.gpd.PWNormGPD, g.cff.MellinBarnesCFF):
    # pass


# def test_aux_H():
    # """Calculate some CFFs."""
    # cff = MyTheory(residualt='exp')
    # cff.parameters.update(par_DM12)
    # assert cff.ImH(pt0) == approx(1172.8817791633037)

# def test_aux_E():
    # """Calculate some CFFs."""
    # cff = MyTheory(residualt='exp')
    # cff.parameters.update(par_DM12)
    # assert cff.ImE(pt0) == approx(1728.4783865956076)

def test_GPDzero():
    """Calculate GPDs on forward trajectory eta=0."""
    gpd = g.gpd.PWNormGPD(residualt='exp')
    gpd.parameters.update(par_DM12)
    assert gpd.Hx(ptzero) == approx(
            np.array([1836.47, 6.9642, 0]), rel=1e-5)
    assert gpd.Ex(ptzero) == approx(
            np.array([3098.61, 0.0756646, 0]), rel=1e-5)


def test_GPDtraj():
    """Calculate GPDs on border/cross-over trajectory eta=x."""
    gpd = g.gpd.PWNormGPD(residualt='exp')
    gpd.parameters.update(par_DM12)
    assert gpd.Hx(pttraj) == approx(
            np.array([1344.02, 2.6866, 0]), rel=1e-5)
    assert gpd.Ex(pttraj) == approx(
            np.array([1980.69, 0.0358365, 0]), rel=1e-5)
