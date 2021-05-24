"""Testing code for GPD models."""

import numpy as np
from pytest import approx, mark

import gepard as g

j_test = 0.5+0.3j
t_test = -0.1
xi_test = 0.001
Q2_test = 4.

pt0 = g.DummyPoint({'xi': xi_test, 't': t_test, 'Q2': Q2_test})
pt0.xB = 2.*pt0.xi/(1.+pt0.xi)

par_test = {'ns': 2./3. - 0.4, 'al0s': 1.1, 'alps': 0.25, 'ms': 1.1,
            'ng': 0.4, 'al0g': 1.2, 'alpg': 0.25, 'mg': 1.2}

par_fit = {'ns':  0.152039, 'al0s': 1.15751, 'alps': 0.15, 'ms': 0.478391,
           'secs': -0.15152, 'this': 0.,  # 'ng': 0.4,  # provided by ns
           'al0g': 1.24732, 'alpg': 0.15, 'mg': 0.7, 'secg': -0.81217, 'thig': 0.}

#  -- Testing actual ansatze i.e. shapes of functions at single j-point


def test_gpdj_toy():
    """Test simplest singlet GPD (no params)."""
    assert g.gpdj.toy(j_test) == approx(((0.05957785753421334-0.07079531134708178j),
                                        (0.4061603004928009-0.8446849472313395j),
                                        0j, 0j))


def test_gpdj_test():
    """Test testing singlet GPD."""
    assert g.gpdj.test(j_test, t_test, par_test) == approx((
                             (0.2816440151875255-0.827096278856685j),
                             (0.32840732368726266-1.240107972227522j),
                             0j, 0j))


def test_gpdj_fit():
    """Test fitting nl-SO(3) singlet GPD."""
    assert g.gpdj.fit(j_test, t_test, par_fit) == approx((
                            (0.0952245874481982-0.5039416636123855j),
                            (0.14952730665086755-1.5915106823348937j),
                            0j, 0j))


# -- Testing complete GPDs on MB contour


def test_ConformalMoment_gpdH_test():
    """Test ConformalMoment GPD class with testing GPD H."""
    test_gpd = g.model.Test()
    test_gpd.parameters.update(par_test)
    assert test_gpd.gpd_H(0.1, -0.2)[:1, :2] == approx(
            np.array([[1.65552601-0.00182673j, 3.1300559-0.00436946j]]))


def test_ConformalMoment_gpdH_fit():
    """Test ConformalMoment GPD class with testing GPD H."""
    fit_gpd = g.model.Fit()
    fit_gpd.parameters.update(par_fit)
    assert fit_gpd.gpd_H(0.1, -0.2)[:1, :2] == approx(
            np.array([[1.1665696086-0.00161121675988j, 5.59105109-0.0109293227j]]))


@mark.slow
def test_ConformalMoment_gpd_parallel():
    """Test parallel evaluation of GPDs on the Mellin-Barnes contour."""
    fit_gpd = g.model.Fit()
    fit_gpd.parameters.update(par_fit)
    assert fit_gpd.gpd_H_para(0.1, -0.2)[:1, :2] == approx(
            np.array([[1.1665696086-0.00161121675988j, 5.59105109-0.0109293227j]]))


# def test_GPDtraj():
#     """Calculate GPDs on border trajectory xi=x."""
#     assert True
#     # assert m.gpdHtrajQ(pt0)/1000. == approx(1.34402)
#     # assert_almost_equal(m.gpdHtrajG(pt0), 2.6866, 4)
#     # assert_almost_equal(m.gpdEtrajQ(pt0)/1000., 1.98069, 5)
#     # assert_almost_equal(m.gpdEtrajG(pt0)*100, 3.58365, 5)


# def test_GPDzero():
#     """Calculate GPDs on forward trajectory xi=0."""
#     assert True
#     # assert m.gpdHzeroQ(pt0)/1000. == approx(1.83647)
#     # assert_almost_equal(m.gpdHzeroG(pt0), 6.9642, 4)
#     # assert_almost_equal(m.gpdEzeroQ(pt0)/1000., 3.09861, 5)
#     # assert_almost_equal(m.gpdEzeroG(pt0)*100, 7.56646, 5)