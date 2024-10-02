"""Testing dependence on renormalization and factorization scales.

These tests correspond to Tables 1 and 2, column xi=0.1 from TowardsDVCS
paper and are also reproduced in notebook npb07.ipynb available at
https://github.com/openhep/dvcs-old

(Numbers for Table 2 are slightly different from the published paper,
which had bug in the non-diagonal term of singlet evolution operator.
See notebook npb07.ipynb for correct table.)

"""

import gepard as g
import numpy as np
from pytest import approx, mark


class MyTheory(g.ConformalSpaceGPD, g.MellinBarnesCFF):
    """GPD/CFF model from the Towards DVCS paper."""
    
    def __init__(self, type, **kwargs) -> None:
        # defaults from old Fortran's GEPARD.INI
        kwargs.setdefault('p', 1)
        kwargs.setdefault('scheme', 'csbar')
        kwargs.setdefault('nf', 4)
        kwargs.setdefault('Q02', 2.5)
        kwargs.setdefault('asp', np.array([0.05, 0.05, 0.05]))
        kwargs.setdefault('phi', 1.9)
        self.type = type
        super().__init__(**kwargs)
        R = 0.5  # ratio sbar/ubar
        self.frot = np.array([[1, 0, 1, 1],
                              [0, 1, 0, 0],
                              [-R/(2+R), 0, 1, -1]])

    def H(self, eta: float, t: float) -> np.ndarray:
        return g.gpd.ansatz07_fixed(self.jpoints, t, self.type).transpose()


@mark.slow
def test_mu_renorm_NS():
    "Testing NLO non-singlet DVCS dependence on renormalization scale."
    pt = g.data.DataPoint({'xi': 0.1, 'Q2': 4, 't': -0.25})
    th0 = MyTheory('softNS', p=1, scheme='msbar', rr2=0.5)
    th0.dvcs_charges = (0, 0, 1/6)
    th1 = MyTheory('softNS', p=1, scheme='msbar', rr2=1)
    th1.dvcs_charges = (0, 0, 1/6)
    th2 = MyTheory('softNS', p=1, scheme='msbar', rr2=2)
    th2.dvcs_charges = (0, 0, 1/6)
    H0 = np.sqrt(th0.ReH(pt)**2 + th0.ImH(pt)**2)
    H1 = np.sqrt(th1.ReH(pt)**2 + th1.ImH(pt)**2)
    H2 = np.sqrt(th2.ReH(pt)**2 + th2.ImH(pt)**2)
    deli = 100*(H0-H2)/H1
    assert  deli == approx(7.574170899296725)

@mark.slow
def test_mu_fact_NS():
    "Testing NLO non-singlet DVCS dependence on factorization scale."
    pt = g.data.DataPoint({'xi': 0.1, 'Q2': 4, 't': -0.25})
    th0 = MyTheory('softNS', p=1, scheme='msbar', rf2=0.5)
    th0.dvcs_charges = (0, 0, 1/6)
    th1 = MyTheory('softNS', p=1, scheme='msbar', rf2=1)
    th1.dvcs_charges = (0, 0, 1/6)
    th2 = MyTheory('softNS', p=1, scheme='msbar', rf2=2)
    th2.dvcs_charges = (0, 0, 1/6)
    H0 = np.sqrt(th0.ReH(pt)**2 + th0.ImH(pt)**2)
    H1 = np.sqrt(th1.ReH(pt)**2 + th1.ImH(pt)**2)
    H2 = np.sqrt(th2.ReH(pt)**2 + th2.ImH(pt)**2)
    deli = 100*(H0-H2)/H1
    assert  deli == approx(0.6696785048040836)


@mark.slow
def test_mu_renorm_SI():
    "Testing NLO singlet DVCS dependence on renormalization scale."
    pt = g.data.DataPoint({'xi': 0.1, 'Q2': 4, 't': -0.25})
    th0 = MyTheory('soft', p=1, scheme='msbar', rr2=0.5)
    th0.dvcs_charges = (5/18, 5/18, 0)
    th1 = MyTheory('soft', p=1, scheme='msbar', rr2=1)
    th1.dvcs_charges = (5/18, 5/18, 0)
    th2 = MyTheory('soft', p=1, scheme='msbar', rr2=2)
    th2.dvcs_charges = (5/18, 5/18, 0)
    H0 = np.sqrt(th0.ReH(pt)**2 + th0.ImH(pt)**2)
    H1 = np.sqrt(th1.ReH(pt)**2 + th1.ImH(pt)**2)
    H2 = np.sqrt(th2.ReH(pt)**2 + th2.ImH(pt)**2)
    deli = 100*(H0-H2)/H1
    assert  deli == approx(10.851472760455286)

@mark.slow
def test_mu_fact_SI():
    "Testing NLO singlet DVCS dependence on factorization scale."
    pt = g.data.DataPoint({'xi': 0.1, 'Q2': 4, 't': -0.25})
    th0 = MyTheory('soft', p=1, scheme='msbar', rf2=0.5)
    th0.dvcs_charges = (5/18, 5/18, 0)
    th1 = MyTheory('soft', p=1, scheme='msbar', rf2=1)
    th1.dvcs_charges = (5/18, 5/18, 0)
    th2 = MyTheory('soft', p=1, scheme='msbar', rf2=2)
    th2.dvcs_charges = (5/18, 5/18, 0)
    H0 = np.sqrt(th0.ReH(pt)**2 + th0.ImH(pt)**2)
    H1 = np.sqrt(th1.ReH(pt)**2 + th1.ImH(pt)**2)
    H2 = np.sqrt(th2.ReH(pt)**2 + th2.ImH(pt)**2)
    deli = 100*(H0-H2)/H1
    assert  deli == approx(-3.437367000034853)

