"""Tests for DVCS. Corresponds to Towards DVCS paper hep-ph/0703179"""

# see also dvcs-old on github openhep repo for more detailed comparison with paper

import gepard as g
import numpy as np
from pytest import approx, fixture, mark


class MyTheory(g.ConformalSpaceGPD, g.MellinBarnesCFF):
    """GPD/CFF model from the paper - singlet and nonsinglet part."""
    
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

th_ns = MyTheory('hardNS', p=1, scheme='msbar')
qns = 1/6
th_ns.dvcs_charges = (0, 0, qns)   # NS part only

th_si = MyTheory('hard', p=1, scheme='msbar')
qs = 5/18
th_si.dvcs_charges = (qs, qs, 0)   # NS part only


def DeltaK(th, xi, Q2):
    Habs = {}
    for qq in [2.5, Q2]:
        pt = g.data.DataPoint({'xi': xi, 'Q2': qq, 't': -0.25})
        Habs[qq] = np.sqrt(th.ReH(pt)**2 + th.ImH(pt)**2)
    return 100*(Habs[Q2]/Habs[2.5] - 1)


@mark.slow
def test_TWDdvcs_Fig8_NLO_MS_NS():
    """Calculate evolved MS-bar NLO NS CFF."""
    reslo = DeltaK(th_ns, 0.01, 10)
    reshi = DeltaK(th_ns, 0.5, 10)
    # following agrees with DM to best than percent
    assert reslo == approx(-0.5478313, rel=1e-7)
    assert reshi == approx(-11.37569, rel=1e-7)


@mark.slow
def test_TWDdvcs_Fig10_NLO_MS_SI():
    """Calculate evolved MS-bar NLO SI CFF."""
    reslo = DeltaK(th_si, 1e-5, 5)
    reshi = DeltaK(th_si, 0.5, 5)
    # following agrees with DM to best than percent
    assert reslo == approx(59.09931981765, rel=1e-7)
    assert reshi == approx(2.787284080654, rel=1e-7)
