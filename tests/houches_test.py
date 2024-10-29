"""Testing forward limit w.r.t. Les Houches benchmarks."""

import gepard as g
import numpy as np
from pytest import approx, mark
from scipy.special import gamma


class HouchesGPD(g.PWNormGPD):
    """Les Houches benchmark model from hep-ph/0204316"""

    def __init__(self, **kwargs):
        kwargs.setdefault('scheme', 'msbar')                                                              
        kwargs.setdefault('nf', 4)                                                                        
        kwargs.setdefault('Q02', 2.0)                                                                     
        kwargs.setdefault('r20', 2.0)                                                                     
        kwargs.setdefault('asp', np.array([0.35, 0.35, 0.35])/2/np.pi)
        super().__init__(**kwargs)

    def pw_strengths(self):
        res = np.zeros((self.npws, len(self.flavor_basis)))
        res[0, :] = np.ones_like(res[0, :])
        return res          

    def H(self, eta, t):
        """GPD H at input scale."""
        betas = np.array([3, 4, 5, 6, 7])
        norms = np.array([5.1072, 3.06432, 1.7, 0.1939875,  0.1939875])
        alphas = np.array([0.2, 0.2, 1.1, 1.1, 1.1])
        res = np.zeros((self.npts, len(self.flavor_basis)), dtype=complex)
        for k, (b, n, a) in enumerate(zip(betas, norms, alphas)):
            res[:, k] = n * gamma(b + 1) / g.special.pochhammer(self.jpoints + 1 - a, b + 1)
        res[:, 5] = 0.2 * (res[:, 3] + res[:, 4])  # sbar = (ubar + dbar)/5
        return res




def test_LH_LO():
    """Les Houches LO evolution."""
    m = HouchesGPD(p=0, accuracy=5, phi=1.9, extended=True, nf=4)
    m.frot_pdf = np.array([[1, 1, 0, 2, 2, 2, 2],
                     [0, 0, 1, 0, 0, 0, 0],
                     [1, 1, 0, 0, 0, 0, 0],
                     [1,-1, 0,-2, 2, 0, 0],
                     [1,-1, 0, 0, 0, 0, 0],
                     [1, 1, 0, 2, 2,-4, 0],
                     [1, 1, 0, 2, 2, 2,-6]
                          ])
    m.antifrot_pdf = np.array([[0, 0, 1/2, 0, 1/2, 0, 0],
                 [0, 0, 1/2, 0,-1/2, 0, 0],
                 [0, 1, 0, 0, 0, 0, 0],
                 [1/8, 0, -1/4, -1/4, 1/4, 1/12, 1/24],
                 [1/8, 0, -1/4,  1/4,-1/4, 1/12, 1/24],
                 [1/8, 0, 0, 0, 0,-1/6, 1/24],
                 [1/8, 0, 0, 0, 0, 0, -1/8]])
    m.frot_j2x = m.frot_pdf
    m.flavor_basis = ['uv', 'dv', 'g', 'dbar', 'ubar', 'sbar', 'cbar']
    m.evolution_basis = ['S', 'G', 'V', 'V3p', 'V3m', 'V8p', 'V15p']
    ##
    x = 0.1
    Q2 = 1.e4
    xHx = x * m.Hx(x, 0, 0, Q2)
    ## Numbers below are from LesHouches table and Pegasus software
    assert (xHx[0], xHx[1], xHx[3]-xHx[4], 2*(xHx[3]+xHx[4]),
            2*xHx[-2], 2*xHx[-1], xHx[2]/x) == approx((0.57267,
            0.28413, 0.010470, 0.40832, 0.11698, 0.058864, 0.88766), abs=1e-5)


@mark.slow
def test_LH_NLO():
    """Les Houches NLO evolution."""
    m = HouchesGPD(p=1, scheme='msbar', accuracy=7, phi=1.9, extended=True, nf=4)
    m.frot_pdf = np.array([[1, 1, 0, 2, 2, 2, 2],
                     [0, 0, 1, 0, 0, 0, 0],
                     [1, 1, 0, 0, 0, 0, 0],
                     [1,-1, 0,-2, 2, 0, 0],
                     [1,-1, 0, 0, 0, 0, 0],
                     [1, 1, 0, 2, 2,-4, 0],
                     [1, 1, 0, 2, 2, 2,-6]
                          ])
    m.antifrot_pdf = np.array([[0, 0, 1/2, 0, 1/2, 0, 0],
                 [0, 0, 1/2, 0,-1/2, 0, 0],
                 [0, 1, 0, 0, 0, 0, 0],
                 [1/8, 0, -1/4, -1/4, 1/4, 1/12, 1/24],
                 [1/8, 0, -1/4,  1/4,-1/4, 1/12, 1/24],
                 [1/8, 0, 0, 0, 0,-1/6, 1/24],
                 [1/8, 0, 0, 0, 0, 0, -1/8]])
    m.frot_j2x = m.frot_pdf
    m.flavor_basis = ['uv', 'dv', 'g', 'dbar', 'ubar', 'sbar', 'cbar']
    m.evolution_basis = ['S', 'G', 'V', 'V3p', 'V3m', 'V8p', 'V15p']
    ##
    x = 0.1
    Q2 = 1.e4
    xHx = x * m.Hx(x, 0, 0, Q2)
    ## Numbers below are from Pegasus software with IMODEV=3 option.
    assert (xHx[0], xHx[1], xHx[3]-xHx[4], 2*(xHx[3]+xHx[4]),
            2*xHx[-2], 2*xHx[-1], xHx[2]/x) == approx((0.55332,
            0.27275, 0.0099885, 0.39241, 0.11462, 0.060207, 0.90281), abs=1e-3)
