
from nose.tools import *
import numpy as np

import Model


m = Model.ComptonGepard(ansatz='TEST', scheme='CSBAR')

# Setting gepard to test values
m.g.parint.p = 1
m.g.parint.nf = 4
m.g.astrong.asp = np.array([0.05, 0.05, 0.05])
m.g.parflt.q02 = 1.0

nf = 4
z = 0.1+2.3j  # testpoint

def test_WgammaLO():
    """Calculate anomalous dimensions up to NNLO."""
    m.g.parint.p = 1
    m.g.harmonic.s1 = m.g.hs1(z)
    m.g.harmonic.s2 = m.g.hs2(z)
    m.g.harmonic.s3 = m.g.hs3(z)
    m.g.harmonic.s4 = m.g.hs4(z)
    qq0 = m.g.wgammavqq0f(nf, z)
    gg0 = m.g.wgammavgg0f(nf, z)
    qq1 = m.g.wgammavqq1f(nf, z)
    qq2 = m.g.wgammavqq2f(nf, z)
    nsp1 = m.g.wgammavnsp1f(nf, z)
    nsm1 = m.g.wgammavnsm1f(nf, z)
    assert_almost_equal(qq0, 4.0622581815675+7.2098651918316j)
    assert_almost_equal(gg0/100., (11.4686+16.4593j)/100., 5)
    assert_almost_equal(qq1/100., (22.6712+19.7996j)/100., 5)
    assert_almost_equal(qq2/100., (84.4257+87.0365j)/100., 5)
    assert_almost_equal(nsp1/100., (24.1347+20.6634j)/100., 6)
    assert_almost_equal(nsm1/100., (24.1+20.714j)/100., 6)


