
import copy, sys
from nose.tools import *
import numpy as np

import utils, Model, Approach, Fitter
from results import DMGLO1, DMepsGLO1, DM12  #use some testpars here?

data = utils.loaddata('data/ep2epgamma', approach=Approach.hotfixedBMK)  
data.update(utils.loaddata('data/en2engamma', approach=Approach.hotfixedBMK))

m = Model.ComptonModelDRFlavored()
m.parameters.update(DMGLO1)
t = Approach.hotfixedBMK(m, optimization = False)

# Gepard model
## mGepard = Model.Gepard(ansatz='EFLEXP')
## th = Approach.BMK(mGepard)
## mGepard.parameters.update(DM12)

pt0 = data[56][0].copy()
pt0.in2particle = 'p'
pt0.phi = 2.*np.pi/5.
pt0.prepare(Approach.hotfixedBMK)

ptn = pt0.copy()
ptn.in2particle = 'n'

ptsNN = data[140][:8]


def test_CFF_Flavored():
    """Calculate u-d separated ImH."""
    assert_almost_equal(m.ImH(pt0, f='u'), 3.207118444556537)
    assert_almost_equal(m.ImH(pt0, f='d'), 1.6035592222782684)
    assert_almost_equal(m.ImH(pt0), 1.6035592222782684)
    assert_almost_equal(m.ImH(ptn), 1.0690394815188455)


def test_CFF_FlavoredDR():
    """Calculate ReH via DR from u-d separated ImH."""
    assert_almost_equal(m.ReH(pt0), -5.243332550866902)

test_CFF_FlavoredDR.long = 1


def test_fitNNFlavored():
    """Testing Flavored Neural Net fitting by FitterBrain."""
    np.random.seed(68)
    mNN = Model.ModelNNKelly(output_layer=['ImH', 'ImHu', 'ImHd',
              'ReH', 'ReHu', 'ReHd'], flavored=['ImH', 'ReH'])
    tNN = Approach.BMK(mNN)
    fNN = Fitter.FitterBrain(ptsNN, tNN, nnets=1, nbatch=30)
    fNN.fitgood()
    chisq = tNN.chisq(ptsNN)[0]
    assert_almost_equal(chisq, 11.184309422965365)

test_fitNNFlavored.long = 1

1.6035592222782684
