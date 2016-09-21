
from nose.tools import *
import shelve, numpy

import utils, Model, Approach, Fitter

data = utils.loaddata('data/ep2epgamma', approach=Approach.hotfixedBMK)  
db = shelve.open('theories.db')


# testing data set
testpoints = [data[31][12]] + [data[8][1]] + [data[29][2]] + [data[30][3]]
# testing data set for fits
fitpoints = data[31] + data[8]

# test-network 
testNN = db['KMS11-NN']
testNN.m.useDR = None


def test_basicNN():
    """Testing basic Neural Net framework."""
    chisq = 0.
    for pt in testpoints:
        chisq = chisq + (
                (testNN.predict(pt) - pt.val)**2 / pt.err**2 )
    assert_almost_equal(chisq, 20.094357450964957)

test_basicNN.long=1

def test_fitNN():
    """Testing Neural Net fitting by FitterBrain."""
    numpy.random.seed(81)
    mNN = Model.ModelNN(output_layer=['ImH', 'ReH'])
    tNN = Approach.BMK(mNN)
    fNN = Fitter.FitterBrain(fitpoints, tNN, nnets=1, nbatch=6)
    fNN.fit()
    chisq = tNN.chisq(fitpoints)[0]
    assert_almost_equal(chisq, 34.73003122832)

# Actually this test was broken by Approach.observables consolidation
test_fitNN.long = 1

db.close()
