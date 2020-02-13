
from nose.tools import *
import shelve, numpy

import utils, Model, Approach, Fitter

data = utils.loaddata('data/ep2epgamma', approach=Approach.hotfixedBMK)  
db = shelve.open('theories.db', 'r')
dbndr = shelve.open('nndr.db', 'r')


# testing data set
testpoints = [data[31][12]] + [data[8][1]] + [data[29][2]] + [data[30][3]]
# testing data set for fits
fitpoints = data[31] + data[8]
fitpointsC = data[101][::12] + data[102][::12]

# test-network 
testNN = db['KMS11-NN']
testNN.m.useDR = None

# test-network using dispersion relations
testNNDR = dbndr['NNDR-C15-8']


def test_basicNN():
    """Testing basic Neural Net framework."""
    chisq = 0.
    for pt in testpoints:
        chisq = chisq + (
                (testNN.predict(pt) - pt.val)**2 / pt.err**2 )
    assert_almost_equal(chisq, 20.094357450964957)

test_basicNN.long=1


def test_NNDR():
    """Testing Neural Net + dispersion relations."""
    chisq = testNNDR.chisq(fitpointsC)[0]
    assert_almost_equal(chisq, 10.800606014881728)

test_NNDR.long=1
test_NNDR.extendedtesting=1


def test_fitNN():
    """Testing Neural Net fitting by FitterBrain."""
    numpy.random.seed(81)
    mNN = Model.ModelNN(output_layer=['ImH', 'ReH'])
    tNN = Approach.BMK(mNN)
    fNN = Fitter.FitterBrain(fitpoints, tNN, nnets=1, nbatch=6)
    fNN.fit()
    chisq = tNN.chisq(fitpoints)[0]
    assert_almost_equal(chisq, 34.73003122832)

test_fitNN.long = 1


def test_fitNNDR():
    """Testing Neural Net fitting by FitterBrain (DR+endpoint)."""
    numpy.random.seed(81)
    mNN = Model.ModelNN(output_layer=['ImH', 'ReH'],
            useDR=['ReH'], endpointpower=1)
    tNN = Approach.BMK(mNN)
    fNN = Fitter.FitterBrain(fitpointsC, tNN, nnets=1, nbatch=12)
    fNN.fit()
    chisq = tNN.chisq(fitpointsC)[0]
    assert_almost_equal(chisq, 12.131643143589963)

# Actually this test was broken by Approach.observables consolidation
test_fitNNDR.long = 1

