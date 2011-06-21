
from nose.tools import *
import shelve, numpy

import utils, Model, Approach, Fitter

data = utils.loaddata('data/ep2epgamma', approach=Approach.hotfixedBMK)  
db = shelve.open('theories.db')


# testing data set
testpoints = [data[31][12]] + [data[8][1]] + [data[29][2]] + [data[30][3]]
# testing data set for fits
fitpoints = data[31][12:14] + data[8][1:3] + data[30][2:4]

# test-network 
testNN = db['KMS11-NN']


def test_basicNN():
    """Testing basic Neural Net framework."""
    chisq = 0.
    for pt in testpoints:
        chisq = chisq + (
                (testNN.predict(pt) - pt.val)**2 / pt.err**2 )
    assert_almost_equal(chisq, 20.094357450964957)


def test_fitNN():
    """Testing Neural Net fitting by FitterBrain."""
    numpy.random.seed(68)
    mNN = Model.ModelNN()
    tNN = Approach.hotfixedBMK(mNN)
    fNN = Fitter.FitterBrain(2*fitpoints, tNN, nnets=1, nbatch=6)
    fNN.fit()
    chisq = tNN.chisq(2*fitpoints)[0]
    assert_almost_equal(chisq, 22.754514666801569)

#test_fitNN.long = 1

db.close()
