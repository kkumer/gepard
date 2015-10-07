
from nose.tools import *
import shelve, numpy

import utils, Model, Approach, Fitter

data = utils.loaddata('data/ep2epgamma', approach=Approach.hotfixedBMK)  
#db = shelve.open('theories.db')

from abbrevs import *

# data to be fitted to 
pts = HERMES_U + HERMES_L + HERMES_T + CLAS_U + CLAS_L + BSDwpoints + BSSwpoints
hipts = utils.select(pts, criteria=['Q2<=2.45'])  # npts=68

# test-network 
#testNN = db['KMS11-NN']
#testNN.m.useDR = None

numpy.random.seed(38)
mNN = Model.ModelNN(output_layer=['ImH', 'ReH', 'ImE', 'ReE', 'ImHt',
      'ReHt', 'ImEt', 'ReEt'])
tNN = Approach.BM10(mNN)
fNN = Fitter.FitterBrain(hipts, tNN, nnets=8, nbatch=80)
fNN.fit()
tNN.m.parameters['nnet'] = 'ALL'
