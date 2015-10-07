
from nose.tools import *
import shelve, numpy

import utils, Model, Approach, Fitter

data = utils.loaddata('data/ep2epgamma', approach=Approach.hotfixedBMK)  
#db = shelve.open('theories.db')

from abbrevs import *

# data to be fitted to 
pts = HERMES_U + HERMES_L + HERMES_T + CLAS_U + CLAS_L
lowpts = utils.select(pts, criteria=['Q2<=2.45'])  # npts=48

# test-network 
#testNN = db['KMS11-NN']
#testNN.m.useDR = None

numpy.random.seed(68)
mNN = Model.ModelNN(output_layer=['ImH', 'ReH', 'ImE', 'ReE', 'ImHt',
      'ReHt', 'ImEt', 'ReEt'])
tNN = Approach.BM10(mNN)
fNN = Fitter.FitterBrain(lowpts, tNN, nnets=8, nbatch=150)
fNN.fit()
tNN.m.parameters['nnet'] = 'ALL'
