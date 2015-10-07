
from nose.tools import *
import shelve, numpy

import utils, Model, Approach, Fitter

data = utils.loaddata('data/ep2epgamma', approach=Approach.hotfixedBMK)  
db = shelve.open('theories.db')

from abbrevs import *

# data to be fitted to (just HERMES BSA+BCA):
pts = ( L4_ALUI+ L4_AC_0+ L4_AC_1+ L4_ALL_1+ L4_AUTI_1+ L4_AUTI_0+ L4_AUTI_m1+
L4_AUTDVCS+ L4_ALTI_1+ L4_ALTI_0+ L4_ALTI_m1+ L4_ALTBHDVCS_0)

# test-network 
#testNN = db['KMS11-NN']
#testNN.m.useDR = None

#numpy.random.seed(68)
mNN = Model.ModelNN(output_layer=['ImH', 'ReE'])
tNN = Approach.BM10(mNN)
fNN = Fitter.FitterBrain(GLOfixfull, tNN, nnets=50, nbatch=50)
fNN.fit()
tNN.m.parameters['nnet'] = 'ALL'

mNN2 = Model.ModelNN(output_layer=['ImH', 'ReE'])
tNN2 = Approach.BM10(mNN2)
fNN2 = Fitter.FitterBrain(GLOfixfullKK, tNN2, nnets=50, nbatch=50)
fNN2.fit()
tNN2.m.parameters['nnet'] = 'ALL'
