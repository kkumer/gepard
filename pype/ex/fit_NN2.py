
from nose.tools import *
import shelve, numpy

import utils, Model, Approach, Fitter

data = utils.loaddata('data/ep2epgamma', approach=Approach.hotfixedBMK)  
#db = shelve.open('theories.db')

from abbrevs import *

# data to be fitted to (just HERMES BSA+BCA):
# pts = ( L4_ALUI+ L4_AC_0+ L4_AC_1+ L4_ALL_1+ L4_AUTI_1+ L4_AUTI_0+ L4_AUTI_m1+
# L4_AUTDVCS+ L4_ALTI_1+ L4_ALTI_0+ L4_ALTI_m1+ L4_ALTBHDVCS_0)
# data to be fitted to (total Hall A)
#HApts =  data[116]+data[117]+data[136]+data[135]
#HApts =  data[136]+data[135]
# highQ = utils.select(data[136], ['Q2 > 1.8'])
# highQD = utils.select(data[135], ['Q2 > 1.8'])
cos0 = utils.select(data[136], ['FTn == 0'])
cos1 = utils.select(data[136], ['FTn == 1'])
#pts =  data[116]+data[117]
#pts = C_BSSw0pts[::3] + C_BSSw1pts[::3] + C_BSDwpts[::3]
#pts = utils.select(data[135], ['Q2 < 1.7'])
pts =  cos0
# test-network 
#testNN = db['KMS11-NN']
#testNN.m.useDR = None

numpy.random.seed(68)
mNN = Model.ModelNN(output_layer=['ImH', 'ReH', 'ImHt', 'ReHt', 'ImE', 'ReE', 'ImEt', 'ReEt',
     'ImHeff', 'ReHeff', 'ImHteff', 'ReHteff', 'ImEeff', 'ReEeff', 'ImEteff', 'ReEteff'])
tNN = Approach.BM10(mNN)
fNN = Fitter.FitterBrain(pts, tNN, nnets=10, nbatch=300, batchlen=5)
fNN.fit()
tNN.m.parameters['nnet'] = 'ALL'

