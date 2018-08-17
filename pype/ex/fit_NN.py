
import sys, shelve, logging
logging.basicConfig(level=logging.ERROR)
GEPARD_DIR = '/home/kkumer/gepard'
sys.path.append(GEPARD_DIR+'/pype')
import Model, Approach, Fitter, Data, utils, plots
from results import *
from abbrevs import *

HEpts = ( L4_ALUI+ L4_AC_0+ L4_AC_1+ L4_ALL_1+ L4_AUTI_1+ L4_AUTI_0+ L4_AUTI_m1+
        L4_AUTDVCS+ L4_ALTI_1+ L4_ALTI_0+ L4_ALTI_m1+ L4_ALTBHDVCS_0)  # HERMES as in KMM12

CLAS14pts = CLAS14TSApts + CLAS14BTSApts + CLAS14BSApts
CLAS15pts = C_BSDwpts + C_BSSw0pts + C_BSSw1pts
HA15pts = H_BSDwpts + H_BSSw0pts + H_BSSw1pts
HA17pts = H17_BSDwpts + H17_BSSw0pts + H17_BSSw1pts

pts = CLAS15pts[::4]
#numpy.random.seed(68)

numn = 6

# Investigate bias-variance
#  1. vary the number of CFFs
#  2. vary the hidden layers

mNN = Model.ModelNN(hidden_layers = [3], output_layer=['ImH', 'ReH'])
tNN = Approach.BM10tw2(mNN)
tNN.name = "NN-H3"
fNN = Fitter.FitterBrain(pts, tNN, nnets=numn, nbatch=30)
fNN.fit()
tNN.m.parameters['nnet'] = 'ALL'
fNN.prune(minprob=0.1)
tNN.description = 'CLAS15 2-3-2 {:.1f}/{} p={:.3g}'.format(*tNN.chisq(pts))


mNN2 = Model.ModelNN(hidden_layers = [7], output_layer=['ImH', 'ReH'])
tNN2 = Approach.BM10tw2(mNN2)
tNN2.name = "NN-H7"
fNN2 = Fitter.FitterBrain(pts, tNN2, nnets=numn, nbatch=30)
fNN2.fit()
tNN2.m.parameters['nnet'] = 'ALL'
fNN2.prune(minprob=0.1)
tNN2.description = 'CLAS15 2-7-2 {:.1f}/{} p={:.3g}'.format(*tNN2.chisq(pts))

mNN3 = Model.ModelNN(hidden_layers = [13], output_layer=['ImH', 'ReH'])
tNN3 = Approach.BM10tw2(mNN3)
tNN3.name = "NN-H13"
fNN3 = Fitter.FitterBrain(pts, tNN3, nnets=numn, nbatch=30)
fNN3.fit()
tNN3.m.parameters['nnet'] = 'ALL'
fNN3.prune(minprob=0.1)
tNN3.description = 'CLAS15 2-13-2 {:.1f}/{} p={:.3g}'.format(*tNN3.chisq(pts))

mNN4 = Model.ModelNN(hidden_layers = [21], output_layer=['ImH', 'ReH'])
tNN4 = Approach.BM10tw2(mNN4)
tNN4.name = "NN-H21"
fNN4 = Fitter.FitterBrain(pts, tNN4, nnets=numn, nbatch=30)
fNN4.fit()
tNN4.m.parameters['nnet'] = 'ALL'
fNN4.prune(minprob=0.1)
tNN4.description = 'CLAS15 2-21-2 {:.1f}/{} p={:.3g}'.format(*tNN4.chisq(pts))

mNN5 = Model.ModelNN(hidden_layers = [7, 11], output_layer=['ImH', 'ReH'])
tNN5 = Approach.BM10tw2(mNN5)
tNN5.name = "NN-H7-11"
fNN5 = Fitter.FitterBrain(pts, tNN5, nnets=numn, nbatch=30)
fNN5.fit()
tNN5.m.parameters['nnet'] = 'ALL'
fNN5.prune(minprob=0.1)
tNN5.description = 'CLAS15 2-7-11-2 {:.1f}/{} p={:.3g}'.format(*tNN5.chisq(pts))
