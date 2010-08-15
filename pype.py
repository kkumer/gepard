#!/usr/bin/env python

import shelve

import numpy as np

import Model, Approach, Fitter
import utils 
import plots

import Data
from constants import Mp, Mp2

from results import *

# [1] Load experimental data and theoretical models

data = utils.loaddata('data/ep2epgamma', approach=Approach.hotfixedBMK)  
data.update(utils.loaddata('data/gammastarp2gammap', approach=Approach.hotfixedBMK))
db = shelve.open('theories.db')



## [2] Choose subset of datapoints for fitting

testpoints = data[31][12:14] + data[8][1:3] + data[30][2:4]  # test set
GLOpoints = data[31][12:] + data[8] + data[29]  # DM's GLO set
GLO1points = data[31][12:] + data[8] + data[29] + data[30]  # DM's GLO1 set
HERMESpoints = data[31][12:] +  data[29]
BSApoints = data[8] + data[29]
HAD17 = utils.select(data[33], criteria=['Q2 == 1.5', 't == -0.17'])
HA17 = utils.select(data[34], criteria=['t == -0.17'])
HA28 = utils.select(data[34], criteria=['t == -0.28'])
HA33 = utils.select(data[34], criteria=['t == -0.33'])
DVCSpoints = data[36] + data[37] + data[38] + data[39] + \
  data[40] + data[41] + data[42] + data[43] + data[44] + \
  data[45]
ALTGLOpoints = data[5] + data[25] + data[32][18:]
ALTGLO1points = data[5] + data[25] + data[32] + HAD17 + HA17
ALTGLO2points = data[5] + data[25] + data[32][18:] + HAD17[::2] + HA17[::2]
ALTGLO3points = data[5] + data[25] + data[32][18:] + data[30] + HA17
ALTGLO4points = data[25] + data[32][18:]


## [3] Create a theory

# Gepard only
#mGepard = Model.ComptonGepard(cutq2=0.5)
#tGepard = Approach.hotfixedBMK(mGepard)


# DR only
mDRonly = Model.ModelDR()
tDR = Approach.hotfixedBMK(mDRonly)
tDR.name = 'DR model'

mDRonly1 = Model.ModelDR()
tDR1 = Approach.hotfixedBMK(mDRonly1)
tDR1.name = 'DR model 1'


# Hybrid: Gepard+DR (can reuse above Gepard)
#mDRsea = Model.ComptonModelDRsea()
#m = Model.Hybrid(mGepard, mDRsea)
#t = Approach.hotfixedBMK(m)
#t.name = 'DR + Gepard sea'
#g = t.m.g


tDR.m.parameters.update(DMepsGLO)
tDR1.m.parameters.update(DMepsGLO1)

#t.m.parameters.update(allTHI)

# NN
#mNN = Model.ModelNN(hidden_layers=[26], output_layer=['ImH', 'ReH', 'ImE', 'ReE', 'ImHt', 'ReHt', 'ImEt', 'ReEt'])
#mNN = Model.ModelNN(hidden_layers=[9], endpointpower=3.0)
#tNN = Approach.hotfixedBMK(mNN)
#tNN.name = 'NNtest'
#tNN.description = '(xB,t)-9-2 nets. Fit to ALTGLO4points'

## [4] Do the fit


#f = Fitter.FitterBrain(ALTGLO4points, tNN, nnets=30, nbatch=20, verbose=1)
#f.fit()
#f.prune(minprob=0.5)
#tNN.save(db)
#db.close()

#t.m.release_parameters('M02S','SECG', 'THIS', 'THIG', 'rv', 'bv', 'Mv', 'C', 'MC', 'trv', 'tbv', 'tMv')
#f = Fitter.FitterMinuit(DVCSpoints+data[48]+ALTGLO2points, t)

