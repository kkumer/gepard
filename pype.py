#!/usr/bin/env python

import shelve

import numpy as np

import Model, Approach, Fitter, Data
import utils 
import plots

from constants import Mp, Mp2

from results import *

# [1] Load experimental data and theoretical models

data = utils.loaddata('data/ep2epgamma', approach=Approach.hotfixedBMK)  
#data.update(utils.loaddata('data/gammastarp2gammap', approach=Approach.hotfixedBMK))
#db = shelve.open('theories.db')

## [2] Choose subset of datapoints for fitting

testpoints = data[31][12:14] + data[8][1:3] + data[30][2:4]  # test set
GLOpoints = data[31][12:] + data[8] + data[29]  # DM's GLO set
GLO1points = data[31][12:] + data[8] + data[29] + data[30]  # DM's GLO1 set
#HERMESpoints = data[31][12:] +  data[29]
#BSApoints = data[8] + data[29]
#HAD17 = utils.select(data[33], criteria=['Q2 == 1.5', 't == -0.17'])
#HA17 = utils.select(data[34], criteria=['t == -0.17'])
#HA28 = utils.select(data[34], criteria=['t == -0.28'])
#HA33 = utils.select(data[34], criteria=['t == -0.33'])
#DVCSpoints = data[36] + data[37] + data[38] + data[39] + \
#  data[40] + data[41] + data[42] + data[43] + data[44] + \
#  data[45]
#ALTGLOpoints = data[5] + data[25] + data[32][18:]
#ALTGLO1points = data[5] + data[25] + data[32] + HAD17 + HA17
#ALTGLO2points = data[5] + data[25] + data[32][18:] + HAD17[::2] + HA17[::2]
#ALTGLO3points = data[5] + data[25] + data[32][18:] + data[30] + HA17
#ALTGLO4points = data[25] + data[32][18:]
#BSDw2Cpoints = utils.select(data[26], criteria=['Q2 == 2.3'])
#BSDw2CDpoints = utils.select(data[50], criteria=['Q2 == 2.3'])
#BSSwpoints = utils.select(data[51], criteria=['FTn != 2'])
TSA1points = utils.select(data[52], criteria=['FTn == -1'])



## [3] Create a theory

# Gepard only
#mGepard = Model.ComptonGepard(cutq2=0.5)
#tGepard = Approach.hotfixedBMK(mGepard)


# DR only
mDRonly = Model.ModelDR()
tDR = Approach.hotfixedBMK(mDRonly)
tDR.name = 'DR model GLO'
tDR.m.parameters.update(DMepsGLO)

mDRonly1 = Model.ModelDR()
tDR1 = Approach.hotfixedBMK(mDRonly1)
tDR1.name = 'DR model GLO1'
tDR1.m.parameters.update(DMepsGLO1)


## Hybrid: Gepard+DR (can reuse above Gepard)
#mDRsea = Model.ComptonModelDRsea()
#m = Model.Hybrid(mGepard, mDRsea)
#t = Approach.hotfixedBMK(m)
#t.name = 'DR + Gepard sea'
#g = t.m.g

mDRonly2 = Model.ModelDR()
t = Approach.BM10(mDRonly2)
t.name = 'BM10'
t.m.parameters.update(DMepsGLO1)

#t.m.parameters.update(hy1THI)

# NN
#mNN = Model.ModelNN(hidden_layers=[15], output_layer=['ImH', 'ReH', 'ImE', 'ReE', 'ImHt', 'ReHt', 'ImEt', 'ReEt'])
#mNN = Model.ModelNN(hidden_layers=[9], endpointpower=3.0)
#tNN = Approach.hotfixedBMK(mNN)
#tNN.name = 'NNtest'
#tNN.description = 'x (xB,t)-13-8 nets trained on ALTGLO+BSDw2CD+BSSw for 500 batches'

## [4] Do the fit

#f = Fitter.FitterBrain(BSDw2CDpoints, tNN, nnets=20, nbatch=50, verbose=1)
#f.fit()
#f.prune(minprob=0.5)
#tNN.save(db)
#db.close()

## Fitting to both small- and large-x data
#t.m.release_parameters('M02S','SECG', 'THIS', 'THIG', 'rv', 'bv', 'Mv', 'C', 'MC', 'trv', 'tbv', 'tMv')
#f = Fitter.FitterMinuit(DVCSpoints+data[48]+ALTGLO2points, t)

## DR fit
#t.m.release_parameters('rv', 'bv', 'Mv', 'C', 'MC', 'trv', 'tbv', 'tMv')
#f = Fitter.FitterMinuit(GLOpoints, t)


