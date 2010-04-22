#!/usr/bin/env python

import shelve

import numpy as np
#np.random.seed(68)

import Model, Approach, Fitter
import utils 
import plots

#from results import *

# [1] Load experimental data and theoretical models

data = utils.loaddata('data/ep2epgamma')  # dictionary {1 : DataSet instance, ...}
db = shelve.open('theories.db')


## [2] Choose subset of datapoints for fitting

GLOpoints = data[32][12:] + data[8] + data[29]  # DM's GLO set
GLO1points = data[31][12:] + data[8] + data[29] + data[30]  # DM's GLO1 set
fitpoints = data[31][12:14] + data[8][1:3] + data[30][2:4]  # test set
HA17 = utils.select(data[34], criteria=['t == -0.17'])
HA28 = utils.select(data[34], criteria=['t == -0.28'])
HA33 = utils.select(data[34], criteria=['t == -0.33'])
#fitpoints = GLO1points + 6*data[30]
fitpoints = GLOpoints + HA17[::4] + HA33[::4] + data[30]
#fitpoints = data[26]
#fitpoints = GLOpoints
ptSS = HA17[11]

## [3] Create a theory

m = Model.ModelNN()
m = Model.ModelNN(hidden_layers = [11, 9], output_layer=['ImH', 'ReH', 'ImHt', 'ReHt', 'ImE', 'ReE', 'ImEt', 'ReEt'])
#m = Model.ModelNN(output_layer=['ImH', 'ReH', 'ImHt', 'ReHt', 'ImE', 'ReE'])
#m = Model.ModelNN(output_layer=['ImH'])
t = Approach.hotfixedBMK(m)

#m = Model.ModelDR()
#m.parameters.update(DMGLO)
#t = Approach.hotfixedBMK(m)
#t.name = 'DMGLO'
#t.description = 'DM fit to HERMES and CLAS BSA/BCA data. Only CFF H.'
#t.save(db)
#del m, t
#
#m = Model.ModelDR()
#m.parameters.update(DMGLO1)
#t = Approach.hotfixedBMK(m)
#t.name = 'DMGLO1'
#t.description = 'DM fit to HERMES and CLAS BSA/BCA data, and HALL-A. With large CFF Htilde'
#t.save(db)
#db['DMGLO1'] = t
#del m, t
#
#db.close()

## [4] Do the fit

#t = db['DMGLO1'].copy()
#t.model.release_parameters('bS', 'rv', 'bv', 'C', 'MC', 'trv')
#f = Fitter.FitterMinuit(fitpoints, t)

f = Fitter.FitterBrain(fitpoints, t, nnets=6, nbatch=50, verbose=1)
#f.fit()
