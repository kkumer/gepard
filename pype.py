#!/usr/bin/env python

import shelve, shutil

import numpy as np
#np.random.seed(68)

import Model, Approach, Fitter
import utils 
import plots

#from results import *

# [1] Load experimental data and theoretical models

data = utils.loaddata('data/ep2epgamma')  # dictionary {1 : DataSet instance, ...}
data.update(utils.loaddata('data/gammastarp2gammap'))
db = shelve.open('theories.db')

shutil.copy2('test/GEPARD.INI.FIT', 'GEPARD.INI')

## [2] Choose subset of datapoints for fitting

GLOpoints = data[32][12:] + data[8] + data[29]  # DM's GLO set
#GLO1points = data[31][12:] + data[8] + data[29] + data[30]  # DM's GLO1 set
#fitpoints = data[31][12:14] + data[8][1:3] + data[30][2:4]  # test set
#HA17 = utils.select(data[34], criteria=['t == -0.17'])
#HA28 = utils.select(data[34], criteria=['t == -0.28'])
#HA33 = utils.select(data[34], criteria=['t == -0.33'])
#fitpoints = GLO1points + 6*data[30]
#fitpoints = GLOpoints + HA17[::4] + HA33[::4] + data[30]
#fitpoints = GLOpoints
fitpoints = data[36][:4] + data[32][12:16]
DVCSpoints = data[36] + data[37] + data[38] + data[39] + \
  data[40] + data[41] + data[42] + data[43] + data[44] + \
  data[45]

allpoints = GLOpoints + DVCSpoints

#ptSS = HA17[11]

## [3] Create a theory

mGepard = Model.ComptonGepard()
mDR = Model.ComptonModelDR()

m = Model.ComptonGepardDR(mGepard, mDR)
t = Approach.hotfixedBMK(m)

def setpar(i, val):
    mGepard.parameters[mGepard.parameters_index[i]] = val

## [4] Do the fit


#f.fit()

#f = Fitter.FitterBrain(fitpoints, t, nnets=12, nbatch=150, verbose=1)
#f.fit()

setpar(21, 0.5)   
setpar(11, 0.15)  
setpar(12, 1.16)  
setpar(13, 0.15)  
setpar(14, 0.4)   
setpar(16, 2.0)   
setpar(19, -10.0) 
setpar(22, 1.25)  
setpar(23, 0.15)  
setpar(24, 0.5)   
setpar(26, 2.0)   
setpar(29, -32.0) 
setpar(15, 0.0)   
setpar(25, 0.0)   
setpar(17, 0.0)   
setpar(27, 0.0)   
setpar(18, 0.0)   
setpar(28, 0.0)   
t.m.g.parint.p = 0
t.m.g.init()
#t.model.release_parameters('M02S','SKEWS','SKEWG', 'bS', 'Mv')
#f = Fitter.FitterMinuit(allpoints, t)
t.model.release_parameters('M02S','SKEWS', 'Mv')
f = Fitter.FitterMinuit(fitpoints, t)
#f.fit()
