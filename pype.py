#!/usr/bin/env python

import shelve, shutil

import numpy as np
#np.random.seed(68)

import Model, Approach, Fitter
import utils 
import plots

#from results import *

# [1] Load experimental data and theoretical models

#data = utils.loaddata('data/ep2epgamma')  # dictionary {1 : DataSet instance, ...}
data = utils.loaddata('data/gammastarp2gammap')  
#db = shelve.open('theories.db')

shutil.copy2('test/GEPARD.INI.FIT', 'GEPARD.INI')

## [2] Choose subset of datapoints for fitting

#GLOpoints = data[32][12:] + data[8] + data[29]  # DM's GLO set
#GLO1points = data[31][12:] + data[8] + data[29] + data[30]  # DM's GLO1 set
#fitpoints = data[31][12:14] + data[8][1:3] + data[30][2:4]  # test set
#HA17 = utils.select(data[34], criteria=['t == -0.17'])
#HA28 = utils.select(data[34], criteria=['t == -0.28'])
#HA33 = utils.select(data[34], criteria=['t == -0.33'])
#fitpoints = GLO1points + 6*data[30]
#fitpoints = GLOpoints + HA17[::4] + HA33[::4] + data[30]
#fitpoints = GLOpoints
fitpoints = data[36]

#ptSS = HA17[11]

## [3] Create a theory

mGepard = Model.ComptonGepard()
mDR = Model.ComptonModelDR()

m = Model.ComptonGepardDR(mGepard, mDR)
t = Approach.hotfixedBMK(m)

t.m.g.parint.p = 0
t.m.g.init()

## [4] Do the fit

t.model.release_parameters('M02S', 'SKEWS')
f = Fitter.FitterMinuit(fitpoints, t)

#f.fit()

#f = Fitter.FitterBrain(fitpoints, t, nnets=12, nbatch=150, verbose=1)
#f.fit()
