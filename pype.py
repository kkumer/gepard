#!/usr/bin/env python

import shelve

import numpy as np
#np.random.seed(68)

import Model, Approach, Fitter
import utils 
import plots

from results import *

# [1] Load experimental data and theoretical models

data = utils.loaddata('data/ep2epgamma')  # dictionary {1 : DataSet instance, ...}
data.update(utils.loaddata('data/gammastarp2gammap'))
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


## [3] Create a theory

# Gepard only
mGepard = Model.ComptonGepard(cutq2=0.5)
tGepard = Approach.hotfixedBMK(mGepard)

# DR only
mDRonly = Model.ModelDR()
tDR = Approach.hotfixedBMK(mDRonly)
tDR.name = 'DR model'

mDRonly1 = Model.ModelDR()
tDR1 = Approach.hotfixedBMK(mDRonly1)
tDR1.name = 'DR model 1'

# Hybrid: Gepard+DR (can reuse above Gepard)
mDRsea = Model.ComptonModelDRsea()
m = Model.Hybrid(mGepard, mDRsea)
t = Approach.hotfixedBMK(m)
t.name = 'DR + Gepard sea'
g = t.m.g

tDR.m.parameters.update(DMepsGLO)
tDR1.m.parameters.update(DMepsGLO1)


def setpar(i, val):
    mGepard.parameters[mGepard.parameters_index[i]] = val

## [4] Do the fit


#f = Fitter.FitterBrain(fitpoints, t, nnets=12, nbatch=150, verbose=1)
#f.fit()
setpar(11,  0.15203911208796006)
setpar(12,  1.1575060246398083)
#setpar(31,  8.0)
setpar(13,  0.15)
setpar(14,  0.478391)
setpar(15,  0.)
setpar(16,  2.)
setpar(17, -0.15152)
setpar(18,  0.)
setpar(19,  0.)
setpar(22,  1.247316701070471)
#setpar(41,  6.0)
setpar(23,  0.15)
setpar(24,  0.7)
setpar(25,  0.)
setpar(26,  2.)
setpar(27, -0.81217)
setpar(28,  0.)
setpar(29,  0.)
setpar(32,  0.)
setpar(42,  -0.9)
t.m.g.parint.p = 0


t.m.parameters.update(ALTGLO2)

t.m.g.init()

tDR.m.release_parameters('bS', 'rv', 'bv', 'C', 'MC')
fDR = Fitter.FitterMinuit(GLOpoints, tDR)

tDR1.m.release_parameters('bS', 'rv', 'bv', 'C', 'MC', 'trv', 'tbv')
fDR1 = Fitter.FitterMinuit(GLO1points, tDR1)

#t.m.release_parameters('rv', 'bv', 'Nres', 'bres', 'C', 'MC')
#f = Fitter.FitterMinuit(GLOpoints, t)
#f = Fitter.FitterMinuit(BSApoints, t)

#t.m.release_parameters('rv', 'bv', 'C', 'MC', 'trv', 'tbv')
#f = Fitter.FitterMinuit(GLOpoints, t)

#t.m.parameters['tNv'] = 0.6
#t.m.parameters['tMv'] = 0.8
#t.m.parameters['Mv'] = 0.8
t.m.release_parameters('M02S','SECS','SECG', 'THIS', 'THIG', 'rv', 'bv', 'Mv', 'C', 'MC', 'trv', 'tbv', 'tMv')
#t.m.release_parameters('M02S','SECS','SECG', 'THIS', 'THIG', 'rv', 'bv', 'C', 'MC', 'trv', 'tbv')
#f = Fitter.FitterMinuit(DVCSpoints+data[48]+ALTGLO1points, t)
#f = Fitter.FitterMinuit(DVCSpoints+data[48]+ALTGLO2points, t)
f = Fitter.FitterMinuit(DVCSpoints+data[48]+ALTGLO3points, t)
#f = Fitter.FitterMinuit(DVCSpoints+data[48]+GLO1points+HA17[::2], t)

