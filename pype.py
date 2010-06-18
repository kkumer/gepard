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

GLOpoints = data[31][12:] + data[8] + data[29]  # DM's GLO set
GLO1points = data[31][12:] + data[8] + data[29] + data[30]  # DM's GLO1 set
testpoints = data[31][12:14] + data[8][1:3] + data[30][2:4]  # test set
#HA17 = utils.select(data[34], criteria=['t == -0.17'])
#HA28 = utils.select(data[34], criteria=['t == -0.28'])
#HA33 = utils.select(data[34], criteria=['t == -0.33'])
DVCSpoints = data[36] + data[37] + data[38] + data[39] + \
  data[40] + data[41] + data[42] + data[43] + data[44] + \
  data[45]

testpoints = GLOpoints[-6:]
testpoints = GLOpoints[16:24]
allpoints = GLOpoints + DVCSpoints

#ptSS = HA17[11]
ptt = testpoints[1]
pt0 = GLOpoints[0]
pt4 = GLOpoints[4]
pt5 = GLOpoints[5]
pt9 = GLOpoints[9]
pt26 = GLOpoints[26]

## [3] Create a theory

# Gepard only
mGepard = Model.ComptonGepard()
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

#tDR.m.parameters.update(DMGLO)
#tDR1.m.parameters.update(DMGLO1)

t.m.parameters.update(DMepsGLO)

#tDR.m.parameters['NS'] = 0.6
#tDR.m.parameters['alS'] = 1.25
#tDR.m.parameters['MS'] = 0.846

#tDR.m.parameters['Mv'] = 1.5


def setpar(i, val):
    mGepard.parameters[mGepard.parameters_index[i]] = val

## [4] Do the fit


#f = Fitter.FitterBrain(fitpoints, t, nnets=12, nbatch=150, verbose=1)
#f.fit()
setpar(11,  0.15203911208796006)
setpar(12,  1.1575060246398083)
setpar(13,  0.15)
setpar(14,  0.41768649210641556)
setpar(15,  0.)
setpar(16,  2.)
setpar(17,  0.)
setpar(18,  0.)
setpar(19,  -10.5)
setpar(22,  1.247316701070471)
setpar(23,  0.15)
setpar(24,  0.7)
setpar(25,  0.)
setpar(26,  2.)
setpar(27,  0.)
setpar(28,  0.)
setpar(29,  -31.89)
t.m.g.parint.p = 0

#setpar(11, 0.1678)
#setpar(12, 1.12835)
#setpar(14, 0.565539)
#setpar(19, 0.)
#setpar(22, 1.099)
#setpar(24, 0.7)
#setpar(29, 0.)
#t.m.g.parint.p = 1



t.m.g.init()

tDR.m.release_parameters('bS', 'rv', 'bv', 'C', 'MC')
fDR = Fitter.FitterMinuit(GLOpoints, tDR)

tDR1.m.release_parameters('bS', 'rv', 'bv', 'C', 'MC', 'trv', 'tbv')
fDR1 = Fitter.FitterMinuit(GLO1points, tDR1)

#t.m.parameters['bv'] = 0.5
#t.m.release_parameters('rv', 'Mv', 'C', 'MC', 'tNv')
#f = Fitter.FitterMinuit(GLOpoints, t)

