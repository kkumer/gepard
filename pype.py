#!/usr/bin/env python

import shelve

import utils 
import models
import Approach
import fit
import plots

from constants import *
from results import *

# [1] Load experimental data and theoretical models

data = utils.loaddata('data/ep2epgamma')  # dictionary {1 : DataSet instance, ...}
db = shelve.open('theories.db')

# some shortcuts
DMGLO = db['DMGLO']
DMGLO1 = db['DMGLO1']
KKGLO = db['KKGLO']
KKGLO1 = db['KKGLO1']
NN1 = db['NN1']


## [3] Choose subset of datapoints for fitting

GLOpoints = data[32][12:] + data[8] + data[29]  # DM's GLO set
GLO1points = data[31][12:] + data[8] + data[29] + data[30]  # DM's GLO1 set
fitpoints = data[31][12:14] + data[8][1:3] + data[30][2:4]  # test set

## [4] Do the fit

#t = DMGLO1.copy()
#t.model.release_parameters('bS', 'rv', 'bv', 'C', 'MC', 'trv')
#f = fit.FitterMinuit(GLOpoints, t)
#newDMGLO = f.fit()

m = models.ModelNN()
t = Approach.hotfixedBMK(m)
f = fit.FitterBrain(data[5], t)

