#!/usr/bin/env python

import shelve

import utils 
import models
import Approach
import fit
import plots

from constants import *
from results import *

# [1] Choose theoretical approach

mDR = models.ModelDR()
mNN = models.ModelNN()
#bmk = Approach.BMK()
bDR = Approach.hotfixedBMK(mDR, optimization = False)  # no optimizations
#bNN = Approach.hotfixedBMK(mNN, optimization = False)  # no optimizations

# [2] Load data and models

data = utils.loaddata('data/ep2epgamma', approach=bDR)  # dictionary {1 : DataSet instance, ...}
db = shelve.open('models.db')

## [3] Choose subset of datapoints for fitting

GLOpoints = data[31][12:] + data[8] + data[29]  # DM's GLO set
GLO1points = data[31][12:] + data[8] + data[29] + data[30]  # DM's GLO1 set
fitpoints = data[31][12:14] + data[8][1:3] + data[30][2:4]  # test set

## [4] Do the fit

#ff.parameters.update(DMGLO1)
#ff.release_parameters('bS', 'Mv')
#f = fit.FitterMinuit(fitpoints, b, ff)
##toy = f.fit()

