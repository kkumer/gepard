#!/usr/bin/env python

import pickle

import utils 
import models
import Approach
import fit
import plots

from constants import *
from results import *

# [1] Choose theoretical approach

#ff = models.ModelDR()
ff = models.ModelNN()
#bmk = Approach.BMK()
b = Approach.hotfixedBMK(ff, optimization = False)  # no optimizations

# [2] Load data

data = utils.loaddata('data/ep2epgamma', approach=b)  # dictionary {1 : DataSet instance, ...}

## [3] Choose subset of datapoints for fitting

#GLOpoints = data[31][12:] + data[8] + data[29]  # DM's GLO set
#GLO1points = data[31][12:] + data[8] + data[29] + data[30]  # DM's GLO1 set
#fitpoints = data[31][12:14] + data[8][1:3] + data[30][2:4]  # test set

## [4] Do the fit

#ff.parameter_dict.update(DMGLO1)
#ff.release_parameters('bS', 'Mv')
#f = fit.FitterMinuit(fitpoints, b, ff)
##toy = f.fit()

## Load some results
#KKGLO = pickle.load(open('KKGLO.pkl', 'r'))
#KKGLO1 = pickle.load(open('KKGLO1.pkl', 'r'))
