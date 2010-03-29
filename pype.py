#!/usr/bin/env python

import utils 
import models
import Approach
import fit

from constants import *
from results import *

# [1] Choose theoretical approach

ff = models.ModelDR()
#ff = models.ModelNN()
#bmk = Approach.BMK()
b = Approach.hotfixedBMK(ff, optimization = False)  # no optimizations

# [2] Load data

data = utils.loaddata('data/ep2epgamma', approach=b)  # dictionary {1 : DataSet instance, ...}

# [3] Choose subset of datapoints for fitting

#fitpoints = data[31][12:] + data[8] + data[29]  # DM's GLO set
#fitpoints = data[31][12:] + data[8] + data[29] + data[30]  # DM's GLO1 set
fitpoints = data[31][12:14] + data[8][1:3] + data[30][2:4]  # test set

# [4] Do the fit

ff.parameter_dict.update(DMGLO1)
ff.parameter_dict.update({'fix_C':True, 'fix_MC':True})
f = fit.FitterMinuit(fitpoints, b, ff)
#f.fit()
f.printres()
