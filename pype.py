#!/usr/bin/env python

import sys, os

try: # if you have ROOT you might want minuit2
    from minuit2 import Minuit2 as Minuit
except:
    from minuit import Minuit

from scipy.special import gammainc

import plots
import utils 
import models
import Approach
import fit

from constants import *
from results import *

# [1] Load data

data = utils.loaddata('data/ep2epgamma')   # dictionary {1 : DataSet instance, ...}

# [2] Choose datapoints for fitting

fitpoints = data[31][12:] + data[8] + data[29]  # DM's GLO set
#fitpoints = data[31][12:] + data[8] + data[29] + data[30]  # DM's GLO1 set
#fitpoints = [data[31][12]] + [data[8][1]] + [data[29][2]] + [data[30][3]]  # test set


# [3] Choose theoretical approach

ff = models.ModelDR()
#ff = models.ModelNN()
#bmk = Approach.BMK()
bnoopt = Approach.hotfixedBMK(ff, optimization = False)  # no optimizations
bopt = Approach.hotfixedBMK(ff, optimization = True)     # optimized formulas
b = bnoopt

# Transform data into conventions of approach used, and 
# pre-calculate what can be pre-calculated

## Prepare just fitpoints ...
# [pt.to_conventions(b) for pt in fitpoints]
# [pt.prepare(b) for pt in fitpoints]
## ... or prepare ALL available datapoints
[[pt.to_conventions(b) for pt in set] for set in data.values()]
[[pt.prepare(b) for pt in set] for set in data.values()]


f = fit.FitterMinuit(fitpoints, b, ff)
#f.fit()
f.printres()
