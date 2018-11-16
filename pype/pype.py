#!/usr/bin/env python

import shelve, copy, sys, logging, builtins
import numpy as np

logging.basicConfig(level=logging.WARNING)
#logging.basicConfig(level=logging.INFO)
#logging.basicConfig(level=logging.DEBUG)

import Model, Approach, Fitter, Data, utils, plots
from results import *
from utils import listdb
from abbrevs import *

db = shelve.open('theories.db')

th = db['KM15']
Model.ComptonGepard.gepardPool.pop()
#thb = db['KMM12']
#Model.ComptonGepard.gepardPool.pop()

#mGK = Model.GK12()
#thGK = Approach.BM10(mGK)
#thGK.name = 'GK12'


# Some combination of datasets used for global fitting
pts = GLO15b

utils.describe_data(pts)

# Print chi-square of model w.r.t. this dataset
th.print_chisq(pts)


## --- New fit ---

## --- Creating new model: ---
## Gepard sea part
#mGepard = Model.ComptonGepard(p=0, q02=4.0)
#Model.ComptonGepard.gepardPool.pop()
#thGepard = Approach.BM10(mGepard)
## DR part
#mDRsea = Model.ComptonModelDRPPsea()
## Hybridization
#m = Model.HybridDipole(mGepard, mDRsea)
#th = Approach.BM10(m)
#th.model.fix_parameters('ALL')
#th.model.release_parameters(
#   'ALPS', 'M03S', 'SECS', 'THIS', 'ALPG', 'M02G', 'SECG', 'THIG')

#datcut = utils.select(H109XL+H109WdepXL+H1ZEUScut, criteria=['Q2 >= 4.0'])
#f = Fitter.FitterMinuit(datcut, th)
#f.minuit.tol = 80
#f.minuit.printMode = 1
#f.minuit.maxcalls = 1000

# allpars = ['Mv', 'rv', 'bv', 'C', 'MC', 
# 'tMv', 'trv', 'tbv', 'rpi', 'Mpi', 'M02S', 'SECS', 'THIS', 'SECG', 'THIG']

#f.fit()
