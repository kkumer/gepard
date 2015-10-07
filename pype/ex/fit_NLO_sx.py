#!/usr/bin/env python

"""Fitting just to small-x DIS and DVCS data."""

import shelve, logging
logging.basicConfig(level=logging.INFO)

import Model, Approach, Fitter, Data, utils, plots
from results import *
from utils import listdb
from abbrevs import *

db = shelve.open('/home/kkumer/pype/theories.db')

#thKMM12 = db['KMM12']
#Model.ComptonGepard.gepardPool.pop()

# Gepard sea part
mGepard = Model.ComptonGepard(p=1, q02=2.0, newgepard=True)
th = Approach.BMK(mGepard)
#th = Approach.BM10(mGepard)

## Do the DIS fit
th.model.fix_parameters('ALL')
th.model.release_parameters('NS', 'AL0S', 'AL0G')
f = Fitter.FitterMinuit(DISpoints, th)
f.fit()

th.print_chisq(f.fitpoints)
mGepard.print_parameters_errors()

## Do the DVCS fit
th.model.fix_parameters('ALL')
th.model.release_parameters('M02S', 'SECS', 'THIS', 'SECG', 'THIG')
f = Fitter.FitterMinuit(H1ZEUS, th)
f.minuit.tol = 80
f.minuit.printMode = 1
f.minuit.maxcalls = 1000
f.fit()


