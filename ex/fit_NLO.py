#!/usr/bin/env python

import shelve, logging
logging.basicConfig(level=logging.INFO)

import Model, Approach, Fitter, Data, utils, plots
from results import *
from utils import listdb
from abbrevs import *

db = shelve.open('/home/kkumer/pype/theories.db')

#thKMM12 = db['KMM12']
#Model.ComptonGepard.gepardPool.pop()
thprelim1 = db['KM14NLO prelim1']
Model.ComptonGepard.gepardPool.pop()

# Gepard sea part
mGepard = Model.ComptonGepard(p=1, q02=2.)
Model.ComptonGepard.gepardPool.pop()
thGepard = Approach.BM10(mGepard)

## Do the DIS fit
thGepard.model.fix_parameters('ALL')
thGepard.model.release_parameters('NS', 'AL0S', 'AL0G')
f = Fitter.FitterMinuit(DISpoints, thGepard)
f.fit()

thGepard.print_chisq(f.fitpoints)
mGepard.print_parameters_errors()

## Do the preparatory DVCS fit to collider data
thGepard.model.fix_parameters('ALL')
thGepard.model.release_parameters('M02S', 'SECS', 'THIS', 'SECG', 'THIG')
f = Fitter.FitterMinuit(H1ZEUS[::3], thGepard)
f.minuit.tol = 80
f.minuit.printMode = 1
f.minuit.maxcalls = 1000
f.fit()

thGepard.print_chisq(f.fitpoints)
mGepard.print_parameters_errors()

# DR part
mDRsea = Model.ComptonModelDRPPsea()
# Hybridization
m = Model.HybridDipole(mGepard, mDRsea)
th = Approach.BM10(m)
th.name = 'KM14NLO prelim2'
th.description = '15p GLOnoBSS2+BSSwpoints'

## Do the total DVCS fit
th.model.fix_parameters('ALL')
th.model.release_parameters('Mv', 'rv', 'bv', 'C', 'MC', 
'tMv', 'trv', 'tbv', 'rpi', 'Mpi', 'M02S', 'SECS', 'THIS', 'SECG', 'THIG')
# Next line is cheating, of course, to save cca. 1 hour by
# going immediately to minimum
m.parameters.update(thprelim1.m.parameters)
f = Fitter.FitterMinuit(GLOnoBSS2[::2] + BSSwpoints + data[32], th)
f.minuit.tol = 80
f.minuit.printMode = 1
f.minuit.maxcalls = 3000
f.fit()

th.print_chisq(f.fitpoints)
m.print_parameters_errors()

