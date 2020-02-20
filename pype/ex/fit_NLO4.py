#!/usr/bin/env python

import shelve, logging
logging.basicConfig(level=logging.INFO)

import Model, Approach, Fitter, Data, utils, plots
from results import *
from utils import listdb
from abbrevs import *

db = shelve.open('/home/kkumer/pype/theories.db')

thold = db['KM14NLO prelim3']
Model.ComptonGepard.gepardPool.pop()

# Gepard sea part
mGepard = Model.ComptonGepard(ansatz='FIT14', p=1, q02=2.)
Model.ComptonGepard.gepardPool.pop()
thGepard = Approach.BM10(mGepard)
# DR part
mDRsea = Model.ComptonModelDRPPsea()
# Hybridization
m = Model.HybridKelly(mGepard, mDRsea)
th = Approach.BM10(m)
th.name = 'KM14NLO prelim4'
th.description = '20p GLOnoBSS2'
# Initial point:
th.m.parameters.update(thold.m.parameters)

## Do the DIS fit
th.m.fix_parameters('ALL')
th.m.release_parameters('NS', 'NG', 'AL0S', 'AL0G')
f = Fitter.FitterMinuit(DISpoints, th)
f.fit()

th.print_chisq(f.fitpoints)
m.print_parameters_errors()

## Do the preparatory DVCS fit to collider data
th.m.fix_parameters('ALL')
th.m.release_parameters('M02S', 'SECS', 'THIS', 'M02G', 'SECG', 'THIG')
f = Fitter.FitterMinuit(H1ZEUS[::3], th)
f.minuit.tol = 80
#f.minuit.printMode = 1
f.minuit.maxcalls = 1000
f.fit()

th.print_chisq(f.fitpoints)
th.m.print_parameters_errors()

## Do the total DVCS fit
th.m.fix_parameters('ALL')
th.m.release_parameters('Mv', 'rv', 'bv', 'C', 'MC', 
'tMv', 'trv', 'tbv', 'rpi', 'Mpi', 'M02S', 'M02G', 
'SECS', 'THIS', 'SECG', 'THIG',
'ND', 'AL0D', 'ALPD', 'M02D')
f = Fitter.FitterMinuit(GLOnoBSS2+BSSwpoints, th)
f.minuit.tol = 80
f.minuit.printMode = 1
f.minuit.maxcalls = 3000
f.fit()

th.print_chisq(f.fitpoints)
m.print_parameters_errors()
