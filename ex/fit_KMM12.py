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

# Gepard sea part
mGepard = Model.ComptonGepard(p=0, newgepard=True)
thGepard = Approach.BM10(mGepard)
# DR part
mDRsea = Model.ComptonModelDRPPsea()
# Hybridization
m = Model.HybridDipole(mGepard, mDRsea)
th = Approach.BM10(m)

## Do the DIS fit
th.model.fix_parameters('ALL')
th.model.release_parameters('NS', 'AL0S', 'AL0G')
f = Fitter.FitterMinuit(DISpoints, thGepard)
f.fit()

thGepard.print_chisq(f.fitpoints)
mGepard.print_parameters_errors()

## Do the DVCS fit
th.model.fix_parameters('ALL')
th.model.release_parameters('Mv', 'rv', 'bv', 'C', 'MC', 
'tMv', 'trv', 'tbv', 'rpi', 'Mpi', 'M02S', 'SECS', 'THIS', 'SECG', 'THIG')
# Next line is cheating, of course, to save cca. 1 hour by
# going immediately to minimum
m.parameters.update(KMM12)
f = Fitter.FitterMinuit(GLOnoBSS2 + BSSwpoints, th)
f.minuit.tol = 80
f.minuit.printMode = 1
f.minuit.maxcalls = 1000
f.fit()

th.print_chisq(f.fitpoints)
m.print_parameters_errors()

