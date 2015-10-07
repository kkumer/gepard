#!/usr/bin/env python

import shelve, logging
logging.basicConfig(level=logging.INFO)

import Model, Approach, Fitter, Data, utils, plots
from results import *
from utils import listdb
from abbrevs import *

# Gepard sea part
mGepard = Model.ComptonGepard(p=0)
Model.ComptonGepard.gepardPool.pop()
# DR part
mDRsea = Model.ComptonModelDRPPsea()
# Hybridization
m = Model.HybridDipole(mGepard, mDRsea)
th = Approach.BM10(m)

## Do the total DVCS fit
th.model.fix_parameters('ALL')
th.model.release_parameters('Mv', 'rv', 'bv', 'C', 'MC', 
'tMv', 'trv', 'tbv', 'rpi', 'Mpi', 'M02S', 'SECS', 'THIS', 'SECG', 'THIG')
# Next line helps shorten time
m.parameters.update(KMM12)
f = Fitter.FitterMinuit(GLO15 + data[116] + data[117] + data[101] + data[102], th)
f.minuit.tol = 80
f.minuit.printMode = 1
f.minuit.maxcalls = 8000
f.fit()

th.print_chisq(f.fitpoints)
m.print_parameters_errors()

