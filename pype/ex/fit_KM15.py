#!/usr/bin/env python

import shelve, logging
logging.basicConfig(level=logging.INFO)

import Model, Approach, Fitter, Data, utils, plots
from results import *
from utils import listdb
from abbrevs import *

db = shelve.open('theories.db')

# Gepard sea part
mGepard = Model.ComptonGepard(p=0)
Model.ComptonGepard.gepardPool.pop()
# DR part
mDRsea = Model.ComptonModelDRPPsea()
# Hybridization
m = Model.HybridKelly(mGepard, mDRsea)
th = Approach.BM10tw2(m)

## Do the total DVCS fit
th.model.fix_parameters('ALL')
th.model.release_parameters('Mv', 'rv', 'bv', 'C', 'MC', 
'tMv', 'trv', 'tbv', 'rpi', 'Mpi', 'M02S', 'SECS', 'THIS', 'SECG', 'THIG')
pars_Gepard = db['KMM12'].m.Gepard.parameters
pars_DR = db['KMM12'].m.DR.parameters
#del pars_Gepard['EKAPG']
#del pars_Gepard['EKAPS']
m.parameters.update(pars_Gepard)
m.parameters.update(pars_DR)
f = Fitter.FitterMinuit(GLO15b , th)
f.minuit.tol = 80
f.minuit.printMode = 1
f.minuit.maxcalls = 15000
f.fit()

th.print_chisq(f.fitpoints)
m.print_parameters_errors()

