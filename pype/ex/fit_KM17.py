#!/usr/bin/env python

import shelve, logging
logging.basicConfig(level=logging.INFO)

import Model, Approach, Fitter, Data, utils, plots
from results import *
from utils import listdb
from abbrevs import *

db = shelve.open('theories.db')
#dba = shelve.open('my.db')

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
#pars_Gepard = dba['HA17B'].m.Gepard.parameters
#pars_DR = dba['HA17B'].m.DR.parameters
pars_Gepard = db['KM15'].m.Gepard.parameters
pars_DR = db['KM15'].m.DR.parameters
m.parameters.update(pars_Gepard)
m.parameters.update(pars_DR)
highQD = utils.select(data[117], ['Q2 > 1.8'])
lowQD = utils.select(data[117], ['Q2 < 1.8'])
#HApts =  data[116]+data[117] +  highQ + highQD
#1 pts = data[135][::2] + data[136][::2] + H1ZEUS
#2 starts from 233 ends at 66.6/41:
#pts = data[135][::2] + data[136][::2] + CLASKKpts + ALUIpts
#3  starts from 204 goes to 31.:
#pts = data[135][::2] + data[136][::2] + BCApts
#4 start from 217 ends at 37.1/52
#pts = data[135][::2] + data[136][::2] + data[101][::4] + data[102][::4]
#5 start from 221
pts = data[135][::2] + data[136][::2] + AULpts + ALLpts + AUTIpts   
f = Fitter.FitterMinuit(pts, th)
f.minuit.tol = 80
f.minuit.printMode = 1
f.minuit.maxcalls = 6000
f.fit()

th.print_chisq(f.fitpoints)
m.print_parameters_errors()

