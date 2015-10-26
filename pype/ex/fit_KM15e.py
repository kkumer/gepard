#!/usr/bin/env python2

import shelve, logging
logging.basicConfig(level=logging.INFO)

import Model, Approach, Fitter, Data, utils, plots
from results import *
from utils import listdb
from abbrevs import *

dba = shelve.open('my.db')

th = dba['KM15prelimD']
th.name = 'KM15prelimE'

f = Fitter.FitterMinuit(GLO15b , th)
f.minuit.tol = 80
f.minuit.printMode = 1
f.minuit.maxcalls = 3000
f.fit()

th.print_chisq(f.fitpoints)
m.print_parameters_errors()

