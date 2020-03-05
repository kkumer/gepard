#!/usr/bin/env python

import shelve, copy, sys, builtins
import numpy as np

import logging, logzero
_lg = logzero.logger
logzero.loglevel(logging.WARNING)
logzero.logfile("/tmp/{}.log".format(__name__), loglevel=logging.INFO, maxBytes=1000000, backupCount=5)

import Model, Approach, Fitter, Data, utils, plots
from results import *
from utils import listdb
from abbrevs import *


# Some combination of datasets to be fitted to
pts = GLO15b


## --- New fit ---

## --- Creating new model: ---
## Gepard sea part
mGepard = Model.ComptonGepard(p=0, q02=4.0)
Model.ComptonGepard.gepardPool.pop()
thGepard = Approach.BM10(mGepard)
# DR part
mDRsea = Model.ComptonModelDRPPsea()
# Hybridization
m = Model.HybridDipole(mGepard, mDRsea)
_lg.info('New model created: {}'.format(m.__class__))
th = Approach.BM10(m)
th.model.fix_parameters('ALL')
th.model.release_parameters('M02S', 'SECS', 'ALPG', 'M02G', 'SECG')

## --- Fitting: ---
pts = utils.select(H1ZEUS, criteria=['Q2 >= 4.0'])
f = Fitter.FitterMinuit(pts, th)
f.minuit.tol = 80
f.minuit.print_level = 1

_lg.info('Start fit to {} data points'.format(len(f.fitpoints)))
f.fit()
