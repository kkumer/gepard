#!/usr/bin/env python

import os, shelve, copy, sys, builtins
import numpy as np

import logging, logzero
_lg = logzero.logger
logzero.loglevel(logging.INFO)
basename = os.path.splitext(os.path.basename(__file__))[0]
logfilename = "/home/kkumer/tmp/{}.log".format(basename)
logzero.logfile(logfilename,
        loglevel=logging.INFO, maxBytes=1000000, backupCount=5)

import Model, Approach, Fitter, Data, utils, plots, shelve
from results import *
from utils import listdb
from abbrevs import *


# Some combination of datasets to be fitted to
pts = GLO15b

db = shelve.open('theories.db')

th = db['KM15']

