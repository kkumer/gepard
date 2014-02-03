#!/usr/bin/env python

import shelve, copy, sys, logging, __builtin__

logging.basicConfig(level=logging.INFO)

import Model, Approach, Fitter, Data, utils, plots
from results import *
from utils import listdb
from abbrevs import *

db = shelve.open('/home/kkumer/pype/theories.db')

thKM10 = db['KM10']
Model.ComptonGepard.gepardPool.pop()


## [4] Do the fit
#th.model.fix_parameters('ALL')
#th.model.release_parameters(
#   'ALPS', 'M03S', 'SECS', 'THIS', 'ALPG', 'M02G', 'SECG', 'THIG')

#datcut = utils.select(H109XL+H109WdepXL+H1ZEUScut, criteria=['Q2 >= 4.0'])
#f = Fitter.FitterMinuit(datcut, th)
#f.minuit.tol = 80
#f.minuit.printMode = 1
#f.minuit.maxcalls = 1000


#f.fit()
#fl = open('aux.par')
#fl.write(str(th.m.chisq(f.fitpoints)))
#fl.write(str(th.m.parameters))
#fl.write('\n')
#fl.write(str(th.m.covariance))
#fl.close()
#f.fit()


## [5] Some shortcuts ...



def pc(th, Q2cut=4.):
    #exps = ['UNP5points', 'ALTGLO5', 'CLAS', 'CLASDM', 'BSDw', 'BSSw', 'TSA1', 'BTSA', 'TPpoints']
    #ptssets = [UNP5points, ALTGLO5points, data[25], data[8], BSDwpoints, BSSwpoints, TSA1points, BTSApoints, TPpoints]
    #exps = ['H1ZEUS', 'ALUIpts', 'BCApts', 'CLASpts', 'BSDwpoints', 'BSSwpoints', 'AULpts', 'ALLpts', 'AUTIpts' ]
    #ptssets = [H1ZEUS, ALUIpts, BCApts, CLASpts, BSDwpoints, BSSwpoints, AULpts, ALLpts, AUTIpts ]
    exps = ['H1ZEUS DVCS', 'H1-09 XL', "H1-09 W-dep"]
    ptssets = [H1ZEUS, H109XL, H109WdepXL]
    for name, pts in zip(exps,ptssets):
        print '%10s: chi/npts = %6.2f/%d' % (name, th.chisq(pts)[0], len(pts))
        cutpts = utils.select(pts, criteria=['Q2>=%f' % Q2cut])
        print '%10s: chi/npts = %6.2f/%d (cut)' % (name, th.chisq(cutpts)[0], len(cutpts))

