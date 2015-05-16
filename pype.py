#!/usr/bin/env python

import shelve, copy, sys, logging, __builtin__
import numpy as np

logging.basicConfig(level=logging.INFO)
#logging.basicConfig(level=logging.DEBUG)

import Model, Approach, Fitter, Data, utils, plots
from results import *
from utils import listdb
from abbrevs import *

db = shelve.open('theories.db')

th = db['KMM12']
Model.ComptonGepard.gepardPool.pop()
thb = db['KMM12']
Model.ComptonGepard.gepardPool.pop()
thb.m.parameters.update(KMM12b)
thC = db['KM15prelimC']
thC.name = 'KM15pC'
Model.ComptonGepard.gepardPool.pop()

#mGK = Model.GK12()
#thGK = Approach.BM10(mGK)
#thGK.name = 'GK12'

pt0 = DISpoints[0]
pt0.Q2 = 4.
pt0.xB = 1e-4
pt0.xi = pt0.xB/(2.-pt0.xB)

## Gepard sea part
#mGepard = Model.ComptonGepard(p=0, q02=4.0)
#Model.ComptonGepard.gepardPool.pop()
#thGepard = Approach.BM10(mGepard)
## DR part
#mDRsea = Model.ComptonModelDRPPsea()
## Hybridization
#m = Model.HybridDipole(mGepard, mDRsea)
#th = Approach.BM10(m)

pts = GLOnoBSS2 + BSSwpoints

GLOfix = (ALUIpts + BCApts + CLASpts + BSDwpoints + 
            AULpts + ALLpts + AUTIpts + BSSwpoints)

## [4] Do the fit
#th.model.fix_parameters('ALL')
#th.model.release_parameters(
#   'ALPS', 'M03S', 'SECS', 'THIS', 'ALPG', 'M02G', 'SECG', 'THIG')

#datcut = utils.select(H109XL+H109WdepXL+H1ZEUScut, criteria=['Q2 >= 4.0'])
#f = Fitter.FitterMinuit(datcut, th)
#f.minuit.tol = 80
#f.minuit.printMode = 1
#f.minuit.maxcalls = 1000

allpars = ['Mv', 'rv', 'bv', 'C', 'MC', 
'tMv', 'trv', 'tbv', 'rpi', 'Mpi', 'M02S', 'SECS', 'THIS', 'SECG', 'THIG']

#f.fit()
#fl = open('aux.par')
#fl.write(str(th.m.chisq(f.fitpoints)))
#fl.write(str(th.m.parameters))
#fl.write('\n')
#fl.write(str(th.m.covariance))
#fl.close()
#f.fit()


## [5] Some shortcuts ...

def thvsth(th1, th2, pts):
    for pt in pts:
        try:
            p1 = th1.predict(pt)
            p2 = th2.predict(pt)
            print '%5s  %6.1e  %4.1f  =  %+8.4e vs %+8.4e is %6.3f' % (
                    pt.y1name, pt.xB, pt.Q2,
                    p1, p2, abs((p2-p1)/p1))
        except:
            print '%4s %.1g = failed' % (pt.y1name, pt.xB)


def gvsDR(th, pts, attr='ImH'):
    for pt in pts:
        try:
            print '%4s %.1g = %.1g vs %.1f' % (
                    pt.y1name, pt.xB, 
                    getattr(th.m.Gepard, attr)(pt), 
                    getattr(th.m.DR, attr)(pt))
        except:
            print '%4s %.1g = failed' % (pt.y1name, pt.xB)


def pc(th, Q2cut=2., wcut=False):
    #exps = ['UNP5points', 'ALTGLO5', 'CLAS', 'CLASDM', 'BSDw', 'BSSw', 'TSA1', 'BTSA', 'TPpoints']
    #ptssets = [UNP5points, ALTGLO5points, data[25], data[8], BSDwpoints, BSSwpoints, TSA1points, BTSApoints, TPpoints]
    exps = ['H1ZEUS', 'ALUIpts', 'BCApts', 'CLASpts', 'BSDw06', 'BSSw06', 'AULpts', 'ALLpts', 'AUTIpts', 'BSDw_C', 'BSS0w_C', 'BSS1w_C', 'BSDw_HA', 'BSS0w_HA', 'BSS1w_HA' ]
    ptssets = [H1ZEUS, ALUIpts, BCApts, CLASpts, BSDwpoints, BSSwpoints, AULpts, ALLpts, AUTIpts,
            data[101], data[102][:48], data[102][48:], data[117], data[116][:10], data[116][10:] ]
    #exps = ['H1ZEUS DVCS', 'H1-09 XL', "H1-09 W-dep"]
    #ptssets = [H1ZEUS, H109XL, H109WdepXL]
    for name, pts in zip(exps,ptssets):
        print '%10s: chi/npts = %6.2f/%d' % (name, th.chisq(pts)[0], len(pts))
        if wcut:
            cutpts = utils.select(pts, criteria=['Q2>=%f' % Q2cut])
            print '%10s: chi/npts = %6.2f/%d (cut)' % (name, th.chisq(cutpts)[0], len(cutpts))

def _derpt(th, p, pt, f=False, h=0.05):
    """Compute derivative of f w.r.t. model parameter p at point pt.
    
    Simple difference is used (f(p+h/2)-f(p-h/2))/h.
    f is string representing appropriate method of th, or
    observable will be taken as yaxis of pt

    """
    if f:
        fun = th.__getattribute__(f) 
    else:
        fun = th.__getattribute__(pt.yaxis)
    mem = th.m.parameters[p]
    th.m.parameters[p] = mem+h/2.
    up = fun(pt)
    th.m.parameters[p] = mem-h/2.
    down = fun(pt)
    th.m.parameters[p] = mem
    return (up-down)/h


def der(th, pars, pts, f=False,  h=0.05):
    """Compute average derivative at points for each par in pars."""

    for par in pars:
        ders = np.array([_derpt(th, par, pt, f, h) for pt in pts])
        print '%4s  |  %7.4f' % (par, ders.mean())
