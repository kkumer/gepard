#!/usr/bin/env python

#import pylab
#import matplotlib.pyplot as plt

import shelve, copy, sys, logging, builtins

import numpy as np
import scipy.stats

logging.basicConfig(level=logging.INFO)
#logging.basicConfig(level=logging.DEBUG)

import Model, Approach, Fitter, Data, utils, plots
from constants import Mp, Mp2

from results import *
from math import sqrt

from dispersion import *
from quadrature import rthtquadrature
from utils import listdb

## [1] Load experimental data and theoretical models

from abbrevs import *

H1ZEUScut = utils.select(H1ZEUS, criteria=['Q2 >= 2.0'])

# transform total X data into XL=X/(eps+1/R) data
def Rfun(Q2):
    p = 0.447; delp = 0.064
    a = 2.1; dela = 1.
    x = Q2/0.776**2
    R = x/(1+a*x)**p
    delR = sqrt( p**2 * x**2 * R**2 * dela**2 / (1.+a*x)**2 +
	   log(1.+a*x)**2 * R**2 * delp**2) 
    return (R, delR)

H109WdepXL = []
for pt in utils.select(data[79], criteria=['Q2 >= 4.0']):
    ptxl = copy.deepcopy(pt)
    y = (pt.W**2 + pt.Q2 - Mp2)/(pt.s - Mp2)
    eps = (1-y)/(1-y+y**2/2)
    R, delR = Rfun(pt.Q2)
    ptxl.val = pt.val/(eps + 2./R)
    errsig = pt.err / (eps+1./R)
    errR = pt.val * delR / (1.+eps*R)**2
    ptxl.err = sqrt( errsig**2 + errR**2)
    #print ptxl.err/ptxl.val*100, errR/errsig
    ptxl.in1polarization = 1
    ptxl.in1polarizationvector = 'L'
    ptxl.y1namelong = 'differential cross section XL'
    H109WdepXL.append(ptxl)

builtins.H109WdepXL = H109WdepXL

## [3] Create a theory

db = shelve.open('/home/kkumer/pype/theories.db')
thdvcs = db['KMM12nlo2']
Model.ComptonGepard.gepardPool.pop()

#thLO = db['dvmp']
#thLO.name = 'LO'
#Model.ComptonGepard.gepardPool.pop()
#thNLO = db['dvmpnlo']
#Model.ComptonGepard.gepardPool.pop()
#thNLO.name = 'NLO'
#th = thAFKM12
#Model.ComptonGepard.gepardPool.pop()
#thKM10 = db['KM10']
#Model.ComptonGepard.gepardPool.pop()
#theories = [thAFKM12, thKM10]

m = Model.ComptonGepard(p=1, q02=2.)
Model.ComptonGepard.gepardPool.pop()
th = Approach.BMK(m)
#th.m.parameters.update(KMM12) #  LO
th.m.parameters.update(thdvcs.m.parameters)
th.name = 'NLO dvem prelim1'

ptc = Data.DummyPoint()
ptc.process = 'gammastarp2rho0p'
ptc.W = 75.
ptc.Q2 = 4.
ptc.t = -0.0
ptc.s = 320*320
utils.fill_kinematics(ptc)

astrong = 0.3403676
prefac = (4./3.) * 0.209 * astrong / 3 / sqrt(ptc.Q2)


## [4] Do the fit
th.model.fix_parameters('ALL')
#th.model.release_parameters(
#   'ALPS', 'M03S', 'SECS', 'THIS', 'ALPG', 'M02G', 'SECG', 'THIG')
#th.model.release_parameters(
#   'EAL0S', 'EALPS', 'EM02S', 'ESECS', 'ETHIS', 'KAPS',
#   'EAL1G', 'EM02G',  'ESECG')
#th.model.release_parameters(
#  'KAPS',  'M02S', 'M02G', 'SECS', 'SECG', 'EAL0S', 'EM02S', 'ESECS', 'EAL0G')
th.model.release_parameters('M02S', 'M02G', 'SECS', 'SECG', 'THIS', 'THIG')
#th.model.release_parameters(
#   'rv', 'Mv', 'bv', 'C', 'MC', 'trv', 'tbv')
#th.model.release_parameters('M02S', 'SECS', 'SECG', 'THIS', 'THIG', 
#   'rv', 'bv', 'Mv', 'C', 'MC', 'trv', 'tbv', 'tMv', 'rpi', 'Mpi')
#f = Fitter.FitterMinuit(GLOnoBSS2+BSSwpoints, th)

datcut = utils.select(H109XL+H109WdepXL+H1ZEUScut, criteria=['Q2 >= 2.0'])
f = Fitter.FitterMinuit(datcut, th)
f.minuit.tol = 80
f.minuit.printMode = 1
f.minuit.maxcalls = 3000
f.fit()


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
        print('%10s: chi/npts = %6.2f/%d' % (name, th.chisq(pts)[0], len(pts)))
        cutpts = utils.select(pts, criteria=['Q2>=%f' % Q2cut])
        print('%10s: chi/npts = %6.2f/%d (cut)' % (name, th.chisq(cutpts)[0], len(cutpts)))


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
        print('%4s  |  %5.2f' % (par, ders.mean()))

#def rth(m, pt, tht):
#    """Calculate <r(tht)> for model m."""
#    norm = rthtquadrature(lambda r: m.gpdHQbpolpol(pt, r, tht), 0.0, 2.5)
#    aux = rthtquadrature(lambda r: r*m.gpdHQbpolpol(pt, r, tht), 0.0, 2.5)
#    return aux/norm

#def bsq(m, pt):
#    """Calculate <b^2> for model m."""
#    norm = rthtquadrature(lambda b: b*m.gpdHQbIntg(pt, b), 0.0, 2.5)
#    aux = rthtquadrature(lambda b: b**3*m.gpdHQbIntg(pt, b), 0.0, 2.5)
#    return aux/norm
