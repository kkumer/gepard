#!/usr/bin/env python

#import pylab
#import matplotlib.pyplot as plt

import shelve, copy, sys

import numpy as np
import scipy.stats

import Model, Approach, Fitter, Data, utils, plots
from constants import Mp, Mp2

from results import *
from math import sqrt

## [1] Load experimental data and theoretical models

data = utils.loaddata('/home/kkumer/pype/data/ep2epgamma', approach=Approach.hotfixedBMK)  
data.update(utils.loaddata('/home/kkumer/pype/data/gammastarp2gammap', approach=Approach.hotfixedBMK))
db = shelve.open('/home/kkumer/pype/theories.db')
dell = shelve.open('/home/kkumer/pype/dellB.db')

## [2] Choose subset of datapoints for fitting

#testpoints = data[31][12:14] + data[8][1:3] + data[30][2:4]  # test set
GLOpoints = data[31][12:] + data[8] + data[29]  # DM's GLO set
GLO1points = data[31][12:] + data[8] + data[29] + data[30]  # DM's GLO1 set
##HERMESpoints = data[31][12:] +  data[29]
##BSApoints = data[8] + data[29]
##HAD17 = utils.select(data[33], criteria=['Q2 == 1.5', 't == -0.17'])
##HA17 = utils.select(data[34], criteria=['t == -0.17'])
##HA28 = utils.select(data[34], criteria=['t == -0.28'])
##HA33 = utils.select(data[34], criteria=['t == -0.33'])
DVCSpoints = data[36] + data[37] + data[38] + data[39] + \
  data[40] + data[41] + data[42] + data[43] + data[44] + \
  data[45]
ALTGLOpoints = data[5] + data[25] + data[32][18:]  # KK's CLAS BSA
##ALTGLO1points = data[5] + data[25] + data[32] + HAD17 + HA17
##ALTGLO2points = data[5] + data[25] + data[32][18:] + HAD17[::2] + HA17[::2]
##ALTGLO3points = data[5] + data[25] + data[32][18:] + data[30] + HA17
#ALTGLO4points = data[25] + data[32][18:]
ALTGLO5points = data[5] + data[8] + data[32][18:]   # DM's CLAS BSA
Hpoints = data[5] + data[32][18:]
##BSDw2Cpoints = utils.select(data[26], criteria=['Q2 == 2.3'])
##BSDw2CDpoints = utils.select(data[50], criteria=['Q2 == 2.3'])
BSDwpoints = utils.select(data[50], criteria=['FTn == -1'])
BSSwpoints = utils.select(data[51], criteria=['FTn>=0', 'FTn <= 1'])
#BSDwDMpoints = utils.select(data[55], criteria=['FTn == -1'])
#BSSwDMpoints = utils.select(data[56], criteria=['FTn>=0', 'FTn <= 1'])
TSA1points = utils.select(data[52], criteria=['FTn == -1'])
#DMpoints = data[5] + data[32][18:] + data[8] + data[30]
#DMTSApoints = data[5] + data[32][18:] + TSA1points + data[8] + data[30]
#TSApoints = TSA1points + data[54]
BTSApoints = utils.select(data[53], criteria=['FTn==0'])
UNPpoints = ALTGLOpoints + BSSwpoints + BSDwpoints
UNP5points = ALTGLO5points + BSSwpoints + BSDwpoints
H1ZEUSpoints = DVCSpoints + data[48]


## [3] Create a theory

# Gepard only
mGepard = Model.ComptonGepard(cutq2=0.5)
#tGepard = Approach.hotfixedBMK(mGepard)

# DR only
#mDRonly = Model.ModelDR()
#tDR = Approach.hotfixedBMK(mDRonly)
#tDR.name = 'KM09a'
#tDR.m.parameters.update(DMepsGLO)

#tDRBM = Approach.BM10(mDRonly)
#tDRBM.name = 'KM09a+BM10'
#tDRBM.m.parameters.update(DMepsGLO)

#mDRonly1 = Model.ModelDR()
#tDR1 = Approach.hotfixedBMK(mDRonly1)
#tDR1.name = 'KM09b'
#tDR1.m.parameters.update(DMepsGLO1)

#tDR1BM = Approach.BM10(mDRonly1)
#tDR1BM.name = 'KM09b+BM10'
#tDR1BM.m.parameters.update(DMepsGLO1)

## Hybrid: Gepard+DR (can reuse above Gepard)
#mDRsea = Model.ComptonModelDRsea()
#m = Model.HybridDipole(mGepard, mDRsea)
#th = Approach.hotfixedBMK(m)
#th.name = 'KM10a'
#th.m.parameters.update(KM10a)  
#g = th.m.g

#thBM = Approach.BM10(m)
#thBM.m.name = "KM10b+BM10"


mDRPPsea = Model.ComptonModelDRPPsea()
#m = Model.HybridDipole(mGepard, mDRPPsea)
m = Model.HybridKelly(mGepard, mDRPPsea)
th = Approach.BM10(m)
th.name = 'new'
g = th.m.g
th.m.parameters.update(zeroHt)
#th.m.covariance = KM10cov

#thHF = Approach.hotfixedBMK(m)
#thHF.name = 'KM10 - hotfixed'
#thHF.m.parameters.update(KM10)


# NN
#mNN = Model.ModelNN(hidden_layers=[15], output_layer=['ImH', 'ReH', 'ImE', 'ReE', 'ImHt', 'ReHt', 'ImEt', 'ReEt'])
#mNN = Model.ModelNN(hidden_layers=[9], endpointpower=3.0)
#mNN = Model.ModelNN(hidden_layers=[11])
#mNN = Model.ModelNN(hidden_layers=[11], output_layer=['ImH', 'ReH', 'ReE','ImHt', 'ReHt', 'ReEt'])
#tNN = Approach.hotfixedBMK(mNN)
#tNN.name = 'NNtest'
#tNN.description = 'x (xB,t)-11-2 nets trained on HERMES-like mock data for 200 batches'

#tDR2 = db['DR2-H']
#tDR3 = db['DR3-H']
#tNN = db['NN-H']
#tNNf = dell['NN-H-final']

#tNNf.name = 'Neural nets'
#tDR3.name = 'Model fit'

## [4] Do the fit
th.m.fix_parameters('ALL')

#f = Fitter.FitterBrain(Hpoints, tNNf, nnets=50, nbatch=400, verbose=1)
#f.fit()
#f.prune(minprob=0.5)
#tNN.m.parameters['nnet'] = 'ALL'
#tNN.save(db)
#db.close()

## Fitting to both small- and large-x data
#th.m.release_parameters('M02S', 'SECS', 'SECG', 'THIS', 'THIG', 'rv', 'bv', 'Mv', 'C', 'MC', 'trv', 'tbv', 'tMv')#, 'rpi', 'Mpi')
th.m.release_parameters('rv', 'bv', 'Mv', 'C', 'MC', 'trv', 'tbv', 'tMv', 'rpi', 'Mpi')
#f = Fitter.FitterMinuit(H1ZEUSpoints+UNP5points+TSApoints+BTSApoints, th)
#f = Fitter.FitterMinuit(DVCSpoints+GLOpoints, th)
#f.minuit.tol = 80
#f.minuit.maxcalls = 100
#f.minuit.printMode = 1
#for par in th.m.parameter_names:
#    if not th.m.parameters['fix_'+par]:
#        print '---- '+par+' ---'
#        th.scan(par, f.fitpoints)


## DR fit
#tDR.m.release_parameters('rv', 'Mv')
#tDR.m.release_parameters('rv', 'Mv')
#f = Fitter.FitterMinuit(mockH, tDR)

## [5] Some shortcuts ...

def ld(db):
    utils.listdb(db)

def pc(th):
    exps = ['UNP5points', 'ALTGLO5', 'CLAS', 'CLASDM', 'BSDw', 'BSSw', 'TSA1', 'BTSA']
    ptssets = [UNP5points, ALTGLO5points, data[25], data[8], BSDwpoints, BSSwpoints, TSA1points, BTSApoints]
    for name, pts in zip(exps,ptssets):
        print '%10s: chi/npts = %6.2f/%d' % (name, th.chisq(pts)[0], len(pts))
        #cutpts = utils.select(pts, criteria=['Q2>=1.6'])
        #print '%10s: chi/npts = %6.2f/%d (cut)' % (name, th.chisq(cutpts)[0], len(cutpts))

    
# fixed target datapoint FIXME: WRONG!
def ptfix(th, Q=1, pol=-1, Ee=160., xB=0.1, Q2=2.2, t=-0.1, phi=3.5, FTn=None):
    ptf = Data.DummyPoint()
    ptf.in1energy = Ee
    ptf.s = 2 * Mp * ptf.in1energy + Mp2
    ptf.in1charge = Q
    ptf.in1polarization = pol
    ptf.xB = xB
    ptf.Q2 = Q2
    ptf.t = t
    ptf.xi = ptf.xB/(2.-ptf.xB)
    ptf.phi = phi
    ptf.frame = 'Trento'
    ptf.units = {'phi': 'radian'}
    utils.fill_kinematics(ptf)
    th.to_conventions(ptf)
    th.prepare(ptf)
    return ptf


# collider datapoint
def ptcol(th, Q=-1, pol=1, Ee=5, Ep=250, xB=0.001, Q2=4., t=-0.2, phi=1., FTn=None):
    ptc = Data.DummyPoint()
    ptc.in1energy = Ee
    ptc.in2energy = Ep
    ptc.s = 2 * ptc.in1energy * (ptc.in2energy + sqrt(
                ptc.in2energy**2 - Mp2)) + Mp2
    ptc.in1charge =  Q
    ptc.in1polarization = pol
    ptc.in2polarization = 1 # relevant only for XLP and TSA
    ptc.xB = xB
    ptc.Q2 = Q2
    ptc.t = t
    ptc.xi = ptc.xB/(2.-ptc.xB)
    if phi:
        ptc.phi = phi
        ptc.units = {'phi': 'radian'}
    elif FTn:
        ptc.FTn = FTn
    ptc.frame = 'Trento'
    utils.fill_kinematics(ptc)
    th.to_conventions(ptc)
    th.prepare(ptc)
    return ptc

#ptc = ptcol(th, Q=1, pol=0, Ee=5, Ep=350, phi=3.14/2., xB=1e-2)

def ccals(th, pt):
    cals = ['DVCSunp', 'INTunp', 'INTunpV', 'INTunpA']
    for c in cals:
        print '%8s =  %10.5f + %10.5f * I ' % (c, getattr(th, 'CCAL'+c)(pt), 
                getattr(th, 'CCAL'+c)(pt, im=1))


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
        print '%4s  |  %5.2f' % (par, ders.mean())

