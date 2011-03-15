#!/usr/bin/env python

#import pylab
#import matplotlib.pyplot as plt

import shelve, copy

import numpy as np
import scipy.stats

import Model, Approach, Fitter, Data, utils, plots
from constants import Mp, Mp2

from results import *
from math import sqrt

## [1] Load experimental data and theoretical models

data = utils.loaddata('data/ep2epgamma', approach=Approach.hotfixedBMK)  
data.update(utils.loaddata('data/gammastarp2gammap', approach=Approach.hotfixedBMK))
db = shelve.open('aux.db')
#ths = shelve.open('theories.db')
dell = shelve.open('dell.db')
dellB = shelve.open('dellB.db')

## [2] Choose subset of datapoints for fitting

testpoints = data[31][12:14] + data[8][1:3] + data[30][2:4]  # test set
GLOpoints = data[31][12:] + data[8] + data[29]  # DM's GLO set
GLO1points = data[31][12:] + data[8] + data[29] + data[30]  # DM's GLO1 set
#HERMESpoints = data[31][12:] +  data[29]
#BSApoints = data[8] + data[29]
#HAD17 = utils.select(data[33], criteria=['Q2 == 1.5', 't == -0.17'])
#HA17 = utils.select(data[34], criteria=['t == -0.17'])
#HA28 = utils.select(data[34], criteria=['t == -0.28'])
#HA33 = utils.select(data[34], criteria=['t == -0.33'])
DVCSpoints = data[36] + data[37] + data[38] + data[39] + \
  data[40] + data[41] + data[42] + data[43] + data[44] + \
  data[45]
ALTGLOpoints = data[5] + data[25] + data[32][18:]
#ALTGLO1points = data[5] + data[25] + data[32] + HAD17 + HA17
#ALTGLO2points = data[5] + data[25] + data[32][18:] + HAD17[::2] + HA17[::2]
#ALTGLO3points = data[5] + data[25] + data[32][18:] + data[30] + HA17
ALTGLO4points = data[25] + data[32][18:]
ALTGLO5points = data[5] + data[8] + data[32][18:]   # DM's CLAS BSA
#BSDw2Cpoints = utils.select(data[26], criteria=['Q2 == 2.3'])
#BSDw2CDpoints = utils.select(data[50], criteria=['Q2 == 2.3'])
BSDwpoints = utils.select(data[50], criteria=['FTn == -1'])
BSSwpoints = utils.select(data[51], criteria=['FTn>=0', 'FTn <= 1'])
BSDwDMpoints = utils.select(data[55], criteria=['FTn == -1'])
BSSwDMpoints = utils.select(data[56], criteria=['FTn>=0', 'FTn <= 1'])
TSA1points = utils.select(data[52], criteria=['FTn == -1'])
DMpoints = data[5] + data[32][18:] + data[8] + data[30]
DMTSApoints = data[5] + data[32][18:] + TSA1points + data[8] + data[30]
TSApoints = TSA1points + data[54]
BTSApoints = utils.select(data[53], criteria=['FTn==0'])
UNPpoints = ALTGLOpoints + BSSwpoints + BSDwpoints
UNP5points = ALTGLO5points + BSSwpoints + BSDwpoints
H1ZEUSpoints = DVCSpoints + data[48]

mockH = data[1001][::27] + data[1002][::27]
oldmockH = []
np.random.seed(42)  # For reproducibility
for pt in mockH:
    oldmockH.append(copy.deepcopy(pt))
    pt.err = 0.05
    pt.val = pt.val + scipy.stats.norm.rvs(loc=0, scale=0.05)
#mockHHa = mockH + data[1003]  # Adding Hall A

## [3] Create a theory

# Gepard only
mGepard = Model.ComptonGepard(cutq2=0.5)
tGepard = Approach.hotfixedBMK(mGepard)

# DR only
mDRonly = Model.ModelDR()
tDR = Approach.hotfixedBMK(mDRonly)
tDR.name = 'KM09a'
tDR.m.parameters.update(DMepsGLO)

tDRBM = Approach.BM10(mDRonly)
tDRBM.name = 'KM09a+BM10'
tDRBM.m.parameters.update(DMepsGLO)

mDRonly1 = Model.ModelDR()
tDR1 = Approach.hotfixedBMK(mDRonly1)
tDR1.name = 'KM09b'
tDR1.m.parameters.update(DMepsGLO1)

tDR1BM = Approach.BM10(mDRonly1)
tDR1BM.name = 'KM09b+BM10'
tDR1BM.m.parameters.update(DMepsGLO1)

## Hybrid: Gepard+DR (can reuse above Gepard)
#mDRsea = Model.ComptonModelDRsea()
#m = Model.HybridDipole(mGepard, mDRsea)
#th = Approach.hotfixedBMK(m)
#th.name = 'KM10b'
#th.m.parameters.update(KM10b)  
#g = th.m.g

#thBM = Approach.BM10(m)
#thBM.m.name = "KM10b+BM10"


mDRPPsea = Model.ComptonModelDRPPsea()
#m = Model.HybridDipole(mGepard, mDRPPsea)
m = Model.HybridKelly(mGepard, mDRPPsea)
th = Approach.BM10(m)
th.name = 'KM10b'
g = th.m.g
th.m.parameters.update(KM10b)



# NN
#mNN = Model.ModelNN(hidden_layers=[15], output_layer=['ImH', 'ReH', 'ImE', 'ReE', 'ImHt', 'ReHt', 'ImEt', 'ReEt'])
#mNN = Model.ModelNN(hidden_layers=[9], endpointpower=3.0)
mNN = Model.ModelNN(hidden_layers=[11])
#mNN = Model.ModelNN(hidden_layers=[11], output_layer=['ImH', 'ReH', 'ReE','ImHt', 'ReHt', 'ReEt'])
tNN = Approach.hotfixedBMK(mNN)
tNN.name = 'NNtest'
tNN.description = 'x (xB,t)-11-6 nets trained on HERMES+CLAS+HallA mock data for 100 batches'

## [4] Do the fit
#th.m.fix_parameters('ALL')

#f = Fitter.FitterBrain(traindata, tNN, nnets=20, nbatch=200, verbose=1)
#f = Fitter.FitterBrain(traindata, tNN, nnets=8, nbatch=400, verbose=2)
#f.fit()
#f.prune(minprob=0.5)
#tNN.m.parameters['nnet'] = 'ALL'
#tNN.save(db)
#db.close()

## Fitting to both small- and large-x data
#th.m.release_parameters('M02S', 'SECS', 'SECG', 'THIS', 'THIG', 'rv', 'bv', 'Mv', 'C', 'MC', 'trv', 'tbv', 'tMv')#, 'rpi', 'Mpi')
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
th.m.release_parameters('rv', 'bv', 'Mv', 'C', 'MC')
f = Fitter.FitterMinuit(mockH, th)

## [5] Some shortcuts ...

def ld(db):
    utils.listdb(db)

def pc(th):
    exps = ['UNP5points', 'H1ZEUS', 'ALTGLO5', 'CLAS', 'CLASDM', 'BSDw', 'BSSw', 'TSA1']
    ptssets = [UNP5points, H1ZEUSpoints, ALTGLO5points, data[25], data[8], BSDwpoints, BSSwpoints, TSA1points]
    for name, pts in zip(exps,ptssets):
        print '%6s: chi/npts = %5.2f/%d' % (name, th.chisq(pts)[0], len(pts))
        cutpts = utils.select(pts, criteria=['Q2>=1.6'])
        print '%6s: chi/npts = %5.2f/%d (cut)' % (name, th.chisq(cutpts)[0], len(cutpts))

    
# fixed target datapoint FIXME: WRONG!
def ptfix(th, Q=-1, pol=1, Ee=30., xB=0.1, Q2=4., t=-0.2, phi=1., FTn=None):
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

ptc = ptcol(th, Q=1, pol=0, Ee=5, Ep=350, phi=3.14/2., xB=1e-2)

