#!/usr/bin/env python

import shelve

import numpy as np

import Model, Approach, Fitter, Data
import utils 
import plots

from constants import Mp, Mp2

from results import *

# [1] Load experimental data and theoretical models

data = utils.loaddata('data/ep2epgamma', approach=Approach.hotfixedBMK)  
data.update(utils.loaddata('data/gammastarp2gammap', approach=Approach.hotfixedBMK))
db = shelve.open('theories.db')

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
#ALTGLO4points = data[25] + data[32][18:]
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


## [3] Create a theory

# Gepard only
mGepard = Model.ComptonGepard(cutq2=0.5)
tGepard = Approach.hotfixedBMK(mGepard)


# DR only
mDRonly = Model.ModelDR()
tDR = Approach.hotfixedBMK(mDRonly)
tDR.name = '(H1/ZEUS)+HERMES+CLAS'
tDR.m.parameters.update(DMepsGLO)

mDRonly1 = Model.ModelDR()
tDR1 = Approach.hotfixedBMK(mDRonly1)
tDR1.name = '(H1/ZEUS)+HERMES+CLAS+HallA'
tDR1.m.parameters.update(DMepsGLO1)


## Hybrid: Gepard+DR (can reuse above Gepard)
#mDRsea = Model.ComptonModelDRsea()
#m = Model.Hybrid(mGepard, mDRsea)
#t = Approach.hotfixedBMK(m)
#t.name = 'DR + Gepard sea'
#g = t.m.g

mDRPPsea = Model.ComptonModelDRPPsea()
m = Model.Hybrid(mGepard, mDRPPsea)
t = Approach.BM10(m)
t.name = 'H1/ZEUS+HERMES+CLAS+HallA'
g = t.m.g

#t.m.parameters.update(hy1THI)
#t.m.parameters.update(KKunp5)
t.m.parameters.update(KKunpTSA1)
#t.m.parameters.update(KKunpTSA1cut16)


# NN
#mNN = Model.ModelNN(hidden_layers=[15], output_layer=['ImH', 'ReH', 'ImE', 'ReE', 'ImHt', 'ReHt', 'ImEt', 'ReEt'])
#mNN = Model.ModelNN(hidden_layers=[9], endpointpower=3.0)
#tNN = Approach.hotfixedBMK(mNN)
#tNN.name = 'NNtest'
#tNN.description = 'x (xB,t)-13-8 nets trained on ALTGLO+BSDw2CD+BSSw for 500 batches'

## [4] Do the fit
#t.m.fix_parameters('ALL')

#f = Fitter.FitterBrain(UNPpoints+TSApoints+BTSApoints, tNN, nnets=10, nbatch=50, verbose=1)
#f.fit()
#f.prune(minprob=0.5)
#tNN.save(db)
#db.close()

## Fitting to both small- and large-x data
#t.m.release_parameters('M02S', 'SECS', 'SECG', 'THIS', 'THIG', 'rv', 'bv', 'Mv', 'C', 'MC', 'trv', 'tbv', 'tMv', 'rpi', 'Mpi')
#f = Fitter.FitterMinuit(H1ZEUSpoints+UNP5points+TSApoints+BTSApoints, t)
#f.minuit.tol = 80
#f.minuit.maxcalls = 100
#f.minuit.printMode = 1
#for par in t.m.parameter_names:
#    if not t.m.parameters['fix_'+par]:
#        print '---- '+par+' ---'
#        t.scan(par, f.fitpoints)


## DR fit
#t.m.release_parameters('rv', 'bv', 'Mv', 'C', 'MC', 'trv', 'tbv', 'tMv', 'rpi', 'Mpi')
#f = Fitter.FitterMinuit(ALTGLOpoints+BSDwpoints+BSSwpoints, t)


def pc(th):
    exps = ['UNP5points', 'H1ZEUS', 'ALTGLO5', 'CLAS', 'CLASDM', 'BSDw', 'BSSw', 'TSA1']
    ptssets = [UNP5points, H1ZEUSpoints, ALTGLO5points, data[25], data[8], BSDwpoints, BSSwpoints, TSA1points]
    for name, pts in zip(exps,ptssets):
        print '%6s: chi/npts = %5.2f/%d' % (name, th.chisq(pts)[0], len(pts))
        cutpts = utils.select(pts, criteria=['Q2>=1.6'])
        print '%6s: chi/npts = %5.2f/%d (cut)' % (name, th.chisq(cutpts)[0], len(cutpts))

