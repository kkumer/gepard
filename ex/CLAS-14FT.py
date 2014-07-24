#!/usr/bin/env python
# coding: utf-8

"""Give predictions for harmonics of CLAS 2014 data according to KMM12 model."""

# Results sent to Silvia Niccolai produced with SVN version 283

## Introduction and Initialization

# In[1]:

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import shelve, copy
import plots, utils, Approach
from abbrevs import *


UNCERT = False  # do we want uncertainties
# In[2]:

# Example model
db = shelve.open('/home/kkumer/pyper/theories.db')
th = db['KMM12']


### CLAS


dfBSAFT = pd.DataFrame([(pt.Q2, pt.xB, pt.tm, pt.FTn, pt.val, pt.err, pt) for pt in data[85]], columns=('Q2', 'xB', 'tm', 'FTn', 'val', 'err', 'pt'))
dfTSAFT = pd.DataFrame([(pt.Q2, pt.xB, pt.tm, pt.FTn, pt.val, pt.err, pt) for pt in data[86]], columns=('Q2', 'xB', 'tm', 'FTn', 'val', 'err', 'pt'))
dfBTSAFT = pd.DataFrame([(pt.Q2, pt.xB, pt.tm, pt.FTn, pt.val, pt.err, pt) for pt in data[87]], columns=('Q2', 'xB', 'tm', 'FTn', 'val', 'err', 'pt'))


all_binsFT = dfBSAFT[['Q2', 'xB']].drop_duplicates() # vaid also for TSA and BTSA
binsFT = all_binsFT.reset_index(drop=True)


binBSAFT = {}
for k in range(len(binsFT)):
    binBSAFT[k] = dfBSAFT[dfBSAFT[['Q2', 'xB']].values == binsFT.iloc[k].values].drop_duplicates()
    
binTSAFT = {}
for k in range(len(binsFT)):
    binTSAFT[k] = dfTSAFT[dfTSAFT[['Q2', 'xB']].values == binsFT.iloc[k].values].drop_duplicates()
    
binBTSAFT = {}
for k in range(len(binsFT)):
    binBTSAFT[k] = dfBTSAFT[dfBTSAFT[['Q2', 'xB']].values == binsFT.iloc[k].values].drop_duplicates()


binBTSAFT2 = {}
for k in range(len(binsFT)):
    half = len(binBTSAFT[k])/2
    binBTSAFT2[k] = binBTSAFT[k][half:]  # take just FTn=1


def predictbinFT(cpts, error=False, ntms=25):
    """Predict t-dep of FT bin"""
    pt0 = cpts[0]
    tmin = th.tmin(pt0.Q2, pt0.xB, pt0.eps2)
    # start plotting not from tmin but slightly right from it
    tms = np.linspace(-tmin*1.3, pt0.Q2/4., ntms)
    pred = []
    for tm in tms:
        pt = pt0.copy()
        pt.tm = tm
        pt.t = -tm
        th.prepare(pt)
        pred.append(th.predict(pt, error=error))
    return tms, pred, cpts



allpredFT = {}
for nb in range(len(binsFT)-3):
    allpredFT[nb] = tuple([predictbinFT(bin[nb].pt.values, error=False) for bin in (binBSAFT, binTSAFT, binBTSAFT, binBTSAFT2)])



def plotbinFT(n, error=False):
    """Plot BSA, TSA and BTSA for bin n, after allpredFT dict is done."""
    obses = ('$A_{LU}^{\\sin\\phi}$', '$A_{UL}^{\\sin\\phi}$', '$A_{LL}^{\\cos0\\phi}$', '$A_{LL}^{\\cos\\phi}$')
    FTns = [-1, -1, 0, 1]
    fig, axs = plt.subplots(1,4, figsize=[16,6])
    for npn, ax in enumerate(axs):
        tms, pred, cpts = allpredFT[n][npn]
        if error:
            pred = np.array(pred)
            ax.fill_between(tms, pred[:,0]-pred[:,1], pred[:,0]+pred[:,1], color='red', alpha=0.9)
        else:
            ax.plot(tms, pred, color='red')
        pve = np.array([(pt.tm, pt.val, pt.err) for pt in cpts if pt.FTn == FTns[npn]])
        ax.errorbar(pve[:,0], pve[:,1], pve[:,2], linestyle='None')
        ax.set_ylim(0., 0.4)
        #ax.set_xlim(0, 360)
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))  # tickmarks
        ax.text(0.9, 0.9, obses[npn], ha='right', transform = ax.transAxes, fontsize=20)
        pt0 = cpts[0]
        ax.text(0.03, 0.08, 'binFT #{}'.format(n), ha='left', va='bottom', transform=ax.transAxes, fontsize=14)
        ax.text(0.03, 0.03, 'Q2 = {:.3f},  xB = {:.3f}'.format(pt0.Q2, pt0.xB), 
                ha='left', va='bottom', transform=ax.transAxes)
    axs[2].set_ylim(0, 1.0)
    axs[3].set_ylim(-0.4, 0.4)
    fig.savefig('fig1')



f = open('/home/kkumer/CLAS14FT-KMM12.dat', 'w')
sstr = '# {:^4s}  ' + 2*'{:^6s}  '  + 8*'  {:^6s}' + '\n'
sstr2 = '# {:^4s}  ' + 2*'{:^6s}  ' + 4*'  {:^6s}' + '\n'
fstr = 3*'{:.4f}  ' + 8*'  {: .3f}' + '\n'
fstr2 = 3*'{:.4f}  ' + 4*'  {: .3f}' + '\n'
f.write('# 2014-07-24 Predictions of KMM12 model (arXiv:1301.1230 and arXiv:0904.0458)\n')
f.write('#' + 76*'-' + '\n')
if UNCERT:
    f.write(sstr.format('Q2', 'xB', '-t', 'BSA', 'delBSA', 'TSA', 'delTSA', 'BTSA0', 'del0', 'BTSA1', 'del1'))
else:
    f.write(sstr2.format('Q2', 'xB', '-t', 'BSA', 'TSA',  'BTSA0', 'BTSA1'))
f.write('#' + 76*'-' + '\n')
nb = 0
for Q2, xB in binsFT[:1].values:
    P = np.array(allpredFT[nb])
    for tm, BSA, TSA, BTSA, BTSA1 in zip(P[0,0], P[0,1], P[1,1], P[2,1], P[3,1]):
        if UNCERT:
            f.write(fstr.format(Q2, xB, tm, \
                    BSA[0], BSA[1], TSA[0], TSA[1], BTSA[0], BTSA[1], BTSA1[0], BTSA1[1]))
        else:
            f.write(fstr2.format(Q2, xB, tm, BSA, TSA, BTSA, BTSA1))
    nb += 1
f.close()


for nb in range(len(binsFT)-3):
    plotbinFT(nb, error=UNCERT)


print(' .. Done. ..')


