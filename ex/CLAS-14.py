#!/usr/bin/env python
# coding: utf-8

"""Give predictions for CLAS 2014 data according to KMM12 model."""

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


UNCERT = True  # do we want uncertainties
# In[2]:

# Example model
db = shelve.open('/home/kkumer/pyper/theories.db')
th = db['KMM12']


### CLAS

##### Private communication from Silvia Niccolai

# In[4]:

dfBSA = pd.DataFrame([(pt.Q2, pt.xB, pt.tm, pt.phi, pt.val, pt.err, pt) for pt in data[82]], columns=('Q2', 'xB', 'tm', 'phi', 'val', 'err', 'pt'))
dfTSA = pd.DataFrame([(pt.Q2, pt.xB, pt.tm, pt.phi, pt.val, pt.err, pt) for pt in data[83]], columns=('Q2', 'xB', 'tm', 'phi', 'val', 'err', 'pt'))
dfBTSA = pd.DataFrame([(pt.Q2, pt.xB, pt.tm, pt.phi, pt.val, pt.err, pt) for pt in data[84]], columns=('Q2', 'xB', 'tm', 'phi', 'val', 'err', 'pt'))


# In[5]:

all_bins = dfBSA[['Q2', 'xB', 'tm']].drop_duplicates() # vaid also for TSA and BTSA
bins = all_bins[all_bins.tm < all_bins.Q2/4.].reset_index(drop=True)


# In[6]:

binBSA = {}
for k in range(len(bins)):
    binBSA[k] = dfBSA[dfBSA[['Q2', 'xB', 'tm']].values == bins.iloc[k].values].drop_duplicates()
    
binTSA = {}
for k in range(len(bins)):
    binTSA[k] = dfTSA[dfTSA[['Q2', 'xB', 'tm']].values == bins.iloc[k].values].drop_duplicates()
    
binBTSA = {}
for k in range(len(bins)):
    binBTSA[k] = dfBTSA[dfBTSA[['Q2', 'xB', 'tm']].values == bins.iloc[k].values].drop_duplicates()


# In[121]:

phis = np.linspace(-np.pi, np.pi, 25)  # Use 4*n+1 points
def predictbin(cpts, error=False):
    """Predict cpts and return them to Trento frame."""
    pt0 = cpts[0]
    pts = [copy.deepcopy(pt0) for phi in phis]
    for pt, phi in zip (pts, phis): pt.phi = phi
    pred = [th.predict(pt, error=error) for pt in pts]
    Trento_phis = [Approach.BMK.from_conventions(pt.copy()).phi for pt in pts]
    # CLAS points
    cpts = [Approach.BMK.from_conventions(pt.copy()) for pt in cpts]
    return Trento_phis, pred, cpts


# In[154]:

allpred = {}
for nb in range(len(bins)):
    allpred[nb] = tuple([predictbin(bin[nb].pt.values, error=UNCERT) for bin in (binBSA, binTSA, binBTSA)])


# In[155]:

def plotbin(n, error=False):
    """Plot BSA, TSA and BTSA for bin n, after allpred dict is done."""
    obses = ('BSA', 'TSA', 'BTSA')
    fig, axs = plt.subplots(1,3, figsize=[12,6])
    for npn, ax in enumerate(axs):
        Tphis, pred, cpts = allpred[n][npn]
        if error:
            pred = np.array(pred)
            ax.fill_between(Tphis, pred[:,0]-pred[:,1], pred[:,0]+pred[:,1], color='red', alpha=0.9)
        else:
            ax.plot(Tphis, pred, color='red')
        pve = np.array([(pt.phi, pt.val, pt.err) for pt in cpts])
        ax.errorbar(pve[:,0], pve[:,1], pve[:,2], linestyle='None')
        ax.set_ylim(-0.5, 0.5)
        ax.set_xlim(0, 360)
        ax.text(0.9, 0.9, obses[npn], ha='right', transform = ax.transAxes)
        pt0 = cpts[0]
        ax.text(0.03, 0.08, 'bin #{}'.format(n), ha='left', va='bottom', transform=ax.transAxes, fontsize=14)
        ax.text(0.03, 0.03, 'Q2 = {:.3f},  xB = {:.3f},  -t = {:.3f}'.format(pt0.Q2, pt0.xB, pt0.tm), 
                ha='left', va='bottom', transform=ax.transAxes)
    axs[2].set_ylim(0, 1.0)
    fig.savefig('bin-{}.png'.format(n))



# In[156]:

for nb in range(len(bins)):
    plotbin(nb, error=UNCERT)


# In[158]:

# Printing output to file

#fstr = 3*'{:.4f}  ' + '{:5.1f}' + 3*'  {: .3f}'
f = open('CLAS14-KMM12.dat', 'w')
sstr = '# {:^4s}  ' + 2*'{:^6s}  ' + '{:^5s}' + 6*'  {:^6s}' + '\n'
sstr2 = '# {:^4s}  ' + 2*'{:^6s}  ' + '{:^5s}' + 3*'  {:^6s}' + '\n'
fstr = 3*'{:.4f}  ' + '{:5.1f}' + 6*'  {: .3f}' + '\n'
fstr2 = 3*'{:.4f}  ' + '{:5.1f}' + 3*'  {: .3f}' + '\n'
f.write('# 2014-07-24 Predictions of KMM12 model (arXiv:1301.1230 and arXiv:0904.0458)\n')
f.write('#' + 76*'-' + '\n')
if UNCERT:
    f.write(sstr.format('Q2', 'xB', '-t', 'phi', 'BSA', 'delBSA', 'TSA', 'delTSA', 'BTSA', 'delDSA'))
else:
    f.write(sstr2.format('Q2', 'xB', '-t', 'phi', 'BSA', 'TSA', 'BTSA'))
f.write('#' + 76*'-' + '\n')
nb = 0
for Q2, xB, tm in bins[:].values:
    P = np.array(allpred[nb])
    for Tphi, BSA, TSA, BTSA in zip(P[0,0], P[0,1], P[1,1], P[2,1]):
        if UNCERT:
            f.write(fstr.format(Q2, xB, tm, Tphi, BSA[0], BSA[1], TSA[0], TSA[1], BTSA[0], BTSA[1]))
        else:
            f.write(fstr2.format(Q2, xB, tm, Tphi, BSA, TSA, BTSA))
    nb += 1
f.close()


print(' .. Done. ..')


