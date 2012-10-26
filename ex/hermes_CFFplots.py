#!/usr/bin/env python

""" Plotting CFFs resulting from local fits."""

import shelve, copy, sys, logging

import numpy as np
import scipy.stats

logging.basicConfig(level=logging.WARNING)

import Model, Approach, Fitter, Data, utils, plots
from constants import Mp, Mp2

from results import *
from math import sqrt

from dispersion import *
from quadrature import rthtquadrature


## [1] Load experimental data and theoretical models

data = utils.loaddata('/home/kkumer/pype/data/ep2epgamma', approach=Approach.BMK)  

## [2] Choose subset of datapoints for fitting

# Updated data by Morgan and DM
L4_ALUI = utils.select(data[71], criteria=['FTn == -1'])
L4_AC_0 = utils.select(data[70], criteria=['FTn == 0'])
L4_AC_1 = utils.select(data[70], criteria=['FTn == 1'])
# polarized target data
L4_AUL = utils.select(data[52], criteria=['FTn == -1'])
L4_ALL_0 = utils.select(data[53], criteria=['FTn==0'])
L4_ALL_1 = utils.select(data[53], criteria=['FTn==1'])
L4_AUTI_1 = utils.select(data[66], criteria=['FTn==1'])
L4_AUTI_0 = utils.select(data[66], criteria=['FTn==0'])
L4_AUTI_m1 = utils.select(data[66], criteria=['FTn==-1'])
L4_AUTDVCS = data[65]

bins = zip(L4_ALUI, L4_AC_0, L4_AC_1, L4_AUL, L4_ALL_0, 
        L4_ALL_1, L4_AUTI_1, L4_AUTI_0, L4_AUTI_m1, L4_AUTDVCS)
#bins = zip(L4_ALUI, L4_AC_0, L4_AC_1, L4_AUL, L4_ALL_0, 
#        L4_ALL_1, L4_AUTI_1, L4_AUTI_0, L4_AUTI_m1)



## [4] Do the fit bin by bin

#fitCFFs = ['ImH', 'ReH', 'ImHt']

try:
    if sys.argv[1] == 'one':
        fitCFFs = ['ImH', 'ReE']
        fl = open('fits/ReE/ReE.res', 'w')
        fls = [open('fits/ReE/'+f+'.dat', 'w') for f in fitCFFs]
    elif sys.argv[1] == 'two':
        fitCFFs = ['ImH', 'ReH', 'ImHt']
        fl = open('fits/ReH-ImHt/ReH-ImHt.res', 'w')
        fls = [open('fits/ReH-ImHt/'+f+'.dat', 'w') for f in fitCFFs]
    elif sys.argv[1] == 'three':
        fitCFFs = ['ImH', 'ReH', 'ReHt']
        fl = open('fits/ReH-ReHt/ReH-ReHt.res', 'w')
        fls = [open('fits/ReH-ReHt/'+f+'.dat', 'w') for f in fitCFFs]
    else:
        sys.argv[123]  # triggering exception
except IndexError:
    sys.stderr.write('\n Usage: hermes_CFFplots  one | two | three')
    sys.stderr.write('\n one -- Model with ImH and ReE')
    sys.stderr.write('\n two -- Model with ImH, ReH and ImHt')
    sys.stderr.write('\n three -- Model with ImH, ReH and ReHt\n')
    sys.exit(1)

totchi = 0
sys.stdout.write('\n Doing bin ..')
for nbin in range(len(bins)):
    sys.stdout.write('.%i.' % (nbin+1,))
    sys.stdout.flush()
    # Model creation
    m = Model.ModelLocal()
    m.parameters['pImH'] = 10. 
    th = Approach.BM10(m)
    th.model.fix_parameters('ALL')
    for cff in fitCFFs:
        th.model.release_parameters('p%s' % cff)
    f = Fitter.FitterMinuit(bins[nbin], th)
    f.minuit.printMode = 0
    fl.write('\n')
    fl.write('-- Bin %2s: ' % (nbin+1,))
    f.fit()
    pt0 = bins[nbin][3]
    fl.write('t = %6.3f,  xB = %6.3f,  Q2 = %6.3f\n' % (pt0.t, pt0.xB, pt0.Q2))
    fl.write('\n')
    for (p,cf) in zip(th.m.free_parameters(), fls):
        val = th.m.parameters[p]
        err = sqrt(th.m.covariance[p,p])
        cf.write('%i  %f  %f\n' % (nbin+1, val, err))
    #plots.binplot(bands=th, points=f.fitpoints, justbars=True, title=th.name, path='.')
    fl.write('P(chi-square, d.o.f) = P(%1.2f, %2d) = %5.4f\n' % th.chisq(f.fitpoints))
    totchi += th.chisq(f.fitpoints)[0]
    fl.write('\n')
    ndof = len(f.fitpoints) - utils.npars(th.m)
    for p in th.m.free_parameters():
        val = th.m.parameters[p]
        err = sqrt(th.m.covariance[p,p])
        pval = 2*(1.-scipy.stats.t.cdf(abs(val/err), ndof))
        fl.write('%5s = %8.3f +- %5.3f  (p = %.3g)\n' % (p, val, err, pval))

print '\n TOTAL chisq = %f' % (totchi,)
fl.write('====  TOTAL chisq = %f' % (totchi,))
fl.close()
for f in fls:
    f.close()


