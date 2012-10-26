#!/usr/bin/env python

"""Doing local fits to HERMES data."""

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
data.update(utils.loaddata('/home/kkumer/pype/data/gammastarp2gammap', approach=Approach.BMK))
#data.update(utils.loaddata('/home/kkumer/pype/data/gammastarp2gammap/EIC', approach=Approach.BMK))
#data.update(utils.loaddata('/home/kkumer/pype/data/ep2epgamma/EIC', approach=Approach.BMK))
db = shelve.open('/home/kkumer/pype/theories.db')
#dell = shelve.open('/home/kkumer/pype/dellB.db')

## [2] Choose subset of datapoints for fitting

##  #6 bin data
##  B0 = utils.select(data[67], criteria=['FTn == 0'])  # HERMES
##  B1 = utils.select(data[67], criteria=['FTn == 1'])  # HERMES
##  
##  # Local bins, 4 bins per dimension
##  # obsolete unpolarized target data, Morgan will provide update
##  #
##  A = utils.select(data[68], criteria=['FTn == -1'])  
##  # Faking 4 bins by dropping some points:
##  L4_ALUI = A[:1]+A[2:4]+A[5:7]+A[8:10]+A[11:13]+A[14:16]+A[17:]
##  #L4_AC_0 = utils.select(data[31], criteria=['FTn==0'])  
##  L4_AC_0 = B0[:1]+B0[2:4]+B0[5:7]+B0[8:10]+B0[11:13]+B0[14:16]+B0[17:]
##  #L4_AC_1 = utils.select(data[31], criteria=['FTn==1'])  
##  L4_AC_1 = B1[:1]+B1[2:4]+B1[5:7]+B1[8:10]+B1[11:13]+B1[14:16]+B1[17:]

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
bins = zip(L4_ALUI, L4_AC_0, L4_AC_1, L4_AUL, L4_ALL_0, 
        L4_ALL_1, L4_AUTI_1, L4_AUTI_0, L4_AUTI_m1)



## [4] Do the fit bin by bin

fl = open('one.res', 'w')
#for nextCFF in ['pReH', 'pImE', 'pReE', 'pImHt', 'pReHt', 'pImEt', 'pReEt']:
for nextCFF in ['pReE']:
    totchi = 0
    for nbin in range(len(bins)):
        # Model creation
        m = Model.ModelLocal()
        #m.parameters['pImH'] = 10. 
        th = Approach.BM10(m)
        th.name = 'C2-%s-bin%s' % (nextCFF[1:], nbin+1)
        #print 'Fitting %s ....\n' % th.name
        th.model.fix_parameters('ALL')
        th.model.release_parameters('pImH', nextCFF)
        f = Fitter.FitterMinuit(bins[nbin], th)
        #f.minuit.printMode = 2
        fl.write('\n')
        fl.write('-- Bin %2s: ' % (nbin+1,))
        sys.stdout.write('-- Bin %2s: \n' % (nbin+1,))
        f.fit()
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

    print '\n %s: TOTAL chisq = %f' % (nextCFF, totchi)
    fl.write('==== %s: TOTAL chisq = %f' % (nextCFF, totchi))
fl.close()


