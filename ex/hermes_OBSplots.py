#!/usr/bin/env python
""" 
Preparing data filesfor grace plots of HERMES observables and local fits 
to them.
"""


import shelve, copy, sys, logging

import numpy as np
import scipy.stats

logging.basicConfig(level=logging.WARNING)

import Model, Approach, Fitter, Data, utils, plots
from constants import Mp, Mp2

from results import *
from math import sqrt


## [1] Load experimental data and theoretical models

data = utils.loaddata('/home/kkumer/pype/data/ep2epgamma', approach=Approach.BMK)  
db = shelve.open('/home/kkumer/pype/theories.db')
thKM12a = db['KM12a']

## [2] Choose datapoints for fitting

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

obsnames = ['ALUI', 'AC0', 'AC1', 'AUL', 'ALL0', 'ALL1', 
        'AUTI1', 'AUTI0', 'AUTIm1', 'AUTDVCS']

def Hptval(pt):
    """pt.val in HERMES conventions"""
    # cos-varphi and sin-varphi harmonics change sign
    if pt.has_key('varFTn'):
        val = -pt.val
    else:
        val = pt.val
    # All cos-phi (positive FTn) harmonics change sign
    return np.sign(0.2-pt.FTn)*val

def printset(filename, set, nbin=4):
    f = open(filename+'-tm.dat', 'w')
    for pt in set[:nbin]:
        f.write("%f  %f  %f\n" % (pt.tm, Hptval(pt), pt.err))
    f.write('\n')
    f.close()
    f = open(filename+'-xB.dat', 'w')
    for pt in set[nbin:2*nbin]:
        f.write("%f  %f  %f\n" % (pt.xB, Hptval(pt), pt.err))
    f.write('\n')
    f.close()
    f = open(filename+'-Q2.dat', 'w')
    for pt in set[2*nbin:]:
        f.write("%f  %f  %f\n" % (pt.Q2, Hptval(pt), pt.err))
    f.write('\n')
    f.close()

printset('fits/HERMES-ALUI', L4_ALUI)
printset('fits/HERMES-AC_0', L4_AC_0)
printset('fits/HERMES-AC_1', L4_AC_1)

printset('fits/HERMES-AUL', L4_AUL)
printset('fits/HERMES-ALL0', L4_ALL_0)
printset('fits/HERMES-ALL1', L4_ALL_1)

printset('fits/HERMES-AUTI1', L4_AUTI_1)
printset('fits/HERMES-AUTI0', L4_AUTI_0)
printset('fits/HERMES-AUTIm1', L4_AUTI_m1)
printset('fits/HERMES-AUTDVCS', L4_AUTDVCS)


tmshifts = [0.02, 0.03, 0]
xBshifts = [0.01, 0.015, 0]
Q2shifts = [0.2, 0.3, 0]
shifts = zip(tmshifts, xBshifts, Q2shifts)
## [3] Iterate over all models:
for prefix in ['', 'Ht-', 'KM12a-']:
    sys.stdout.write('\n Doing %smodel.' % (prefix,))
    fls = []
    for obs in  obsnames:
        fls.append([open('fits/%s%s-tm.dat' % (prefix, obs), 'w'),
                    open('fits/%s%s-xB.dat' % (prefix, obs), 'w'),
                    open('fits/%s%s-Q2.dat' % (prefix, obs), 'w')])
    ##    Do the fit bin by bin
    totchi = 0
    sys.stdout.write('\n Doing bin ..')
    for nbin in range(len(bins)):
        sys.stdout.write('.%i.' % (nbin+1,))
        sys.stdout.flush()
        pts = bins[nbin]
        if prefix == '':
            # -- Fitting to ImH and ReE
            m = Model.ModelLocal()
            m.parameters['pImH'] = 10. 
            th = Approach.BM10(m)
            th.model.fix_parameters('ALL')
            th.model.release_parameters('pImH', 'pReE')
            f = Fitter.FitterMinuit(pts, th)
            f.minuit.printMode = 0
            f.fit()
            tmshift, xBshift, Q2shift = shifts[0]
        elif prefix == 'Ht-':
            # -- Fitting to ImH and ReH
            m = Model.ModelLocal()
            m.parameters['pImH'] = 10. 
            th = Approach.BM10(m)
            th.model.fix_parameters('ALL')
            th.model.release_parameters('pImH', 'pReH', 'pImHt')
            f = Fitter.FitterMinuit(pts, th)
            f.minuit.printMode = 0
            f.fit()
            tmshift, xBshift, Q2shift = shifts[1]
        else:
            ## -- KM12a global fit
            #mGepard = Model.ComptonGepard(cutq2=0.5)
            #mDRPPsea = Model.ComptonModelDRPPsea()
            #m = Model.HybridDipole(mGepard, mDRPPsea)
            #th = Approach.BM10(m)
            #th.name = 'KM12a'
            #g = th.m.g
            #th.m.parameters.update(GLOnoBSSandBSS)
            #th.m.covariance = GLOnoBSSandBSScov
            th = thKM12a
            tmshift, xBshift, Q2shift = shifts[2]
        totchi += th.chisq(pts)[0]
        xshift = 0
        if nbin<4:
            files = [ft for ft, fx, fQ in fls]
            x = 'tm'
            xshift = tmshift
        elif nbin<8:
            files = [fx for ft, fx, fQ in fls]
            x = 'xB'
            xshift = xBshift
        else:
            files = [fQ for ft, fx, fQ in fls]
            x = 'Q2'
            xshift = Q2shift
        for f, obs, nobs in zip(files, obsnames, range(len(obsnames))):
            val, err = th.predict(pts[nobs], error=True, orig_conventions=True)
            f.write('%f  %f  %f\n' % (getattr(pts[nobs], x)+xshift, val, err))

    for ft, fx, fQ  in fls:
        ft.close()
        fx.close()
        fQ.close()


    print '%sTOTAL chisq = %f' % (prefix, totchi)

