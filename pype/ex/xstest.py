#!/usr/bin/env python
#

## Before running prepare files using, for bin 4,

# xs.exe 6 1 1 27.6 0.938 0.123  3.587 -0.408 33 > b4pp.dat
# xs.exe 6 -1 1 27.6 0.938 0.123  3.587 -0.408 33 > b4mp.dat
# xs.exe 6 1 -1 27.6 0.938 0.123  3.587 -0.408 33 > b4pm.dat
# xs.exe 6 -1 -1 27.6 0.938 0.123  3.587 -0.408 33 > b4mm.dat


import shelve, copy, sys, logging
from scipy.integrate import romb
import numpy as np


import Model, Approach, Fitter, Data, utils, plots
from constants import Mp, Mp2

data = utils.loaddata('/home/kkumer/gepard/pype/data/ep2epgamma', approach=Approach.BMK)  
db = shelve.open('/home/kkumer/gepard/pype/theories.db')

# Local 4-bin fits
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

th = db['KM12a']

nbin = 4
pp = np.loadtxt('b4pp.dat')
pm = np.loadtxt('b4pm.dat')
mp = np.loadtxt('b4mp.dat')
mm = np.loadtxt('b4mm.dat')

phis = pp[:,0]
dphi = 2*np.pi/32

pt3 = L4_AUL[3]  # 4th bin

unp = ((pp + pm) + (mp + mm))[:,1]
unp_positron = (pp + pm)[:,1]

ALUI = ((pp - pm) - (mp - mm))[:,1] / unp
AC = ((pp + pm) - (mp + mm))[:,1] / unp

AUL = -(pp + pm)[:,-1] / unp_positron 
ALL = -(pp - pm)[:,-1] / unp_positron 

AUTIcos = ((pp+pm) - (mp+mm))[:,2] / unp
AUTIsin = ((pp+pm) - (mp+mm))[:,3] / unp
AUTDVCScos = ((pp+pm) + (mp+mm))[:,2] / unp
AUTDVCSsin = ((pp+pm) + (mp+mm))[:,3] / unp

# ALUI
xsexe = romb(np.sin(phis)*ALUI, dx=dphi)/np.pi
true = th.predict(L4_ALUI[nbin-1], orig_conventions=True)
sys.stdout.write('%13s: % .3f, err = % .2f\n' % ('ALUI', xsexe, (xsexe-true)/true))
# AC0
xsexe = romb(AC, dx=dphi)/2./np.pi
true = th.predict(L4_AC_0[nbin-1], orig_conventions=True)
sys.stdout.write('%13s: % .3f, err = % .2f\n' % ('AC0', xsexe, (xsexe-true)/true))
# AC1
xsexe = romb(np.cos(phis)*AC, dx=dphi)/np.pi
true = th.predict(L4_AC_1[nbin-1], orig_conventions=True)
sys.stdout.write('%13s: % .3f, err = % .2f\n' % ('AC1', xsexe, (xsexe-true)/true))

# AUL
xsexe = romb(np.sin(phis)*AUL, dx=dphi)/np.pi
true = th.predict(L4_AUL[nbin-1], orig_conventions=True)
sys.stdout.write('%13s: % .3f, err = % .2f\n' % ('AUL', xsexe, (xsexe-true)/true))

# ALL
xsexe = romb(ALL, dx=dphi)/2./np.pi
true = th.predict(L4_ALL_0[nbin-1], orig_conventions=True)
sys.stdout.write('%13s: % .3f, err = % .2f\n' % ('ALL0', xsexe, (xsexe-true)/true))

xsexe = romb(np.cos(phis)*ALL, dx=dphi)/np.pi
true = th.predict(L4_ALL_1[nbin-1], orig_conventions=True)
sys.stdout.write('%13s: % .3f, err = % .2f\n' % ('ALL1', xsexe, (xsexe-true)/true))


# AUTIsin 
xsexe = romb(AUTIsin, dx=dphi)/2/np.pi
true = th.predict(L4_AUTI_0[nbin-1], orig_conventions=True)
sys.stdout.write('%13s: % .3f, err = % .2f\n' % ('AUTIsin_0', xsexe, (xsexe-true)/true))

xsexe = romb(np.cos(phis)*AUTIsin, dx=dphi)/np.pi
true = th.predict(L4_AUTI_1[nbin-1], orig_conventions=True)
sys.stdout.write('%13s: % .3f, err = % .2f\n' % ('AUTIsin_1', xsexe, (xsexe-true)/true))

# AUTIcos
xsexe = romb(np.sin(phis)*AUTIcos, dx=dphi)/np.pi
true = th.predict(L4_AUTI_m1[nbin-1], orig_conventions=True)
sys.stdout.write('%13s: % .3f, err = % .2f\n' % ('AUTIcos_m1', xsexe, (xsexe-true)/true))

# AUTDVCSsin
xsexe = romb(AUTDVCSsin, dx=dphi)/2/np.pi
true = th.predict(L4_AUTDVCS[nbin-1], orig_conventions=True)
sys.stdout.write('%13s: % .3f, err = % .2f\n' % ('AUTDVCSsin_0', xsexe, (xsexe-true)/true))
