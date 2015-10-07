#!/usr/bin/env python

# Figs for EIC proceedings

import sys

import matplotlib.pyplot as plt
import numpy as np

import Model, Approach, Data
import utils 

from constants import Mp, Mp2

from results import *
from math import sqrt


## Create a theory

# Gepard only
mGepard = Model.ComptonGepard(cutq2=0.5)


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


# Hybrid model

mDRPPsea = Model.ComptonModelDRPPsea()
m = Model.HybridKelly(mGepard, mDRPPsea)
th = Approach.BM10(m)
th.name = 'KM10b'
th.m.parameters.update(KM10b)  
g = th.m.g


    

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




fig = plt.figure()
fig.suptitle('Fig. 8 (hotfixedBMK)')
# Left panel
ax = fig.add_subplot(1,2,1)
phis = np.linspace(0.001, 6.282)
bsas =  [th.BSA(ptcol(th, Ee=5, Ep=250, Q2=10, phi=f, xB=5e-3)) for f in phis]
ax.plot(phis, bsas, color='blue', label='KM10b')
ax.set_xlabel('$\\phi$')
ax.set_ylabel('$A_{BS}$')
ax.set_ylim(-0.42, 0.5)
ax.axhline(y=0, linewidth=1, color='g')
ax.text(1, -0.25, 'Ee=5, Ep=250, Q2=10')
ax.text(1, -0.35, 't=-0.2, xB=5e-3')
ax.legend(handlelen=0.15, handletextsep = 0.04).draw_frame(0)
# Right panel
ax = fig.add_subplot(1,2,2)
xs = np.logspace(-4., -1.8)
xsL = np.logspace(-2., -0.28, 120)
ys =  [th.BSA(ptcol(th, Ee=30, Ep=360, phi=None, FTn=-1, xB=x)) for x in xs]
ysL = [th.BSA(ptcol(th, Ee=5, Ep=150, Q2=50, phi=None, FTn=-1, xB=x)) for x in xsL]
ax.plot(xs, ys, color='blue', label='KM10b')
ax.plot(xsL, ysL, color='red', label='KM10b')
ax.set_xscale('log')
ax.set_xlabel('$x_B$')
ax.set_ylabel('$A^{(1)}_{BS}$')
ax.set_ylim(-0.42, 0.5)
ax.axhline(y=0, linewidth=1, color='g')
ax.text(0.0001, -0.25, 'Ee=30, Ep=360, Q2=4')
ax.text(0.01, -0.35, 'Ee=5, Ep=150, Q2=50')
ax.legend(handlelen=0.15, handletextsep = 0.04).draw_frame(0)
fig.savefig('Fig8.png')

fig = plt.figure()
fig.suptitle('Fig. 10 (hotfixedBMK)')
# Left panel
ax = fig.add_subplot(1,2,1)
phis = np.linspace(0.001, 6.282)
bcas =  [th.BCA(ptcol(th, Q=1, pol=0, Ee=5, Ep=250, Q2=10, phi=f, xB=5e-3)) for f in phis]
bcas4 =  [th.BCA(ptcol(th, Q=1, pol=0, Ee=5, Ep=250, Q2=4, phi=f, xB=5e-3)) for f in phis]
ax.plot(phis, bcas, color='blue', label='Q2 = 10 GeV^2')
ax.plot(phis, bcas4, 'r-', label='Q2 = 4 GeV^2')
ax.set_xlabel('$\\phi$')
ax.set_ylabel('$A_{BC}$')
ax.set_ylim(-0.42, 0.42)
ax.axhline(y=0, linewidth=1, color='g')
ax.text(1, -0.25, 'Ee=5, Ep=250')
ax.text(1, -0.35, 't=-0.2, xB=5e-3')
ax.legend(handlelen=0.15, handletextsep = 0.04).draw_frame(0)
# Right panel
ax = fig.add_subplot(1,2,2)
xs = np.logspace(-4., -1.8)
xsL = np.logspace(-2., -0.28, 120)
ys =  [th.BCA(ptcol(th, Q=1, pol=0, Ee=30, Ep=360, phi=None, FTn=1, xB=x)) for x in xs]
ysL = [th.BCA(ptcol(th, Q=1, pol=0, Ee=5, Ep=150, Q2=50, phi=None, FTn=1, xB=x)) for x in xsL]
ax.plot(xs, ys, color='blue', label='KM10b')
ax.plot(xsL, ysL, color='red', label='KM10b')
ax.set_xscale('log')
ax.set_xlabel('$x_B$')
ax.set_ylabel('$A^{(1)}_{BC}$')
ax.set_ylim(-0.42, 0.42)
ax.axhline(y=0, linewidth=1, color='g')
ax.text(0.0001, -0.25, 'Ee=30, Ep=360, Q2=4')
ax.text(0.01, -0.35, 'Ee=5, Ep=150, Q2=50')
ax.legend(handlelen=0.15, handletextsep = 0.04).draw_frame(0)
fig.savefig('Fig10.png')

fig = plt.figure()
fig.suptitle('Fig. 9 (BM10)')
# Left panel
ax = fig.add_subplot(1,2,1)
phis = np.linspace(0.001, 6.282)
ysKM09a = [tDRBM.TSA(ptcol(tDRBM, pol=0, Ep=150, xB=0.1, phi=f)) for f in phis]
ysKM09b = [tDR1BM.TSA(ptcol(tDR1BM, pol=0, Ep=150, xB=0.1, phi=f)) for f in phis]
ysKM10b = [th.TSA(ptcol(th, pol=0, Ep=150, xB=0.1, phi=f)) for f in phis]
ax.plot(phis, ysKM09a, 'g--', label='KM09a')
ax.plot(phis, ysKM09b, 'b', label='KM09b')
ax.plot(phis, ysKM10b, 'r-.', label='KM10b')
ax.set_xlabel('$\\phi$')
ax.set_ylabel('$A_{TS}$')
ax.set_ylim(-0.22, 0.22)
ax.axhline(y=0, linewidth=1, color='g')
ax.text(1, -0.15, 'Ee=5, Ep=150')
ax.text(1, -0.20, 'Q2=4, t=-0.2, xB=0.1')
ax.legend(handlelen=0.15, handletextsep = 0.04).draw_frame(0)
# Right panel
ax = fig.add_subplot(1,2,2)
ysKM09aR = [tDRBM.TSA(ptcol(tDRBM, pol=0, Ep=350, xB=0.01, phi=f)) for f in phis]
ysKM09bR = [tDR1BM.TSA(ptcol(tDR1BM, pol=0, Ep=350, xB=0.01, phi=f)) for f in phis]
ysKM10bR = [th.TSA(ptcol(th, pol=0, Ep=350, xB=0.01, phi=f)) for f in phis]
ax.plot(phis, ysKM09aR, 'g--', label='KM09a')
ax.plot(phis, ysKM09bR, 'b', label='KM09b')
ax.plot(phis, ysKM10bR, 'r-.', label='KM10b')
ax.set_xlabel('$\\phi$')
ax.set_ylim(-0.22, 0.22)
ax.axhline(y=0, linewidth=1, color='g')
ax.text(1, -0.15, 'Ee=5, Ep=350')
ax.text(1, -0.20, 'Q2=4, t=-0.2, xB=0.01')
ax.legend(handlelen=0.15, handletextsep = 0.04).draw_frame(0)
fig.savefig('Fig9.png')

