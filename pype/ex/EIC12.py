#!/usr/bin/env python

# Figs for AFKM12 paper with EIC simulations

import sys, shelve, matplotlib

import matplotlib.pyplot as plt
import numpy as np

import Model, Approach, Data
import utils 
from constants import Mp, Mp2

from math import sqrt

#import logging
#logging.basicConfig(level=logging.DEBUG)

db = shelve.open('/home/kkumer/pype/theories.db')
thAFKM12 = db['AFKM12']
Model.ComptonGepard.gepardPool.pop()
thKM10 = db['KM10']
Model.ComptonGepard.gepardPool.pop()
theories = [thAFKM12, thKM10]
#theories = [thAFKM12]

## collider datapoint
#def ptcol(th, Q=-1, pol=0, Ee=20, Ep=250, xB=0.001998, Q2=4., t=-0.1, 
#        phi=np.pi, varphi=-np.pi/2., FTn=None):
##def ptcol(th, Q=-1, pol=0, Ee=20, Ep=250, xB=0.002, Q2=7.3, t=-0.275, 
##        phi=np.pi, varphi=-np.pi/2., FTn=None):
#    ptc = Data.DummyPoint()
#    ptc.in1energy = Ee
#    ptc.in2energy = Ep
#    ptc.s = 2 * ptc.in1energy * (ptc.in2energy + sqrt(
#                ptc.in2energy**2 - Mp2)) + Mp2
#    ptc.in1charge = Q
#    ptc.in1polarization = pol
#    ptc.in2polarizationvector = 'T'
#    ptc.in2polarization = 1 # relevant only for XLP and TSA
#    ptc.xB = xB
#    ptc.Q2 = Q2
#    ptc.t = t
#    ptc.xi = ptc.xB/(2.-ptc.xB)
#    if phi:
#        ptc.phi = phi
#        ptc.units = {'phi': 'radian'}
#    elif FTn:
#        ptc.FTn = FTn
#    ptc.varphi = varphi
#    ptc.frame = 'Trento'
#    utils.fill_kinematics(ptc)
#    th.to_conventions(ptc)
#    th.prepare(ptc)
#    return ptc

#ptI = ptcol(thKM10, Q=1, pol=0, Ee=5, Ep=100, phi=np.pi, Q2=4.4, xB=5.1e-3, t=-0.25)
#ptII = ptcol(thKM10, Q=1, pol=0, Ee=20, Ep=250, phi=np.pi, Q2=4.4, xB=5.1e-4, t=-0.25)



fig = plt.figure()
fig.suptitle('Fig. 12 ')
fig.subplots_adjust(bottom=0.1, right=0.9, wspace=0., hspace=0)

pnrow = 0
for kins in [(5, 100, 5.1e-3), (20, 250, 5.1e-4)]:
    in1energy, in2energy, xB = kins
    pt = Data.DummyPoint()
    pt.exptype = 'collider'
    pt.in1particle = 'e'
    pt.in1charge = 1
    pt.in1energy = in1energy
    pt.in1polarization = 0
    pt.in2particle = 'p'
    pt.in2energy = in2energy
    pt.s = 2 * pt.in1energy * (pt.in2energy + sqrt(
        pt.in2energy**2 - Mp2)) + Mp2
    #pt.xB = 5.1e-3
    pt.xB = xB
    pt.t = -0.25
    pt.Q2 = 4.4
    pt.units = {'phi' : 'radian'}
    utils.fill_kinematics(pt)
    thKM10.__class__.to_conventions(pt)
    thKM10.__class__.prepare(pt)
    #
    ###  phi panel
    ax = fig.add_subplot(2,3,pnrow+1)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    phi = np.linspace(0., 2*np.pi, 50)
    linestyles = ['b-', 'r--', 'g-', 'b-.', 'p:']
    dasheslengths = [(None, None), (20,5), (20,5,5,5), (5,5)]
    #labels = [r'$\sigma^{\uparrow}$ ' + t.name for t in lines]
    pn = 0
    # Must go to BKM explicitely here
    ln = 0
    for th in theories:
        line = th.BCA(pt, vars={'phi':np.pi-phi})
        ax.plot(phi, line, linestyles[ln], dashes=dasheslengths[ln], linewidth=2) 
        ln +=1
    ax.set_xlim(0.0, 2*np.pi)
    ax.set_ylim(-0.45, 0.45)
    ax.set_ylabel('$A_{\\rm C}$')
    ax.text(1, 0.35, "$t = %s \\,{\\rm GeV}^2$" % pt.t)
    ax.text(0.8, -0.17, "$E_e = %d \\,{\\rm GeV}$" % pt.in1energy)
    ax.text(0.8, -0.26, "$E_p = %d \\,{\\rm GeV}$" % pt.in2energy)
    ax.text(0.8, -0.35, "$x_B = %s$" % pt.xB)
    ax.text(0.8, -0.44, "$Q^2 = %s \\,{\\rm GeV}^2$" % pt.Q2)
    # axes labels
    if pnrow > 2:
        ax.set_xlabel('$\\phi\\quad {\\rm [rad]}$')
        ax.set_xticklabels(['$0$', '$\\pi/2$', '$\\pi$', '$3 \\pi/2$'])
    else:
        ax.set_xticklabels([])
    ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2])
    #
    ###  t panel
    pt.FTn = 1
    ax = fig.add_subplot(2,3,pnrow+2)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    tms = list(np.linspace(0.,0.1))+list(np.linspace(0.1, 0.85))
    ln = 0
    for th in theories:
        line = []
        for tm in tms:
            pt.t = -tm
            # Change sign to go to Trento
            line.append(-th.BCA(pt))
        ax.plot(tms, line, linestyles[ln], dashes=dasheslengths[ln], linewidth=2, label=th.name) 
        ln += 1
    ax.set_xlim(0., 0.9)
    ax.set_ylim(-0.45, 0.45)
    if pnrow > 2:
        ax.set_xlabel('$-t \\,{\\rm [GeV}^2{\\rm ]}$')
    else:
        ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.2))  # tickmarks
    # Legend
    leg = ax.legend(loc='upper center', handlelength=4.0, fancybox=True)
    frame  = leg.get_frame()
    frame.set_facecolor('0.90')    # set the frame face color to light gray
    for t in leg.get_texts():
        t.set_fontsize(10)    # the legend text fontsize
    #
    ###  xB panel
    pt.t = -0.25
    pt.tm = 0.25
    ax = fig.add_subplot(2,3,pnrow+3)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    xBs = np.logspace(-4., -0.7, 40)
    line = []
    pn = 2
    pt.Q2 = 2.5
    for xB in xBs:
        del pt.W
        pt.xB = xB
        utils._complete_xBWQ2(pt)
        line.append(-th.BCA(pt))
    ax.semilogx(xBs, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2,
            label="$Q^2= %s \\,{\\rm GeV}^2$" % pt.Q2) 
    pn += 1
    line = []
    pt.Q2 = 13.9
    for xB in xBs:
        del pt.W
        pt.xB = xB
        utils._complete_xBWQ2(pt)
        line.append(-th.BCA(pt))
    ax.semilogx(xBs, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2,
            label="$Q^2= %s \\,{\\rm GeV}^2$" % pt.Q2) 
    ax.set_xlim(5.e-5, 0.2)
    ax.set_ylim(-0.45, 0.45)
    if pnrow > 2:
        ax.set_xlabel('$x_{\\rm B}$')
    else:
        ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.text(0.001, 0.35, "$t = %s \\,{\\rm GeV}^2$" % pt.t)
    # Legend
    leg = ax.legend(loc='upper center', handlelength=3.0, fancybox=True)
    frame  = leg.get_frame()
    frame.set_facecolor('0.90')    # set the frame face color to light gray
    for t in leg.get_texts():
        t.set_fontsize(10)    # the legend text fontsize
    for l in leg.get_lines():
        l.set_linewidth(2.0)  # the legend line width
    pnrow += 3
fig.savefig('Fig12.png')


fig = plt.figure()
fig.suptitle('Fig. 14 ')
fig.subplots_adjust(bottom=0.1, top=1., right=0.9, wspace=0.)
pt = Data.DummyPoint()
pt.exptype = 'collider'
pt.in1particle = 'e'
pt.in1charge = -1
pt.in1energy = 20.
pt.in1polarization = 0.0
pt.in2particle = 'p'
pt.in2energy = 250.
pt.in2polarizationvector = 'T'
pt.in2polarization = 1
pt.s = 2 * pt.in1energy * (pt.in2energy + sqrt(
    pt.in2energy**2 - Mp2)) + Mp2
kins = [(14., 8.2e-4), (14., 1.3e-3), (14., 5.1e-3),
        (7.8, 5.1e-4), (7.8, 8.2e-4), (7.8, 5.1e-3),
        (4.4, 3.2e-4), (4.4, 8.2e-4), (4.4, 5.1e-3)]
pt.varphi = -np.pi/2
pt.units = {'phi' : 'radian'}
linestyles = ['b-', 'g-.', 'r--',]
dasheslengths = [(None, None), (20,5,5,5), (20,5)]
pn = 1
for kin in kins:
    if hasattr(pt, 'W'): del pt.W
    pt.Q2, pt.xB = kin
    if hasattr(pt, 'tm'): del pt.tm
    pt.t = -0.25
    ax = fig.add_subplot(3,3,pn)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    phi = np.linspace(0., 2*np.pi, 50)
    utils.fill_kinematics(pt)
    #labels = [r'$\sigma^{\uparrow}$ ' + t.name for t in lines]
    thAFKM12.__class__.to_conventions(pt)
    thAFKM12.__class__.prepare(pt)
    ln = 0
    for kappa in [1.5, 0., -1.5]:
        thAFKM12.m.parameters['KAPS'] = kappa
        # Must go to BKM explicitely here
        line = thAFKM12.TSA(pt, vars={'phi':np.pi-phi})
        ax.plot(phi, line, linestyles[ln], dashes=dasheslengths[ln], linewidth=2, 
                label='$\\kappa^{\\rm sea} = %.1f$' % kappa) 
        ln +=1
    ax.set_xlim(0.0, 2*np.pi)
    ax.set_ylim(-0.9, 0.9)
    # axes labels
    if pn > 6:
        ax.set_xlabel('$\\phi\\quad {\\rm [rad]}$')
    if (pn % 3) == 1:
        ax.set_ylabel('$A_{\\rm UT}^{\\sin(\\phi - \\phi_S)}$')
    #ax.text(1, 0.35, "$t = %s \\,{\\rm GeV}^2$" % pt.t)
    ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2])
    ax.set_xticklabels(['$0$', '$\\pi/2$', '$\\pi$', '$3 \\pi/2$'])
    # Legend
    if pn == 1:
        leg = ax.legend(loc='upper left', handlelength=3.5,
                labelspacing=0.1, fancybox=True)
        frame  = leg.get_frame()
        frame.set_facecolor('0.90')    # set the frame face color to light gray
        for t in leg.get_texts():
            t.set_fontsize(8)    # the legend text fontsize
        for l in leg.get_lines():
            l.set_linewidth(1.5)  # the legend line width
    #
    #ax.text(0.1, -0.17, "$E_e = %d \\,{\\rm GeV}$" % pt.in1energy, fontsize=14)
    #ax.text(0.1, -0.26, "$E_p = %d \\,{\\rm GeV}$" % pt.in2energy, fontsize=14)
    ax.text(0.1, -0.55, "$x_B = %s$" % pt.xB)
    ax.text(0.1, -0.75, "$Q^2 = %s \\,{\\rm GeV}^2$" % pt.Q2)
    pn += 1
fig.savefig('Fig14.png')
