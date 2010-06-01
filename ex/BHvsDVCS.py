"""Compares Bethe-Heitler and DVCS contributions to cross-section
for various experiments."""

import shelve, os, math

import numpy as np
#np.random.seed(68)

import matplotlib
if os.sys.platform == 'win32':
    matplotlib.use('WxAgg')
else: #linux
    matplotlib.use('TkAgg')
import pylab as plt


import Model, Approach, Fitter, Data
import utils 
import plots

from constants import Mp2, Mp

db = shelve.open('theories.db')
th = db['DMGLO1']

def tfromgamma(Q2, xB, gamma):
    """Calculate momentum transfer, given angle between photons."""
    aux1 = (2*Mp2 + Q2 + 2*(Mp + Q2/(2.*Mp*xB))*
            (Q2/(2.*Mp*xB) - math.sqrt(Q2 + Q2**2/(4.*Mp2*xB**2)) *
            math.cos(gamma)))
    aux2 = Mp + Q2/(2.*Mp*xB) - math.sqrt(
            Q2 + Q2**2/(4.*Mp2*xB**2))*math.cos(gamma)
    return  2*Mp2 - Mp * aux1 / aux2 


def tmin(Q2, xB):
    """BMK Eq. (31)"""
    eps2 = (2. * xB * Mp)**2 / Q2
    return -Q2 * ( 2. * (1.-xB)*(1. - math.sqrt(1.+eps2)) + eps2 ) / (
            4. * xB * (1.-xB) + eps2 )

def xsPart(part, xB, Q2, om1=5.75, om2=False, gamma=5., t=False):
    """Calculate BH or DVCS part of cross-section. 

    Either part='BH' or part='DVCS'. 
    gamma is angle between the photons. 
    Alternatively, momentum transfer t can be specified.
    om1 and om2 are energies of colliding particles.
    
    """
    pt = Data.DummyPoint()
    pt.in1particle = 'e'
    pt.in1charge = -1
    pt.xB = xB
    pt.Q2 = Q2
    pt.phi = np.pi
    pt.in1energy = om1
    if om2:
        pt.exptype = 'collider'
        pt.in2energy = om2
        pt.s = 2 * pt.in1energy * (pt.in2energy + math.sqrt(
            pt.in2energy**2 - Mp2)) + Mp2
    else:
        pt.exptype = 'fixed target'
        pt.s = 2 * Mp * pt.in1energy + Mp2
    if t:
        pt.t = t
    else:
        pt.t = tfromgamma(pt.Q2, pt.xB, gamma * np.pi/180.)
    utils.fill_kinematics(pt)
    pt.prepare(th)
    if part == 'BH':
        return th.PreFacSigma(pt) * th.TBH2unp(pt)
    elif part == 'DVCS':
        return th.PreFacSigma(pt) * th.TDVCS2unp(pt)

def gsubplot(ax, kin):
    """Plot BH-DVCS panel (gamma-dependence) for certain kinematics."""
    title, xB, Q2, E1, E2, lx, ly, gx, gy, gmax = kin
    ax.set_title(title)
    ax.set_yscale('log')  # y-axis to be logarithmic
    ax.set_xlabel('$\\theta_{\\gamma\\gamma}$ [deg]', fontsize=15)
    ax.set_ylabel('$\\sigma_{\\rm BH or DVCS}$ [nb/GeV$^4$]', fontsize=18)
    #ax.set_ylabel('$\\sigma_{\\rm DVCS}$ [nb/GeV$^4$]', fontsize=18)
    gammas = np.linspace(0.01, gmax, 150)
    bhs = np.array([xsPart('BH', xB, Q2, om1=E1, om2=E2, gamma=g) for g in gammas])
    dvcss = np.array([xsPart('DVCS', xB, Q2, om1=E1, om2=E2, gamma=g) for g in gammas])
    ax.plot(gammas, bhs, lw=1)
    ax.plot(gammas, dvcss, lw=2)
    if E2:
        ax.text(gx, gy, 'xB = %s, E1 = %s GeV^2, E2 = %s GeV^2' % (xB, E1, E2))
    else:
        ax.text(gx, gy, 'xB = %s, E1 = %s GeV^2' % (xB, E1))


def tsubplot(ax, kin):
    """Plot BH-DVCS panel (t-dependence) for certain kinematics."""
    title, xB, Q2, E1, E2, lx, ly, gx, gy, gmax = kin
    ax.set_title(title)
    ax.set_yscale('log')  # y-axis to be logarithmic
    ax.set_xlabel('$t$ [GeV$^2$]', fontsize=15)
    ax.set_ylabel('$\\sigma_{\\rm BH or DVCS}$ [nb/GeV$^4$]', fontsize=18)
    #ax.set_ylabel('$\\sigma_{\\rm DVCS}$ [nb/GeV$^4$]', fontsize=18)
    ts = np.linspace(tmin(Q2, xB), -1.0, 80)
    bhs = np.array([xsPart('BH', xB, Q2, om1=E1, om2=E2, t=itt) for itt in ts])
    dvcss = np.array([xsPart('DVCS', xB, Q2, om1=E1, om2=E2, t=itt) for itt in ts])
    ax.plot(ts, bhs, lw=1)
    ax.plot(ts, dvcss, lw=2)
    if E2:
        ax.text(lx, ly, 'xB = %s, E1 = %s GeV^2, E2 = %s GeV^2' % (xB, E1, E2))
    else:
        ax.text(lx, ly, 'xB = %s, E1 = %s GeV^2' % (xB, E1))

def ggrid(kin):
    """Produce DVCS grid (gamma-dependence) for certain kinematics."""
    title, xB, Q2, E1, E2, lx, ly, gx, gy, gmax = kin
    f = open(title.split()[0]+'-gamma.dat', 'w')
    gammas = np.linspace(0.01, gmax, 150)
    dvcss = np.array([xsPart('DVCS', xB, Q2, om1=E1, om2=E2, gamma=g) for g in gammas])
    z = zip(gammas, dvcss)
    for it in z:
        f.write('%s  %s\n' % it)
    f.close()

def tgrid(kin):
    """Produce DVCS grid (t-dependence) for certain kinematics."""
    title, xB, Q2, E1, E2, lx, ly, gx, gy, gmax = kin
    f = open(title.split()[0]+'-t.dat', 'w')
    ts = np.linspace(tmin(Q2, xB), -1.0, 100)
    dvcss = np.array([xsPart('DVCS', xB, Q2, om1=E1, om2=E2, t=itt) for itt in ts])
    z = zip(ts, dvcss)
    for it in z:
        f.write('%s  %s\n' % it)
    f.close()

# [expname, xB, Q2, E1, tlabelx, tlabely, glabelx, glabely, gmax]
kin = ( ['HALL-A (JLAB)', 0.35, 2.5, 5.75, 0, -0.7, 6, 20, 44, 46],
        ['HERMES (DESY)', 0.1, 2.5, 27.6, 0,  -0.7, 6, 8, 200, 20],
        ['COMPASS (CERN)', 0.05, 2.5, 160, 0,  -0.7, 0.6, 4.5, 200, 10],
        ['EIC (BNL)', 0.001, 2.5, 5., 250, -0.98, 500, 0.1, 1000, 0.4] )

fig = plt.figure()
title = 'DVCS (at Q2=2.5 GeV^2)'
fig.canvas.set_window_title(title)
fig.suptitle(title)
for n in range(4):
    ax = fig.add_subplot(2,2,n+1)
    tsubplot(ax, kin[n])
    #tgrid(kin[n])
fig.canvas.draw()
fig.show()

fig = plt.figure()
title = 'DVCS (at Q2=2.5 GeV^2)'
fig.canvas.set_window_title(title)
fig.suptitle(title)
for n in range(4):
    ax = fig.add_subplot(2,2,n+1)
    gsubplot(ax, kin[n])
    #ggrid(kin[n])
fig.canvas.draw()
fig.show()
