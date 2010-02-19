""" 
Auxilliary plotting functions written for ad-hoc plotting
of some specific datasets. Note that actual plotting algorithm
is in `utils` module.

plotHALLA -- plots HALL-A data

"""


import os

import numpy as np
import pylab as plt
from matplotlib.ticker import MultipleLocator

import Data
import utils

def plotHERMES(data, fits=[], path=None, fmt='png'):
    """Makes plot of HERMES preliminary BCA and BSA data with fit line defined by pars"""

    #2: "data/ep2epgamma-BCA-HERMES-08-cos1_b.dat", 
    #4: "data/ep2epgamma-BCA-HERMES-08-cos0_b.dat",
    #5: "data/ep2epgamma-ALU-HERMES-08-Isin1.dat"]
    ids = [2, 4, 5]
    title = '' #'HERMES-08'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    xaxes = ['mt', 'xB', 'Q2']
    ylims = [(-0.05, 0.24), (-0.15, 0.05), (-0.30, 0.05)]
    # we have 3x18=54 points to be separated in nine panels six points each:
    for y in range(3):
        for x in range(3):
            panel = 3*y + x + 1  # 1, 2, ..., 9
            ax = fig.add_subplot(3,3,panel)
            ax.yaxis.set_major_locator(MultipleLocator(0.1))  # tickmarks
            utils.subplot(ax, data[ids[y]][x*6:x*6+6], xaxes[x], [], fits)
            apply(ax.set_ylim, ylims[y])
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def plotCLAS(data, fits=[], path=None, fmt='png'):
    """Makes plot of CLAS BSA data with fit line defined by pars"""

    #datafile = "data/ep2epgamma-ALU-CLAS_KK-07.dat" # id = 25
    dataset = data[25]
    title = '' # 'CLAS-07'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    nmax = len(dataset) - 1
    # Each different Q2 has its own panel
    panel = 1
    npt = 0
    panelpoints = []
    for panel in [7, 8, 9, 4, 5, 6, 1, 2, 3]:
        ax = fig.add_subplot(3,3,panel)
        panelpoints = []
        pt = dataset[npt]
        Q2 = pt.Q2
        while Q2 == pt.Q2:
            panelpoints.append(pt)
            npt += 1
            if npt == nmax:
                break
            pt = dataset[npt]
        # now transform list into DataSet instance ...
        panelset = Data.DataSet(panelpoints)
        panelset.__dict__ = dataset.__dict__.copy()
        # ... and plot
        utils.subplot(ax, panelset, 'mt', ['Q2', 'xB'], fits)
        plt.xlim(0.0, 0.6)
        plt.ylim(0.0, 0.4)
        ax.yaxis.set_major_locator(MultipleLocator(0.1))
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def plotHALLA(data, fits=[], path=None, fmt='png'):
    """Makes plot of HALL-A data with fit line defined by pars"""

    ids = [9, 14, 20, 21, 23, 24]
    title = '' # 'HALLA-06'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    panel = 1
    for id in ids:
        ax = fig.add_subplot(2,3,panel)
        utils.subplot(ax, data[id], 'phi', ['Q2', 't'], fits)
        ax.xaxis.set_major_locator(MultipleLocator(120))
        panel += 1
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def plotCOMPASS(ff, fits=[], path=None, fmt='png'):
    """Makes plot of COMPASS BCS asymmetry and difference and summ of xs for FormFactors model ff."""

    title = 'COMPASS BCSA (asymmetry), BCSD (difference) and BCSS (sum) of XS'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    #fig.subplots_adjust(bottom=0.35)
    # Asymmetry panel
    pt = Data.DummyPoint()
    pt.exptype = 'fixed target'
    pt.in1energy = 160.
    pt.xB = 0.05
    pt.t = -0.2
    pt.Q2 = 2.
    ax = fig.add_subplot(2,2,1)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    phi = np.arange(0., np.pi, 0.2)
    linestyles = ['g--', 'b-', 'r-.']
    labels = ['HERMES+CLAS', 'HERMES+CLAS+HALLA', '+HALLA(phi)']
    #labels = ['GLO1 (DM)', 'GLO1 (KK)', '']
    pn = 0
    for (approach, pars) in fits:
        pt.prepare(approach)
        line = approach.BCSA(pt, 0.8, pars, {'phi':np.pi - phi})
        ax.plot(phi, line, linestyles[pn], linewidth=2, label=labels[pn]) 
        pn += 1
    #ax.set_ylim(0.0, 0.5)
    # axes labels
    ax.set_xlabel('$\\phi$')
    ax.set_ylabel('BCS Asymmetry')
    ax.legend()
    # Difference  panel
    pt = Data.DummyPoint()
    pt.exptype = 'fixed target'
    pt.in1energy = 160.
    pt.xB = 0.05
    pt.t = -0.2
    pt.Q2 = 2.
    ax = fig.add_subplot(2,2,2)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    phi = np.arange(0., np.pi, 0.2)
    linestyles = ['g--', 'b-', 'r-.']
    labels = ['HERMES+CLAS', 'HERMES+CLAS+HALLA', '+HALLA(phi)']
    #labels = ['GLO1 (DM)', 'GLO1 (KK)', '']
    pn = 0
    for (approach, pars) in fits:
        pt.prepare(approach)
        # nb converted to pb:
        line = 1000. * approach.BCSD(pt, 0.8, pars, {'phi':np.pi - phi})
        ax.plot(phi, line, linestyles[pn], linewidth=2, label=labels[pn]) 
        pn += 1
    #ax.set_ylim(0.0, 0.5)
    # axes labels
    ax.set_xlabel('$\\phi$')
    ax.set_ylabel('BCS Difference  [pb]')
    ax.legend()
    # Sum pannel
    pt = Data.DummyPoint()
    pt.exptype = 'fixed target'
    pt.in1energy = 160.
    pt.xB = 0.05
    pt.t = -0.2
    pt.Q2 = 2.
    ax = fig.add_subplot(2,2,3)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    phi = np.arange(0., np.pi, 0.2)
    linestyles = ['g--', 'b-', 'r-.']
    labels = ['HERMES+CLAS', 'HERMES+CLAS+HALLA', '+HALLA(phi)']
    #labels = ['GLO1 (DM)', 'GLO1 (KK)', '']
    pn = 0
    for (approach, pars) in fits:
        pt.prepare(approach)
        # nb converted to pb:
        line = 1000. * approach.BCSS(pt, 0.8, pars, {'phi':np.pi - phi})
        ax.plot(phi, line, linestyles[pn], linewidth=2, label=labels[pn]) 
        pn += 1
    #ax.set_ylim(0.0, 0.5)
    # axes labels
    ax.set_xlabel('$\\phi$')
    ax.set_ylabel('BCS Difference  [pb]')
    ax.legend()
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig


def plotHBCSA(ff, fits=[], path=None, fmt='png'):
    """Makes plot of Im(cffH) and COMPASS BCSA for FormFactors model ff."""

    title = '' # 'Fig 15'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    fig.subplots_adjust(bottom=0.65)
    # Left panel
    pt = Data.DummyPoint()
    ax = fig.add_subplot(1,2,1)
    ax.set_xscale('log')  # x-axis to be logarithmic
    xval = np.power(10., np.arange(-3.5, 0, 0.01)) 
    linestyles = ['g--', 'b-', 'r-.']
    pn = 0
    for (approach, pars) in fits:
        pt.t = 0.0
        # kludge alert!
        line1 = xval * ff.ImH(pt, pars, xval) / np.pi
        pt.t = -0.3
        line2 = xval * ff.ImH(pt, pars, xval) / np.pi
        ax.plot(xval, line1, linestyles[pn])
        ax.plot(xval, line2, linestyles[pn], linewidth=2) 
        pn += 1
    ax.set_ylim(0.0, 0.5)
    ax.set_xlim(0.0005, 1.0)
    #plt.ylim(0.0, 0.5)
    # axes labels
    ax.set_xlabel('$x$')
    ax.set_ylabel('$x H(x, x, t)$')
    # Right panel
    pt = Data.DummyPoint()
    pt.exptype = 'fixed target'
    pt.in1energy = 160.
    pt.xB = 0.05
    pt.t = -0.2
    pt.Q2 = 2.
    ax = fig.add_subplot(1,2,2)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    phi = np.arange(0., np.pi, 0.2)
    linestyles = ['g--', 'b-', 'r-.']
    labels = ['HERMES+CLAS', 'HERMES+CLAS+HALLA', '+HALLA(phi)']
    #labels = ['GLO1 (DM)', 'GLO1 (KK)', '']
    pn = 0
    for (approach, pars) in fits:
        pt.prepare(approach)
        line = approach.BCSA(pt, np.pi - phi, pars)
        ax.plot(phi, line, linestyles[pn], linewidth=2, label=labels[pn]) 
        pn += 1
    #ax.set_ylim(0.0, 0.5)
    # axes labels
    ax.set_xlabel('$\\phi$')
    ax.set_ylabel('BCSA')
    ax.legend()
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def plotH(ff, fits=[], path=None, fmt='png'):
    """Makes plot of Im(cffH)."""

    title = '' # 'Fig 15'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    fig.subplots_adjust(bottom=0.65)
    # Left panel
    pt = Data.DummyPoint()
    ax = fig.add_subplot(1,1,1)
    ax.set_xscale('log')  # x-axis to be logarithmic
    xval = np.power(10., np.arange(-3.5, 0, 0.01)) 
    linestyles = ['g--', 'b-', 'r-.']
    labels = ['HERMES+CLAS', 'HERMES+CLAS+HALLA', '+HALLA(phi)']
    pn = 0
    for (approach, pars) in fits:
        pt.t = 0.0
        # kludge alert!
        line1 = xval * ff.ImH(pt, pars, xval) / np.pi
        pt.t = -0.3
        line2 = xval * ff.ImH(pt, pars, xval) / np.pi
        ax.plot(xval, line1, linestyles[pn], linewidth=2)
        ax.plot(xval, line2, linestyles[pn], linewidth=4, label=labels[pn]) 
        pn += 1
    ax.set_ylim(0.0, 0.5)
    ax.set_xlim(0.0005, 1.0)
    #plt.ylim(0.0, 0.5)
    # axes labels
    ax.set_xlabel('$x$', fontsize=15)
    ax.set_ylabel('$x H(x, x, t)$', fontsize=18)
    ax.legend()
    ax.text(0.001, 0.405, "t = 0")
    ax.text(0.001, 0.12, "t = -0.3 GeV^2")
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def plotBCSA(ff, fits=[], path=None, fmt='png'):
    """Makes plot of COMPASS BCSA for FormFactors model ff."""

    title = '' # 'Fig 15'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    fig.subplots_adjust(bottom=0.65)
    pt = Data.DummyPoint()
    pt.exptype = 'fixed target'
    pt.in1energy = 160.
    pt.xB = 0.05
    pt.t = -0.2
    pt.Q2 = 2.
    ax = fig.add_subplot(1,1,1)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    phi = np.arange(0., np.pi, 0.2)
    linestyles = ['g--', 'b-', 'r-.']
    labels = ['HERMES+CLAS', 'HERMES+CLAS+HALLA', '+HALLA(phi)']
    #labels = ['GLO1 (DM)', 'GLO1 (KK)', '']
    pn = 0
    for (approach, pars) in fits:
        pt.prepare(approach)
        line = approach.BCSA(pt, np.pi - phi, pars)
        ax.plot(phi, line, linestyles[pn], linewidth=2, label=labels[pn]) 
        pn += 1
    #ax.set_ylim(0.0, 0.5)
    # axes labels
    ax.set_xlabel('$\\phi$')
    ax.set_ylabel('BCSA')
    ax.legend()
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig
# FIXME: doesn't work!


def plotHALLAphi(pars):
    """Makes plot of HALL-A cross-section data with fit line defined by pars"""
    datafiles = [
             "data/ep2epgamma-XLU-HALLA-06-Q2_15_t_017.dat",
             #"data/ep2epgamma-XLU-HALLA-06-Q2_15_t_023.dat",
             #"data/ep2epgamma-XLU-HALLA-06-Q2_15_t_028.dat",
             #"data/ep2epgamma-XLU-HALLA-06-Q2_15_t_033.dat",
             #"data/ep2epgamma-XLU-HALLA-06-Q2_19_t_017.dat",
             #"data/ep2epgamma-XLU-HALLA-06-Q2_19_t_023.dat",
             "data/ep2epgamma-XLU-HALLA-06-Q2_19_t_028.dat",
             #"data/ep2epgamma-XLU-HALLA-06-Q2_19_t_033.dat",
             #"data/ep2epgamma-XLU-HALLA-06-Q2_23_t_017.dat",
             #"data/ep2epgamma-XLU-HALLA-06-Q2_23_t_023.dat",
             #"data/ep2epgamma-XLU-HALLA-06-Q2_23_t_028.dat",
             "data/ep2epgamma-XLU-HALLA-06-Q2_23_t_033.dat",
             "data/ep2epgamma-XUU-HALLA-06-Q2_23_t_017.dat",
             "data/ep2epgamma-XUU-HALLA-06-Q2_23_t_023.dat",
             #"data/ep2epgamma-XUU-HALLA-06-Q2_23_t_028.dat",
             "data/ep2epgamma-XUU-HALLA-06-Q2_23_t_033.dat"
             ]
    fig = plt.figure()
    panel = 1
    for file in datafiles: 
        datapoints = []
        loadfile(file, datapoints)
        ax = fig.add_subplot(2,3,panel)
        plotpoints(ax, datapoints, pars, 'phi', ['Q2', 't'])
        #plt.ylim(-0.05, 0.24)
        panel += 1
    fig.show()



def minustest(pars):
    """Tabulates A0-r*A1 combination"""
    A0points = []
    A1points = []
    loadfile("data/ep2epgamma-BCA-HERMES-08-cos0_b.dat", A0points)
    loadfile("data/ep2epgamma-BCA-HERMES-08-cos1_b.dat", A1points)
    print "  A0     A1    |    r     |   m_exp  |  m_the"
    print "---------------------------------------------"
    for n in range(len(A0points)):
        p0 = A0points[n]
        p1 = A1points[n]
        print "% 1.3f % 1.3f  |  % 1.3f  |  % 1.3f  | % 1.3f" % (p0.val, p1.val, p1.r, 
                p0.val - p1.r*p1.val, BCA0minusr1(p0, pars))

def pminustest(pars):
    """Tabulates A0-r*A1 combination"""
    A0points = []
    A1points = []
    loadfile("data/ep2epgamma-BCA-HERMES-08-cos0_a.dat", A0points)
    loadfile("data/ep2epgamma-BCA-HERMES-08-cos1_a.dat", A1points)
    print "  A0     A1    |    r     |   m_exp  |  m_the"
    print "---------------------------------------------"
    for n in range(len(A0points)):
        p0 = A0points[n]
        p1 = A1points[n]
        print "% 1.3f % 1.3f  |  % 1.3f  |  % 1.3f  | % 1.3f" % (p0.val, p1.val, p1.r, 
                p0.val - p1.r*p1.val, BCA0minusr1(p0, pars))

