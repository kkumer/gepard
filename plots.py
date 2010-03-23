""" 
Plotting functions written for ad-hoc plotting
of some specific datasets and CFFs. Note that actual dataset
plotting algorithm is utils.subplot() and here only choice
of points and axes decoration etc. is done

plotHALLA -- plots HALL-A data

"""


import os
import numpy as np

import matplotlib
if os.sys.platform == 'win32':
    matplotlib.use('WxAgg')
else: #linux
    matplotlib.use('TkAgg')
import pylab as plt


import Data
import utils
from constants import toTeX, Mp2, Mp


def HERMESBCA(data, path=None, fmt='png'):
    """Plot HERMES BSA and neural net result."""
    id = 32 
    title = 'HERMES BCA'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    xaxes = ['tm', 'xB', 'Q2']
    # we have 3x6 points
    for x in range(3):
        ax = fig.add_subplot(3, 1, x+1)
        ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))  # tickmarks
        utils.subplot(ax, [data[id][x*6+18:x*6+6+18]], xaxes[x], [], [])
        #apply(ax.set_ylim, ylims[(-0.30, 0.05)])
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def HERMES09(data, fits=[], path=None, fmt='png'):
    """Plot HERMES 0909.3587 BCA and BSA data with fit lines."""

    ids = [2, 4, 5]
    title = '' #'HERMES-08'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    xaxes = ['tm', 'xB', 'Q2']
    ylims = [(-0.05, 0.24), (-0.15, 0.05), (-0.30, 0.05)]
    # we have 3x18=54 points to be separated in nine panels six points each:
    for y, id, shift in zip(range(3), [32, 32, 5], [18, 0, 0]):
        for x in range(3):
            panel = 3*y + x + 1  # 1, 2, ..., 9
            ax = fig.add_subplot(3,3,panel)
            ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))  # tickmarks
            utils.subplot(ax, [data[id][x*6+shift:x*6+6+shift]], xaxes[x], [], fits)
            apply(ax.set_ylim, ylims[y])
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig


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
    xaxes = ['tm', 'xB', 'Q2']
    ylims = [(-0.05, 0.24), (-0.15, 0.05), (-0.30, 0.05)]
    # we have 3x18=54 points to be separated in nine panels six points each:
    for y in range(3):
        for x in range(3):
            panel = 3*y + x + 1  # 1, 2, ..., 9
            ax = fig.add_subplot(3,3,panel)
            ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))  # tickmarks
            utils.subplot(ax, [data[ids[y]][x*6:x*6+6]], xaxes[x], [], fits)
            apply(ax.set_ylim, ylims[y])
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def plotHERMESBSA(data, nn, path=None, fmt='png'):
    """Plot HERMES BSA and neural net result."""
    title = 'HERMES vs NeuralNets'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    xaxes = ['tm', 'xB', 'Q2']
    # we have 3x6 points
    for x in range(3):
        ax = fig.add_subplot(3, 1, x+1)
        ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))  # tickmarks
        utils.subplot(ax, [data[5][x*6:x*6+6]], xaxes[x], [], [])
        #apply(ax.set_ylim, ylims[(-0.30, 0.05)])
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
        utils.subplot(ax, [panelset], 'tm', ['Q2', 'xB'], fits)
        plt.xlim(0.0, 0.6)
        plt.ylim(0.0, 0.4)
        ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
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
        utils.subplot(ax, [data[id]], 'phi', ['Q2', 't'], fits)
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(120))
        panel += 1
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def plotCOMPASS(ff, fits=[], path=None, fmt='png'):
    """Plot COMPASS BCS asymmetry and difference and summ of xs for FormFactors model ff.
    
    FIXME: kinematic completion explicit here
    """

    title = 'COMPASS BCSA (asymmetry), BCSD (difference) and BCSS (sum) of XS'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    #fig.subplots_adjust(bottom=0.35)
    # Asymmetry panel
    pt = Data.DummyPoint()
    pt.exptype = 'fixed target'
    pt.in1energy = 160.
    pt.s = 2 * Mp * pt.in1energy + Mp2
    pt.xB = 0.05
    pt.t = -0.2
    pt.Q2 = 2.
    ax = fig.add_subplot(2,2,1)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    phi = np.arange(0., np.pi, 0.2)
    utils.fill_kinematics(pt)
    linestyles = ['g--', 'b-', 'r-.']
    labels = ['HERMES+CLAS', 'HERMES+CLAS+HALLA', '+HALLA(phi)']
    #labels = ['GLO1 (DM)', 'GLO1 (KK)', '']
    pn = 0
    for (approach, pars) in fits:
        approach.to_conventions(pt)
        approach.prepare(pt)
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
    pt.s = 2 * Mp * pt.in1energy + Mp2
    pt.xB = 0.05
    pt.t = -0.2
    pt.Q2 = 2.
    ax = fig.add_subplot(2,2,2)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    phi = np.arange(0., np.pi, 0.2)
    utils.fill_kinematics(pt)
    linestyles = ['g--', 'b-', 'r-.']
    labels = ['HERMES+CLAS', 'HERMES+CLAS+HALLA', '+HALLA(phi)']
    #labels = ['GLO1 (DM)', 'GLO1 (KK)', '']
    pn = 0
    for (approach, pars) in fits:
        approach.to_conventions(pt)
        approach.prepare(pt)
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
    pt.s = 2 * Mp * pt.in1energy + Mp2
    pt.xB = 0.05
    pt.t = -0.2
    pt.Q2 = 2.
    utils.fill_kinematics(pt)
    ax = fig.add_subplot(2,2,3)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    phi = np.arange(0., np.pi, 0.2)
    linestyles = ['g--', 'b-', 'r-.']
    labels = ['HERMES+CLAS', 'HERMES+CLAS+HALLA', '+HALLA(phi)']
    #labels = ['GLO1 (DM)', 'GLO1 (KK)', '']
    pn = 0
    for (approach, pars) in fits:
        approach.to_conventions(pt)
        approach.prepare(pt)
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

def plotCOMPASSt(ff, fits=[], path=None, fmt='png'):
    """Makes plot of COMPASS BCS asymmetry in terms of t for various values of xB"""

    title = 'Prediction for COMPASS BCSA'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    kinpoints = [(0.007, 1.5), (0.014, 2.5), (0.024, 3.7), 
                           (0.045, 3.0), (0.1, 3.0), (0.2, 4.4)]
    panel = 1
    for xB, Q2 in kinpoints:
        ax = fig.add_subplot(2,3,panel)
        ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
        linestyles = ['g--', 'b-', 'r-.']
        labels = ['HERMES+CLAS', 'HERMES+CLAS+HALLA', '+HALLA(phi)']
        pn = 0
        for (approach, pars) in fits:
            line = []
            tmvals = []
            for tm in np.arange(0.05, 0.6, 0.02):
                tmvals.append(tm)
                pt = Data.DummyPoint()
                pt.exptype = 'fixed target'
                pt.in1energy = 160.
                pt.s = 2 * Mp * pt.in1energy + Mp2
                pt.xB = xB
                pt.Q2 = Q2
                pt.phi = 0.
                pt.units = {'phi' : 'radian'}
                pt.frame = 'Trento'
                pt.tm = tm
                utils.fill_kinematics(pt)
                approach.to_conventions(pt)
                approach.prepare(pt)
                line.append(approach.BCSA(pt, 0.8, pars))
            ax.plot(tmvals, line, linestyles[pn], linewidth=2, 
                    label=labels[pn]) 
            pn += 1
        ax.set_xlim(0.0, 0.8)
        ax.set_ylim(-0.1, 0.4)
        # axes labels
        ax.set_xlabel('$-t$')
        ax.text(0.03, 0.35, "%s = %s" % (toTeX['xB'], xB))
        ax.text(0.03, 0.3, "%s = %s" % (toTeX['Q2'], Q2))
        panel += 1
    #ax.set_ylabel('BCS Asymmetry')
    #ax.legend()
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
    pt.s = 2 * Mp * pt.in1energy + Mp2
    pt.xB = 0.05
    pt.t = -0.2
    pt.Q2 = 2.
    ax = fig.add_subplot(1,2,2)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    phi = np.arange(0., np.pi, 0.2)
    utils.fill_kinematics(pt)
    linestyles = ['g--', 'b-', 'r-.']
    labels = ['HERMES+CLAS', 'HERMES+CLAS+HALLA', '+HALLA(phi)']
    #labels = ['GLO1 (DM)', 'GLO1 (KK)', '']
    pn = 0
    for (approach, pars) in fits:
        approach.to_conventions(pt)
        approach.prepare(pt)
        line = approach.BCSA(pt, 1, pars, {'phi':np.pi - phi})
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
    pt.s = 2 * Mp * pt.in1energy + Mp2
    pt.xB = 0.05
    pt.t = -0.2
    pt.Q2 = 2.
    ax = fig.add_subplot(1,1,1)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    phi = np.arange(0., np.pi, 0.2)
    utils.fill_kinematics(pt)
    linestyles = ['g--', 'b-', 'r-.']
    labels = ['HERMES+CLAS', 'HERMES+CLAS+HALLA', '+HALLA(phi)']
    #labels = ['GLO1 (DM)', 'GLO1 (KK)', '']
    pn = 0
    for (approach, pars) in fits:
        approach.to_conventions(pt)
        approach.prepare(pt)
        line = approach.BCSA(pt, 1, pars, {'phi':np.pi - phi})
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


def plotnnBSA(data, nnapproach, path=None, fmt='png'):
    """Plot HERMES BSA and neural net result."""
    title = 'HERMES (blue) vs NeuralNets (black)'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    xaxes = ['tm', 'xB', 'Q2']
    # we have 3x6 points
    for x in range(3):
        datapoints = data[5][x*6:x*6+6]
        ax = fig.add_subplot(3, 1, x+1)
        ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))  # tickmarks
        nnpoints = []
        for pt in datapoints:
            nnpt = Data.DummyPoint(pt.__dict__.copy())
            nnres = nnapproach.ALUIsin1(pt, {})
            nnpt.val = nnres.mean()
            nnpt.err = nnres.std()
            nnpoints.append(nnpt)
        utils.subplot(ax, [datapoints, nnpoints], xaxes[x], [], [])
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def plotnnH(ff, path=None, fmt='png'):
    """Makes plot of Im(cffH)."""

    title = 'GPD H by neural nets (fit to HERMES BSA only!)'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    #fig.subplots_adjust(bottom=0.35)
    pt = Data.DummyPoint()
    ax = fig.add_subplot(2,1,1)
    ax.set_xscale('log')  # x-axis to be logarithmic
    xvals = np.power(10., np.arange(-3.5, 0, 0.01)) 
    pt.t = 0.0
    up = []
    down = []
    for x in xvals:
        pt.xB = x
        nnres = x * ff.ImH(pt, {}) / np.pi
        mean = nnres.mean()
        std = nnres.std()
        up.append(mean + std/2.)
        down.append(mean - std/2.)
    up = np.array(up)
    down = np.array(down)
    x = plt.concatenate( (xvals, xvals[::-1]) )
    y = plt.concatenate( (up, down[::-1]) )
    ax.fill(x, y, facecolor='g', alpha=0.5)
    ####
    pt.t = -0.3
    up = []
    down = []
    for x in xvals:
        pt.xB = x
        nnres = x * ff.ImH(pt, {}) / np.pi
        mean = nnres.mean()
        std = nnres.std()
        up.append(mean + std/2.)
        down.append(mean - std/2.)
    up = np.array(up)
    down = np.array(down)
    x = plt.concatenate( (xvals, xvals[::-1]) )
    y = plt.concatenate( (up, down[::-1]) )
    ax.fill(x, y, facecolor='r', alpha=0.5)
    #ax.set_ylim(0.0, 0.5)
    ax.set_xlim(0.001, 1.0)
    #plt.ylim(0.0, 0.5)
    # axes labels
    ax.set_xlabel('$x$', fontsize=15)
    ax.set_ylabel('$x H(x, x, t)$', fontsize=18)
    #ax.legend()
    ax.text(0.01, 0.405, "t = 0")
    ax.text(0.01, -0.1, "t = -0.3 GeV^2")
    ####  --- SECOND PANEL ---
    ax = fig.add_subplot(2,1,2)
    xvals = np.linspace(0.05, 0.2, 20)
    pt.t = 0.0
    up = []
    down = []
    for x in xvals:
        pt.xB = x
        nnres = x * ff.ImH(pt, {}) / np.pi
        mean = nnres.mean()
        std = nnres.std()
        up.append(mean + std/2.)
        down.append(mean - std/2.)
    up = np.array(up)
    down = np.array(down)
    x = plt.concatenate( (xvals, xvals[::-1]) )
    y = plt.concatenate( (up, down[::-1]) )
    ax.fill(x, y, facecolor='g', alpha=0.5)
    ####
    pt.t = -0.3
    up = []
    down = []
    for x in xvals:
        pt.xB = x
        nnres = x * ff.ImH(pt, {}) / np.pi
        mean = nnres.mean()
        std = nnres.std()
        up.append(mean + std/2.)
        down.append(mean - std/2.)
    up = np.array(up)
    down = np.array(down)
    x = plt.concatenate( (xvals, xvals[::-1]) )
    y = plt.concatenate( (up, down[::-1]) )
    ax.fill(x, y, facecolor='r', alpha=0.5)
    #ax.set_ylim(0.0, 0.5)
    #ax.set_xlim(0.0005, 1.0)
    #plt.ylim(0.0, 0.5)
    # axes labels
    ax.set_xlabel('$x$', fontsize=15)
    ax.set_ylabel('$x H(x, x, t)$', fontsize=18)
    #ax.legend()
    ax.text(0.1, 0.405, "t = 0")
    ax.text(0.1, 0.12, "t = -0.3 GeV^2")
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig
