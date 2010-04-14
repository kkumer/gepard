""" 
Plotting functions written for ad-hoc plotting
of some specific datasets and CFFs. Note that actual dataset
plotting algorithm is utils.subplot() and here only choice
of points and axes decoration etc. is done

plotHALLA -- plots HALL-A data

"""


import os, math
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
    """Plot HERMES BCA."""
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

def CLAS(data, fits=[], path=None, fmt='png'):
    """Makes plot of CLAS BSA data with fit lines"""

    #datafile = "data/ep2epgamma-ALU-CLAS_KK-07.dat" # id = 25
    dataset = data[25]
    title = 'CLAS-07_KK'
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

def HALLA(data, fits=[], path=None, fmt='png'):
    """Makes plot of HALL-A data with fit lines"""

    ids = [9, 14, 20, 21, 23, 24]
    subsets = {}
    subsets[9] = utils.select(data[33], criteria=['Q2 == 1.5', 't == -0.17'])
    subsets[14] = utils.select(data[33], criteria=['Q2 == 1.9', 't == -0.23'])
    subsets[20] = utils.select(data[33], criteria=['Q2 == 2.3', 't == -0.33'])
    subsets[21] = utils.select(data[34], criteria=['t == -0.17'])
    subsets[23] = utils.select(data[34], criteria=['t == -0.28'])
    subsets[24] = utils.select(data[34], criteria=['t == -0.33'])
    title = '' # 'HALLA-06'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    panel = 1
    for id in ids:
        ax = fig.add_subplot(2,3,panel)
        utils.subplot(ax, [subsets[id]], 'phi', ['Q2', 't'], fits)
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(120))
        panel += 1
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def COMPASS(fits=[], path=None, fmt='png'):
    """Plot COMPASS BCS asymmetry and difference and summ of xs for FormFactors model ff.
    
    FIXME: kinematic completion, charge etc. must be explicit here
    """

    title = 'COMPASS BCSA (asymmetry), BCSD (difference) and BCSS (sum) of XS'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    #fig.subplots_adjust(bottom=0.35)
    # Asymmetry panel
    pt = Data.DummyPoint()
    pt.exptype = 'fixed target'
    pt.in1particle = 'e'
    pt.in1charge = 1
    pt.in1energy = 160.
    pt.in1polarizationvector = 'L'
    pt.in1polarization = 0.8
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
    for approach in fits:
        approach.__class__.to_conventions(pt)
        approach.__class__.prepare(pt)
        line = approach.BCSA(pt, vars={'phi':np.pi - phi})
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
    pt.in1particle = 'e'
    pt.in1charge = 1
    pt.in1energy = 160.
    pt.in1polarizationvector = 'L'
    pt.in1polarization = 0.8
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
    for approach in fits:
        approach.__class__.to_conventions(pt)
        approach.__class__.prepare(pt)
        # nb converted to pb:
        line = 1000. * approach.BCSD(pt, vars={'phi':np.pi - phi})
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
    pt.in1particle = 'e'
    pt.in1charge = 1
    pt.in1energy = 160.
    pt.in1polarizationvector = 'L'
    pt.in1polarization = 0.8
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
    for approach in fits:
        approach.__class__.to_conventions(pt)
        approach.__class__.prepare(pt)
        # nb converted to pb:
        line = 1000. * approach.BCSS(pt, vars={'phi':np.pi - phi})
        ax.plot(phi, line, linestyles[pn], linewidth=2, label=labels[pn]) 
        pn += 1
    #ax.set_ylim(0.0, 0.5)
    # axes labels
    ax.set_xlabel('$\\phi$')
    ax.set_ylabel('BCS sum  [pb]')
    ax.legend()
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def COMPASSt(fits=[], path=None, fmt='png'):
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
        for approach in fits:
            line = []
            tmvals = []
            for tm in np.arange(0.05, 0.6, 0.02):
                tmvals.append(tm)
                pt = Data.DummyPoint()
                pt.exptype = 'fixed target'
                pt.in1particle = 'e'
                pt.in1charge = 1
                pt.in1energy = 160.
                pt.in1polarizationvector = 'L'
                pt.in1polarization = 0.8
                pt.s = 2 * Mp * pt.in1energy + Mp2
                pt.xB = xB
                pt.Q2 = Q2
                pt.phi = 0.
                pt.units = {'phi' : 'radian'}
                pt.frame = 'Trento'
                pt.tm = tm
                utils.fill_kinematics(pt)
                approach.__class__.to_conventions(pt)
                approach.__class__.prepare(pt)
                line.append(approach.BCSA(pt))
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

def EIC(fits=[], path=None, fmt='png'):
    """Plot EIC cross section.
    
    FIXME: kinematic completion, charge etc. must be explicit here
    """

    title = 'EIC cross section'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    #fig.subplots_adjust(bottom=0.35)
    # Asymmetry panel
    pt = Data.DummyPoint()
    pt.exptype = 'collider'
    pt.in1particle = 'e'
    pt.in1charge = 1
    pt.in1energy = 5.
    pt.in1polarizationvector = 'L'
    pt.in1polarization = 1.0
    pt.in2particle = 'p'
    pt.in2energy = 250.
    pt.s = 2 * pt.in1energy * (pt.in2energy + math.sqrt(
        pt.in2energy**2 - Mp2)) + Mp2
    pt.xB = 0.001
    pt.t = -0.1
    pt.Q2 = 2.
    pt.units = {'phi' : 'radian'}
    pt.frame = 'Trento'
    ax = fig.add_subplot(1,1,1)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    phi = np.arange(0., 2*np.pi, 0.2)
    utils.fill_kinematics(pt)
    linestyles = ['g--', 'b-']
    labels = ['polarized model1', 'polarized model2']
    pn = 0
    for approach in fits:
        approach.__class__.to_conventions(pt)
        approach.__class__.prepare(pt)
        line = approach.Xunp(pt, vars={'phi':phi})
        ax.plot(phi, line, linestyles[pn], linewidth=1, label=labels[pn]) 
        pn += 1
    linestyles = ['r-.', 'p:']
    labels = ['unpolarized model1', 'unpolarized model2']
    pn = 0
    for approach in fits:
        approach.__class__.to_conventions(pt)
        approach.__class__.prepare(pt)
        lineBSS = approach.BSS(pt, vars={'phi':phi})
        ax.plot(phi, lineBSS, linestyles[pn], linewidth=2, label=labels[pn]) 
        pn += 1
    #ax.set_ylim(0.0, 0.5)
    # axes labels
    ax.set_xlabel('$\\phi$')
    ax.set_ylabel('Cross section')
    ax.legend()
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def H(theories=[], path=None, fmt='png'):
    """Makes plot of Im(cffH).
    
    TODO: This should be a method of Model.
    
    """
    title = '' # 'Fig 15'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    fig.subplots_adjust(bottom=0.45)
    # Left panel
    pt = Data.DummyPoint()
    ax = fig.add_subplot(1,1,1)
    ax.set_xscale('log')  # x-axis to be logarithmic
    xval = np.power(10., np.arange(-3.5, 0, 0.01)) 
    colors = ['red', 'green', 'brown', 'purple']
    styles = ['-', '--', '-.', ':']
    labels = ['model 1', 'model 2', 'model 3', 'model 4']
    pn = 0
    for t in theories:
        pt.t = 0.0
        # kludge alert!
        line1 = xval * t.model.ImH(pt, xval) / np.pi
        pt.t = -0.3
        line2 = xval * t.model.ImH(pt, xval) / np.pi
        ax.plot(xval, line1, color=colors[pn], linestyle=styles[pn], linewidth=2)
        ax.plot(xval, line2, color=colors[pn], linestyle=styles[pn], linewidth=4, label=labels[pn]) 
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

def HBCSA(ff, fits=[], path=None, fmt='png'):
    """Makes plot of Im(cffH) and COMPASS BCSA for FormFactors model ff."""

    title = 'Fig 15'
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
    for approach in fits:
        pt.t = 0.0
        # kludge alert!
        line1 = xval * ff.ImH(pt, xval) / np.pi
        pt.t = -0.3
        line2 = xval * ff.ImH(pt, xval) / np.pi
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
    pt.in1particle = 'e'
    pt.in1charge = 1
    pt.in1energy = 160.
    pt.in1polarizationvector = 'L'
    pt.in1polarization = 0.8
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
    for approach in fits:
        approach.__class__.to_conventions(pt)
        approach.__class__.prepare(pt)
        line = approach.BCSA(pt, vars={'phi':np.pi - phi})
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

def nnBSA(data, nnapproach, path=None, fmt='png'):
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
            nnres = nnapproach.BSA(pt, {})
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

def nnH(ff, path=None, fmt='png'):
    """Make plot of Im(cffH)."""
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
