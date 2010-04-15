""" 
Plotting functions written for ad-hoc plotting
of some specific datasets and CFFs.

"""

import sys, os, math
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


def subplot(ax, sets, lines=[], band=[], xaxis=None, kinlabels=[], plotlines=True):
    """Plot datapoints together with fit/theory line(s).

    ax -- subplot of matplotlib's figure i.e. ax = figure.add_subplot(..)
    sets --  list of `DataSet` instances to be plotted.
             First set defines abscissa. 
    lines -- list of 'theories' describing curves for plotting.
    band  -- neural network 'theory' or a list of 'theories' defining 
             band by their mean and standard deviation
    xaxis -- abscissa variable; if None, last of sets.xaxes is taken
    kinlabels -- list of constant kinematic variables whose values will
                 be put on plot as annotation
    TODO: legend

    """
    # first, fix the input if needed
    if not isinstance(kinlabels, list): kinlabels = [kinlabels]
    if not isinstance(lines, list): lines = [lines]
    if not xaxis: xaxis = sets[-1].xaxes[-1]
    # [1] Data sets (or fits with errorbars)
    setshapes = ['o', 's']  # first circles, then squares ...
    setcolors = ['blue', 'black']  # circles are blue, squares are black, ...
    setn = 0
    for set in sets:
        xval = []; yval = []; yerr = []
        for pt in set:
            xval.append(getattr(pt, xaxis)) 
            yval.append(pt.val)
            yerr.append(pt.err)
        ax.errorbar(xval, yval, yerr, linestyle='None', elinewidth=setn+1, 
                marker=setshapes[setn], color=setcolors[setn])
        setn += 1
    # [2] Theory lines
    lineshapes = ['s', '^', 'd', 'h']  # first squares, then triangles, diamonds, hexagons
    linecolors = ['red', 'green', 'brown', 'purple']  # squares are red, etc.
    linestyles = ['-', '--', '-.', ':']  # solid, dashed, dot-dashed, dotted
    linen = 0
    for theory in lines:
        # take abscissae from the first set
        line = [theory.predict(pt) for pt in sets[-1]]
        if plotlines:
            # join the dots (xval belongs to last set)
            ax.plot(xval, line, color=linecolors[linen], 
                    linestyle=linestyles[linen], linewidth=2)
        else:
            # put symbols on dots
            ax.plot(xval, line, lineshapes[linen], markersize=5,
                    markerfacecolor=linecolors[linen], markeredgecolor='black')
        linen += 1
    # [3] Theory band
    up = []
    down = []
    for pt in sets[-1]:
        if isinstance(band, list):
            # we have list of theories
            res = np.array([theory.predict(pt) for theory in band])
        else:
            # we have neural network.
            try:
                res = band.predict(pt, parameters={'nnet':'ALL'})
            except ValueError:  # shape mismatch so we have to go one by one
                res = np.array([band.predict(pt, parameters={'nnet':nnet}) for 
                    nnet in range(len(band.model.nets))])
        mean = res.mean()
        std = res.std()
        up.append(mean + std/2.)
        down.append(mean - std/2.)
    up = np.array(up)
    down = np.array(down)
    # xval belongs to last set after loop above
    x = plt.concatenate( (xval, xval[::-1]) )
    y = plt.concatenate( (up, down[::-1]) )
    ax.fill(x, y, facecolor='g', alpha=0.5)
    # [4] Axes 
    ax.set_xlabel(toTeX[xaxis], fontsize=15)
    ax.set_ylabel(toTeX[sets[-1][0].yaxis], fontsize=18)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    # [5] Annotations
    # constant kinematic variables positioning
    labx = min(0, min(xval)) + (max(xval) - min(0, min(xval))) * 0.5
    laby = min(0, min(yval)) + (max(yval) - min(0, min(yval))) * 0.05
    labtxt = ""
    for lab in kinlabels:
        try:
            labtxt += toTeX[lab] + ' = ' + str(getattr(sets[0],lab)) + ', '
        except AttributeError:
            # If dataset doesn't have it, all points should have it 
            labtxt += toTeX[lab] + ' = ' + str(getattr(sets[0][0],lab)) + ', '
    ax.text(labx, laby, labtxt[:-2])
    # ax.frame.set_linewidth(5) # ???
    return

def HERMESBCA(data, lines=[], band=[], path=None, fmt='png'):
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
        subplot(ax, [data[id][x*6+18:x*6+6+18]], lines, band, xaxes[x], [])
        #apply(ax.set_ylim, ylims[(-0.30, 0.05)])
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def HERMESBSA(data, lines=[], band=[], path=None, fmt='png'):
    """Plot HERMES BSA."""
    id = 5
    title = 'HERMES BSA'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    xaxes = ['tm', 'xB', 'Q2']
    # we have 3x6 points
    for x in range(3):
        ax = fig.add_subplot(3, 1, x+1)
        ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))  # tickmarks
        subplot(ax, [data[5][x*6:x*6+6]], lines, band, xaxes[x], [])
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def HERMES09(data, lines=[], band=[], path=None, fmt='png'):
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
            subplot(ax, [data[id][x*6+shift:x*6+6+shift]], lines, band, xaxes[x], [])
            apply(ax.set_ylim, ylims[y])
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def CLAS(data, lines=[], band=[], path=None, fmt='png'):
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
        subplot(ax, [panelset], lines, band, 'tm', ['Q2', 'xB'])
        plt.xlim(0.0, 0.6)
        plt.ylim(0.0, 0.4)
        ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def HALLA(data, lines=[], band=[], path=None, fmt='png'):
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
        subplot(ax, [subsets[id]], lines, band, 'phi', ['Q2', 't'])
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(120))
        panel += 1
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def COMPASS(lines=[], band=[], path=None, fmt='png'):
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
    for approach in lines:
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
    for approach in lines:
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
    for approach in lines:
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

def nnH(ff, path=None, fmt='png'):
    """Make plot of Im(cffH)."""
    title = 'GPD H by neural nets (fit to HERMES BSA only!)'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    #fig.subplots_adjust(bottom=0.35)
    old = ff.parameters['nnet']
    ff.parameters['nnet'] = 'ALL'
    pt = Data.DummyPoint()
    ax = fig.add_subplot(2,1,1)
    ax.set_xscale('log')  # x-axis to be logarithmic
    xvals = np.power(10., np.arange(-3.5, 0, 0.01)) 
    pt.t = 0.0
    up = []
    down = []
    for x in xvals:
        pt.xB = x
        nnres = x * ff.ImH(pt) / np.pi
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
        nnres = x * ff.ImH(pt) / np.pi
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
        nnres = x * ff.ImH(pt) / np.pi
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
        nnres = x * ff.ImH(pt) / np.pi
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
    ff.parameters['nnet'] = old
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig
