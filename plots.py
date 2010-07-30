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

# load experimental data
data = utils.loaddata('data/ep2epgamma')  # dictionary {1 : DataSet instance, ...}
data.update(utils.loaddata('data/gammastarp2gammap')) 


def subplot(ax, sets, lines=[], band=[], xaxis=None, kinlabels=[], plotlines=True):
    """Plot datapoints together with fit/theory line(s).

    ax -- subplot of matplotlib's figure i.e. ax = figure.add_subplot(..)
    sets --  list of `DataSet` instances to be plotted.
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
    setshapes = ['o', 's', '^', 'd']  # first circles, then squares ...
    setcolors = ['blue', 'black', 'purple', 'green']  # circles are blue, squares are black, ...
    setn = 0
    xsets = []
    for set in sets:
        xval = []; yval = []; yerr = []
        for pt in set:
            xval.append(getattr(pt, xaxis)) 
            yval.append(pt.val)
            yerr.append(pt.err)
        xsets.append(xval)
        ax.errorbar(xval, yval, yerr, linestyle='None', elinewidth=1, 
                marker=setshapes[setn], color=setcolors[setn])
        setn += 1
    # [2] Theory lines
    lineshapes = ['s', '^', 'd', 'h']  # first squares, then triangles, diamonds, hexagons
    linecolors = ['red', 'green', 'blue', 'purple']  # squares are red, etc.
    linestyles = ['-', '--', '-.', ':']  # solid, dashed, dot-dashed, dotted
    linen = 0
    for theory in lines:
        # Each theory has to predict each set:
        setn = 0
        for set in sets:
            line = [theory.predict(pt) for pt in set]
            if plotlines:
                # join the dots (xval belongs to last set)
                ax.plot(xsets[setn], line, color=linecolors[linen], 
                        linestyle=linestyles[linen], linewidth=2)
            else:
                # put symbols on dots
                ax.plot(xsets[setn], line, lineshapes[linen], markersize=5,
                        markerfacecolor=linecolors[linen], markeredgecolor='black')
            setn += 1
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
                #FIXME: if len(band.model.nets)==5 exception will NOT occur!
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
    labx = min(0, min(xval)) + (max(xval) - min(0, min(xval))) * 0.35
    laby = min(0, min(yval)) + (max(yval) - min(0, min(yval))) * 0.02
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

def HERMESBCA(lines=[], band=[], path=None, fmt='png'):
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

def HERMESBSA(lines=[], band=[], path=None, fmt='png'):
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

def HERMES09(lines=[], band=[], path=None, fmt='png'):
    """Plot HERMES 0909.3587 BCA and BSA data with fit lines."""

    #ids = [2, 4, 5]
    title = 'HERMES 09'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    xaxes = ['tm', 'xB', 'Q2']
    ylims = [(-0.05, 0.3), (-0.15, 0.15), (-0.45, 0.05)]
    # we have 3x18=54 points to be separated in nine panels six points each:
    for y, id, shift in zip(range(3), [32, 32, 5], [18, 0, 0]):
        for x in range(3):
            panel = 3*y + x + 1  # 1, 2, ..., 9
            ax = fig.add_subplot(3,3,panel)
            ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))  # tickmarks
            subplot(ax, [data[id][x*6+shift:x*6+6+shift]], lines, band, xaxes[x], [])
            apply(ax.set_ylim, ylims[y])
            if (panel % 3) != 1:
                # Leave labels only on leftmost panels
                ax.set_ylabel('')
            else:
                ylabels = ['$BCA\\; \\cos \\phi$', '$BCA\\; \\cos 0\\phi$', '$BSA\\; \\sin\\phi$']
                ax.set_ylabel(ylabels[(panel-1)/3], fontsize=18)

            if panel < 7:
                # Leave labels only on lowest panels
                ax.set_xlabel('')
            else:
                xlabels = ['$-t\\; [{\\rm GeV}^2]$', '$x_B$', '$Q^2\\; [{\\rm GeV}^2]$']
                ax.set_xlabel(xlabels[panel-7], fontsize=18)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def CLAS(lines=[], band=[], path=None, fmt='png'):
    """Makes plot of CLAS BSA data with fit lines"""

    #datafile = "data/ep2epgamma-ALU-CLAS_KK-07.dat" # id = 25
    dataset = data[25]
    title = 'CLAS 07'
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
        subplot(ax, [panelset], lines, band, 'tm', ['Q2', 'xB'],
                plotlines=not(panel==1 or panel==3))
        plt.xlim(0.0, 0.6)
        plt.ylim(0.0, 0.4)
        ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
        if (panel % 3) != 1:
            # Leave labels only on leftmost panels
            ax.set_ylabel('')
        if panel < 7:
            ax.set_xlabel('')
        else:
            ax.set_xlabel('$-t\\; [{\\rm GeV}^2]$', fontsize=18)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def HALLA(lines=[], band=[], path=None, fmt='png'):
    """Makes plot of HALL-A data with fit lines"""

    subsets = {}
    subsets[1] = utils.select(data[26], criteria=['Q2 == 1.9'])
    subsets[2] = utils.select(data[35], criteria=['FTn == 1'])
    subsets[3] = utils.select(data[35], criteria=['FTn == 0'])
    subsets[4] = data[30]
    title = 'HALL-A'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    for panel in range(1,5):
        ax = fig.add_subplot(2,3,panel)
        subplot(ax, [subsets[panel]], lines, band, 't')
        if panel == 2:
            ax.set_ylabel(toTeX['ReCI'], fontsize=18)
        elif panel == 3:
            ax.set_ylabel(toTeX['ReCpDCI'], fontsize=18)
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def H1ZEUS(lines=[], band=[], path=None, fmt='png'):
    """Makes plot of H1 DVCS data with fit lines"""

    subsets = {}
    subsets[1] = [utils.select(data[36], criteria=['Q2 == 8.']),
            utils.select(data[36], criteria=['Q2 == 15.5']), 
            utils.select(data[36], criteria=['Q2 == 25.'])]
    subsets[2] = [data[46]] # ZEUS t-dep
    subsets[3] = [data[48],
            utils.select(data[41], criteria=['Q2 == 8.']),
            utils.select(data[41], criteria=['Q2 == 15.5']), 
            utils.select(data[41], criteria=['Q2 == 25.'])]
    subsets[4] = [data[47]] # ZEUS Q2-dep
    xs = ['t', 't', 'W', 'Q2']
    title = 'H1 07 / ZEUS 08'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    for panel in range(1,5):
        ax = fig.add_subplot(2,2,panel)
        ax.set_yscale('log')  # y-axis to be logarithmic
        subplot(ax, subsets[panel], lines, band, xs[panel-1])
        if panel < 3:
            ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.2))
        # y labels
        if panel==1:
            ax.set_ylabel('$d\\sigma/dt\\quad [{\\rm nb/GeV}^2]$', fontsize=18)
        elif panel==3:
            ax.set_ylabel('$\\sigma\\quad [{\\rm nb}]$', fontsize=18)
        else: # panel 4
            ax.set_ylabel('')
        # x labels
        if panel==1 or panel==2:
            ax.set_xlabel('$t\\quad [{\\rm GeV}^2]$', fontsize=18)
        elif panel==3:
            ax.set_xlabel('$W\\quad [{\\rm GeV}]$', fontsize=18)
        else: # panel 4
            ax.set_xlabel('$Q^2\\quad [{\\rm GeV}^2]$', fontsize=18)
        if panel==1:
            ax.text(-0.8, 22, '${\\rm H1}$', fontsize=16)
            ax.text(-0.8, 12, '${\\rm W = 82}\\, {\\rm GeV}$', fontsize=14)
            ax.text(-0.8, 6, '$Q^2\\!= 8,\\, 15.5,\\, 25\\,{\\rm GeV}^2$', fontsize=14)
        if panel==2:
            ax.text(-0.3, 6, '${\\rm ZEUS}$', fontsize=16)
            ax.text(-0.3, 4.5, '${\\rm W = 104}\\, {\\rm GeV}$', fontsize=14)
            ax.text(-0.3, 3.3, '$Q^2\\!= 3.2\\,{\\rm GeV}^2$', fontsize=14)
        if panel==3:
            ax.text(50, 28, '${\\rm ZEUS}\\, (idem):$', fontsize=16)
            ax.text(50, 5, '${\\rm H1}\\, (idem):$', fontsize=16)
        else: # panel==4
            ax.text(30, 3, '${\\rm ZEUS}\\, (idem)$', fontsize=16)
            #ax.set_xlim(0, 80)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def HALLAphi(lines=[], band=[], path=None, fmt='png'):
    """Makes plot of HALL-A data with fit lines"""

    ids = [9, 14, 20, 21, 23, 24]
    subsets = {}
    subsets[9] = utils.select(data[33], criteria=['Q2 == 1.5', 't == -0.17'])
    subsets[14] = utils.select(data[33], criteria=['Q2 == 1.9', 't == -0.23'])
    subsets[20] = utils.select(data[33], criteria=['Q2 == 2.3', 't == -0.33'])
    subsets[21] = utils.select(data[34], criteria=['t == -0.17'])
    subsets[23] = utils.select(data[34], criteria=['t == -0.28'])
    subsets[24] = utils.select(data[34], criteria=['t == -0.33'])
    title = 'HALL A 06'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    panel = 1
    for id in ids:
        ax = fig.add_subplot(2,3,panel)
        subplot(ax, [subsets[id]], lines, band, 'phi', ['Q2', 't'])
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(2.))
        ax.set_xlabel('$\\phi\\; {\\rm [rad]}$', fontsize=18)
        if not(panel==1 or panel==4):
            # Leave labels only on leftmost panels
            ax.set_ylabel('')
        elif panel==1:
            ax.set_ylabel('$d^4\\Sigma/(dQ^2\\! dx_{\\rm B}dt d\\phi)\\quad \
                    [{\\rm nb/GeV}^4]$', fontsize=18)
        else: # panel 4
            ax.set_ylabel('$d^4\\sigma/(dQ^2\\! dx_{\\rm B}dt d\\phi)\\quad  \
                    [{\\rm nb/GeV}^4]$', fontsize=18)
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
    labels = ['', '', '']
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
    labels = ['', '', '']
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
    #labels = ['HERMES+CLAS', 'HERMES+CLAS+HALLA', '+HALLA(phi)']
    labels = ['', '', '']
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
    linestyles = ['b-.', 'p:', 'r-', 'g--'] 
    labels = [r'$\sigma^{\uparrow}$ ' + t.name for t in fits]
    pn = 0
    for approach in fits:
        approach.__class__.to_conventions(pt)
        approach.__class__.prepare(pt)
        line = approach.Xunp(pt, vars={'phi':phi})
        ax.plot(phi, line, linestyles[pn], linewidth=1, label=labels[pn]) 
        pn += 1
    linestyles = ['r-', 'g--', 'b-.', 'p:'] 
    labels = ['$(\\sigma^{\\uparrow} + \\sigma^{\\downarrow})/2$ ' + t.name for t in fits]
    pn = 0
    for approach in fits:
        approach.__class__.to_conventions(pt)
        approach.__class__.prepare(pt)
        lineBSS = approach.BSS(pt, vars={'phi':phi})
        ax.plot(phi, lineBSS, linestyles[pn], linewidth=2, label=labels[pn]) 
        pn += 1
    #ax.set_ylim(0.0, 0.5)
    # axes labels
    ax.set_xlabel('$\\phi\\quad {\\rm [rad]}$', fontsize=20)
    ax.set_ylabel('$\\sigma\\quad {\\rm [nb]}$', fontsize=20)
    #ax.legend(loc=(0.3,0))
    ax.legend()
    ax.text(0.7, 87., "EIC", fontsize=18)
    ax.text(0.7, 82., "$E_e = %d \\,{\\rm GeV}$" % pt.in1energy, fontsize=18)
    ax.text(0.7, 79., "$E_p = %d \\,{\\rm GeV}$" % pt.in2energy, fontsize=18)
    ax.text(0.7, 76., "$x_B = %s$" % pt.xB, fontsize=18)
    ax.text(0.7, 73., "$Q^2 = %s \\,{\\rm GeV}^2$" % pt.Q2, fontsize=18)
    ax.text(0.7, 71., "$t = %s \\,{\\rm GeV}^2$" % pt.t, fontsize=18)
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
    title = 'GPD H at LO' # 'Fig 15'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    fig.subplots_adjust(bottom=0.45)
    # Left panel
    pt = Data.DummyPoint()
    pt.Q2 = 2.
    ax = fig.add_subplot(1,1,1)
    ax.set_xscale('log')  # x-axis to be logarithmic
    xval = np.power(10., np.arange(-3.5, 0, 0.01)) 
    colors = ['red', 'green', 'brown', 'purple']
    styles = ['-', '--', '-.', ':']
    pn = 0
    for t in theories:
        pt.t = -0.1
        # kludge alert!
        ImHvals = []
        for xi in xval:
            if t.model.__dict__.has_key('Gepard'): t.m.g.newcall = 1
            ImHvals.append(t.model.ImH(pt, xi))
        line1 = xval * np.array(ImHvals) / np.pi
        pt.t = -0.3
        ImHvals = []
        for xi in xval:
            if t.model.__dict__.has_key('Gepard'): t.m.g.newcall = 1
            ImHvals.append(t.model.ImH(pt, xi))
        line2 = xval * ImHvals / np.pi
        ax.plot(xval, line1, color=colors[pn], linestyle=styles[pn], linewidth=2)
        ax.plot(xval, line2, color=colors[pn], linestyle=styles[pn], linewidth=4, label=t.name) 
        pn += 1
    ax.set_ylim(0.0, 0.5)
    ax.set_xlim(0.0005, 1.0)
    #plt.ylim(0.0, 0.5)
    # axes labels
    ax.set_xlabel('$x$', fontsize=15)
    ax.set_ylabel('$x H(x, x, t)$', fontsize=18)
    ax.legend()
    #ax.text(0.001, 0.405, "t = 0")
    ax.text(0.001, 0.27, "$t = -0.1\\, {\\rm GeV}^2$", fontsize=15)
    ax.text(0.001, 0.12, "$t = -0.3\\, {\\rm GeV}^2$", fontsize=15)
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
    pn = 0
    for approach in fits:
        approach.__class__.to_conventions(pt)
        approach.__class__.prepare(pt)
        line = approach.BCSA(pt, vars={'phi':np.pi - phi})
        ax.plot(phi, line, linestyles[pn], linewidth=2, label=approach.name) 
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

def _axband(ax, tm, xvals, fun, color='g', avg=True):
    """Make a band of x*fun defined by ndarray-valued fun."""
    pt = Data.DummyPoint()
    pt.t = -tm
    pt.xi = xvals
    pt.xB = 2*pt.xi/(1.+pt.xi)
    up = []
    down = []
    res = xvals * fun(pt)
    mean = res.mean(axis=0)
    std = res.std(axis=0)
    if avg:
        up = mean + std/2.
        down = mean - std/2.
    else:
        up = res.max(axis=0)
        down = res.min(axis=0)
    x = plt.concatenate( (xvals, xvals[::-1]) )
    y = plt.concatenate( (up, down[::-1]) )
    ax.fill(x, y, facecolor=color, alpha=0.5)

def CFF(t, cffs=None, path=None, fmt='png', average=True):
    """Makes plot of cffs given by neural network
    
    t     -- Neural network 'theory'
    cffs  -- List of CFFs to be plotted. Each produces two panels.
    average -- plot 1-sigma band (othewise plot minimum-maximum band)

    """
    if not cffs:
        cffs = t.model.output_layer
    old = t.model.parameters['nnet']
    t.model.parameters['nnet'] = 'ALL'
    title = t.description
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    colors = ['red', 'blue', 'brown', 'purple']
    linestyles = ['solid', 'dashed']
    tms = [0.1, 0.3]
    # Define abscissas
    logxvals = np.power(10., np.arange(-3.0, 0.01, 0.1))  # left panel
    xvals = np.linspace(0.025, 0.2, 20) # right panel
    # ordinatas for  left ...
    allylims = [(-0.3, 1.0), (-0.3, 1.0), (-0.45, 0.05)]
    # ... and right panles
    ylims = [(-0.3, 1.0), (-0.3, 1.0), (-0.45, 0.05)]
    # Plot panels
    for n in range(len(cffs)):
        cff = cffs[n]
        # all-x logarithmic
        ax = fig.add_subplot(len(cffs), 2, 2*n+1)
        ax.set_xscale('log')  # x-axis to be logarithmic
        for tm in tms :
            _axband(ax, tm, logxvals, getattr(t.model, cff), 
                    color=colors[tms.index(tm)], avg=average)
        ax.set_xlabel('$\\xi$', fontsize=15)
        ax.set_ylabel('x %s' % cff, fontsize=18)
        ax.set_xlim(0.001, 1.0)
        apply(ax.set_ylim, allylims[n])
        ax.axhspan(-0.0005, 0.0005, facecolor='g', alpha=0.6)  # horizontal bar
        ax.axvspan(0.025, 0.2, facecolor='g', alpha=0.1)  # vertical band
        ax.text(0.03, -0.27, "data region", fontsize=14)
        # measured x linear
        ax = fig.add_subplot(len(cffs), 2, 2*n+2)
        for tm in tms:
            _axband(ax, tm, xvals, getattr(t.model, cff), 
                    color=colors[tms.index(tm)], avg=average)
        ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
        ax.set_xlabel('$\\xi$', fontsize=15)
        apply(ax.set_ylim, ylims[n])
        ax.set_xlim(0.025, 0.2)
        ax.text(0.03, -0.27, "data region only", fontsize=14)
        #ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.02))  # tickmarks
    t.model.parameters['nnet'] = old
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

