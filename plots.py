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


import Data, Approach, utils
from constants import toTeX, Mp2, Mp

# load experimental data
data = utils.loaddata('data/ep2epgamma', approach=Approach.hotfixedBMK) 
data.update(utils.loaddata('data/gammastarp2gammap', approach=Approach.hotfixedBMK)) 


def subplot(ax, sets, lines=[], band=[], xaxis=None, kinlabels=[], plotlines=True):
    """Plot datapoints together with fit/theory line(s).

    ax -- subplot of matplotlib's figure i.e. ax = figure.add_subplot(..)
    sets --  list of `DataSet` instances to be plotted.
    lines -- list of 'theories' describing curves for plotting.
    band  -- neural network 'theory' or a list of 'theories' defining 
             band by their mean and standard deviation
    xaxis -- abscissa variable; if None, last of sets.xaxes is taken;
             if 'points' all points are just put equidistantly along x-axis
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
    xn = 1
    xsets = []
    for set in sets:
        xval = []; yval = []; yerr = []
        for pt in set:
            if xaxis == 'point':
                xval.append(xn)
                xn += 1
            else:
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
        up.append(mean + std)
        down.append(mean - std)
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

def anydata(sets, lines=[], band=[], path=None, fmt='png'):
    """Plot Mock HERMES-like BSA.
    sets - list of 'Dataset' instances to be ploted
    
    """
    title = 'Any data'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    ax = fig.add_subplot(1, 1, 1)
    xaxes = ['tm', 'xB', 'Q2']
    subplot(ax, sets, lines, band, 'point', [])
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def HERMES10(lines=[], band=[], path=None, fmt='png'):
    """Plot HERMES PRELIMINARY 06/07 BCA and BSA data with fit lines."""

    #ids = [2, 4, 5]
    title = 'HERMES 10 PRELIMINARY'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    xaxes = ['tm', 'xB', 'Q2']
    ylims = [(-0.05, 0.3), (-0.15, 0.15), (-0.45, 0.05)]
    # we have 3x18=54 points to be separated in nine panels six points each:
    for y, id, shift in zip(range(3), [57, 57, 58], [18, 0, 0]):
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
                ylabels = ['$BCA\\; \\cos \\phi$', '$BCA\\; \\cos 0\\phi$', '$ALUI\\; \\sin\\phi$']
                ax.set_ylabel(ylabels[(panel-1)/3], fontsize=18)

            if panel < 7:
                # Leave labels only on lowest panels
                ax.set_xlabel('')
            else:
                xlabels = ['$-t\\; [{\\rm GeV}^2]$', '$x_B$', '$Q^2\\; [{\\rm GeV}^2]$']
                ax.set_xlabel(xlabels[panel-7], fontsize=18)

            if (panel % 3) == 2:
                # Adjust x-axis on middle column
                ax.set_xlim(0.04, 0.25)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def HERMES10BCA(lines=[], band=[], path=None, fmt='png'):
    """Plot HERMES PRELIMINARY 06/07 BCA data with fit lines."""

    #ids = [2, 4, 5]
    title = 'HERMES 10 BCA PRELIMINARY'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    xaxes = ['tm', 'xB', 'Q2']
    ylims = [(-0.15, 0.15), (-0.05, 0.3), (-0.15, 0.15)]
    # we have 3x18=54 points to be separated in nine panels six points each:
    for y, id, shift in zip(range(3), [57, 57, 57], [0, 18, 36]):
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
                ylabels = ['$BCA\\; \\cos 0\\phi$', '$BCA\\; \\cos \\phi$', '$BCA\\; \\cos 2\\phi$']
                ax.set_ylabel(ylabels[(panel-1)/3], fontsize=18)

            if panel < 7:
                # Leave labels only on lowest panels
                ax.set_xlabel('')
            else:
                xlabels = ['$-t\\; [{\\rm GeV}^2]$', '$x_B$', '$Q^2\\; [{\\rm GeV}^2]$']
                ax.set_xlabel(xlabels[panel-7], fontsize=18)

            if (panel % 3) == 2:
                # Adjust x-axis on middle column
                ax.set_xlim(0.04, 0.25)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def HERMES10BSA(lines=[], band=[], path=None, fmt='png'):
    """Plot HERMES PRELIMINARY 06/07 BSA data with fit lines."""

    #ids = [2, 4, 5]
    title = 'HERMES 10 BSA PRELIMINARY'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    xaxes = ['tm', 'xB', 'Q2']
    ylims = [(-0.45, 0.1), (-0.25, 0.25), (-0.15, 0.15)]
    # we have 3x18=54 points to be separated in nine panels six points each:
    for y, id, shift in zip(range(3), [58, 59, 58], [0, 0, 18]):
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
                ylabels = ['$ALUI\\; \\sin \\phi$', '$ALUDVCS\\; \\sin \\phi$', '$ALUI\\; \\sin 2\\phi$']
                ax.set_ylabel(ylabels[(panel-1)/3], fontsize=18)

            if panel < 7:
                # Leave labels only on lowest panels
                ax.set_xlabel('')
            else:
                xlabels = ['$-t\\; [{\\rm GeV}^2]$', '$x_B$', '$Q^2\\; [{\\rm GeV}^2]$']
                ax.set_xlabel(xlabels[panel-7], fontsize=18)

            if (panel % 3) == 2:
                # Adjust x-axis on middle column
                ax.set_xlim(0.04, 0.25)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def HERMES09(lines=[], band=[], path=None, fmt='png'):
    """Plot HERMES 0909.3587 BCA and BSA data with fit lines."""

    #ids = [2, 4, 5]
    title = '' #'HERMES 09'
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

            if (panel % 3) == 2:
                # Adjust x-axis on middle column
                ax.set_xlim(0.04, 0.25)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def HERMES10TP(obs='TSA', lines=[], band=[], path=None, fmt='png'):
    """Plot HERMES 1004.0177 TSA data with fit lines."""

    title = '' # 'HERMES10-'+obs
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    xaxes = ['tm', 'xB', 'Q2']
    ylims = [(-0.05, 0.3), (-0.15, 0.15), (-0.45, 0.05)]
    xticks = [0.2, 0.1, 2]
    if obs=='TSA':
        id = 52
        harmonics = [-1, -2, -3]
        fun = 'sin'
    else:
        id = 53
        harmonics = [0, 1, 2]
        fun = 'cos'
    subsets = {}
    for k in range(3):
        subsets[k] = utils.select(data[id], criteria=['FTn == %i' % harmonics[k]])
    # we have 3x12=36 points to be separated in nine panels four points each:
    for y in range(3):
        for x in range(3):
            panel = 3*y + x + 1  # 1, 2, ..., 9
            ax = fig.add_subplot(3,3,panel)
            ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))  # tickmarks
            ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(xticks[x]))
            subplot(ax, [subsets[y][x*4:x*4+4]], lines, band, xaxes[x], [])
            #apply(ax.set_ylim, ylims[y])
            if (panel % 3) != 1:
                # Leave labels only on leftmost panels
                ax.set_ylabel('')
            else:
                ylabels = ['$'+obs+'\\; \\'+fun+'(%i\\phi)$' % harmonics[k] for k in range(3)]
                ax.set_ylabel(ylabels[(panel-1)/3], fontsize=18)

            if panel < 7:
                # Leave labels only on lowest panels
                ax.set_xlabel('')
            else:
                xlabels = ['$-t\\; [{\\rm GeV}^2]$', '$x_B$', '$Q^2\\; [{\\rm GeV}^2]$']
                ax.set_xlabel(xlabels[panel-7], fontsize=18)

            if (panel % 3) == 2:
                # Adjust x-axis on middle column
                ax.set_xlim(0.04, 0.25)
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
    title = '' #'CLAS 07'
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

def CLASTSA(lines=[], band=[], path=None, fmt='png'):
    """Plot CLAS TSA."""
    id = 54
    title = 'CLAS-TSA'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    xaxes = ['tm', 'xB']
    # we have 2x3 points
    for x in range(2):
        ax = fig.add_subplot(2, 1, x+1)
        #ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))  # tickmarks
        subplot(ax, [data[id][x*3:x*3+3]], lines, band, xaxes[x], [])
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def HALLA(lines=[], band=[], path=None, fmt='png'):
    """Makes plot of HALL-A data with fit lines"""

    subsets = {}
    subsets[1] = utils.select(data[26], criteria=['Q2 == 2.3'])
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

def HALLAalt(lines=[], band=[], path=None, fmt='png'):
    """Makes plot of 'alternative' HALL-A data with fit lines"""

    subsets = {}
    subsets[1] = utils.select(data[50], criteria=['Q2 == 1.5', 'FTn == -1'])
    subsets[2] = utils.select(data[50], criteria=['Q2 == 1.9', 'FTn == -1'])
    subsets[3] = utils.select(data[50], criteria=['Q2 == 2.3', 'FTn == -1'])
    subsets[4] = utils.select(data[51], criteria=['FTn == 0'])
    subsets[5] = utils.select(data[51], criteria=['FTn == 1'])
    #subsets[6] = utils.select(data[51], criteria=['FTn == 2'])

    title = '' #'Hall-A-alt'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    Qs = ['1.5', '1.9', '2.3', '2.3', '2.3']
    for panel in range(1,6):
        ax = fig.add_subplot(2,3,panel)
        subplot(ax, [subsets[panel]], lines, band, 't')
        ax.set_ylabel('%s(FTn = %i)' % (subsets[panel][0].y0name, 
            subsets[panel][0].FTn), fontsize=16)
        if panel<5:
            ax.text(-0.31, 0.002, '$Q^2\\!= %s\\,{\\rm GeV}^2$' % Qs[panel-1], fontsize=12)
        else:
            ax.text(-0.31, -0.008, '$Q^2\\!= %s\\,{\\rm GeV}^2$' % Qs[panel-1], fontsize=12)
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
    # left  = 0.125  # the left side of the subplots of the figure
    # right = 0.9    # the right side of the subplots of the figure
    # bottom = 0.0   # the bottom of the subplots of the figure
    # top = 0.9      # the top of the subplots of the figure
    # wspace = 0.2   # the amount of width reserved for blank space between subplots
    # hspace = 0.2   # the amount of height reserved for white space between subplots
    fig.subplots_adjust(wspace=0.7)

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
    #title = 'H1 07 / ZEUS 08'
    title = ''
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    fig.subplots_adjust(bottom=0.1, hspace=0.3)
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
            ax.text(-0.8, 10, '${\\rm W = 82}\\, {\\rm GeV}$', fontsize=12)
            ax.text(-0.8, 3, '$Q^2\\!= 8,\\, 15.5,\\, 25\\,{\\rm GeV}^2$', fontsize=12)
        if panel==2:
            ax.text(-0.3, 6, '${\\rm ZEUS}$', fontsize=16)
            ax.text(-0.3, 3.5, '${\\rm W = 104}\\, {\\rm GeV}$', fontsize=12)
            ax.text(-0.3, 2.1, '$Q^2\\!= 3.2\\,{\\rm GeV}^2$', fontsize=12)
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
    title = '' # 'HALL A 06'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    fig.subplots_adjust(wspace=0.2, right=0.95)
    panel = 1
    for id in ids:
        ax = fig.add_subplot(2,3,panel)
        subplot(ax, [subsets[id]], lines, band, 'phi', ['Q2', 't'])
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(2.))
        ax.set_xlabel('$\\phi\\; {\\rm [rad]}$', fontsize=16)
        if not(panel==1 or panel==4):
            # Leave labels only on leftmost panels
            ax.set_ylabel('')
        elif panel==1:
            ax.set_ylabel('$d^4\\Sigma$', fontsize=18)
            #ax.set_ylabel('$d^4\\Sigma/(dQ^2\\! dx_{\\rm B}dt d\\phi)\\quad \
            #        [{\\rm nb/GeV}^4]$', fontsize=18)
        else: # panel 4
            ax.set_ylabel('$d^4\\sigma$', fontsize=18)
            #ax.set_ylabel('$d^4\\sigma/(dQ^2\\! dx_{\\rm B}dt d\\phi)\\quad  \
            #        [{\\rm nb/GeV}^4]$', fontsize=18)
        panel += 1
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def COMPASS(lines=[], xB=0.05, Q2=2, path=None, fmt='png', numbers=False):
    """Plot COMPASS BCS asymmetry and difference and summ of xs for FormFactors model ff.
    
    FIXME: kinematic completion, charge etc. must be explicit here
    """

    title = '' # 'COMPASS $Q^2$ = %s GeV$^2$, $x_B$ = %s' % (Q2, xB)
    filename = 'COMPASS-Q2-%s-xB-%s' % (Q2, xB)
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
    pt.in1polarization = -0.8
    pt.s = 2 * Mp * pt.in1energy + Mp2
    pt.xB = xB
    pt.t = -0.2
    pt.Q2 = Q2
    #ax = fig.add_subplot(2,2,1)
    ax = fig.add_subplot(1,1,1)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    phi = np.arange(0., np.pi, 0.1)
    utils.fill_kinematics(pt)
    linestyles = ['r-', 'g--', 'b-.']
    labels = ['without_HALLA', 'with_HALLA']
    pn = 0
    asym = []
    for approach in lines:
        approach.__class__.to_conventions(pt)
        approach.__class__.prepare(pt)
        line = approach.BCSA(pt, vars={'phi':np.pi - phi})
        ax.plot(phi, line, linestyles[pn], linewidth=2, label=approach.name) 
        asym.append(line)
        pn += 1
    #ax.set_ylim(0.0, 0.5)
    # axes labels
    ax.set_xlabel('$\\phi$')
    ax.set_ylabel('BCS Asymmetry')
    ax.legend().draw_frame(0)
    ## Difference  panel
    #ax = fig.add_subplot(2,2,2)
    #ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    #pn = 0
    #diff = []
    #for approach in lines:
    #    approach.__class__.to_conventions(pt)
    #    approach.__class__.prepare(pt)
    #    # nb converted to pb:
    #    line = 1000. * approach.BCSD(pt, vars={'phi':np.pi - phi})
    #    ax.plot(phi, line, linestyles[pn], linewidth=2, label=approach.name) 
    #    diff.append(line)
    #    pn += 1
    ##ax.set_ylim(0.0, 0.5)
    ## axes labels
    #ax.set_xlabel('$\\phi$')
    #ax.set_ylabel('BCS Difference  [pb/GeV^4]')
    #ax.legend().draw_frame(0)
    ## Sum pannel
    #ax = fig.add_subplot(2,2,3)
    #ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    #pn = 0
    #summ = []
    #for approach in lines:
    #    approach.__class__.to_conventions(pt)
    #    approach.__class__.prepare(pt)
    #    # nb converted to pb:
    #    line = 1000. * approach.BCSS(pt, vars={'phi':np.pi - phi})
    #    ax.plot(phi, line, linestyles[pn], linewidth=2, label=approach.name) 
    #    summ.append(line)
    #    pn += 1
    ##ax.set_ylim(0.0, 0.5)
    ## axes labels
    #ax.set_xlabel('$\\phi$')
    #ax.set_ylabel('BCS sum  [pb/GeV^4]')
    #ax.legend().draw_frame(0)
    #if numbers:
    #    pn = 0
    #    for approach in lines:
    #        mat = np.array([phi, asym[pn], diff[pn], summ[pn]]).transpose()
    #        np.savetxt('%s-%s.dat' % (filename, labels[pn]), mat, fmt='% f', 
    #                delimiter='   ')
    #        pn += 1
    if path:
        fig.subplots_adjust(hspace=0.15)
        fig.subplots_adjust(wspace=0.35)
        fig.subplots_adjust(bottom=0.30)
        fig.savefig(os.path.join(path, filename+'.'+fmt), format=fmt)
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
        linestyles = ['r-', 'g--', 'b-.', 'p:']
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
                pt.in1polarization = -0.8
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
    fig.subplots_adjust(bottom=0.35)
    # Asymmetry panel
    pt = Data.DummyPoint()
    pt.exptype = 'collider'
    pt.in1particle = 'e'
    pt.in1charge = 1
    pt.in1energy = 5.
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
    linestyles = ['r-', 'g--', 'b-.', 'p:']
    labels = [r'$\sigma^{\uparrow}$ ' + t.name for t in fits]
    pn = 0
    for approach in fits:
        approach.__class__.to_conventions(pt)
        approach.__class__.prepare(pt)
        line = approach.Xunp(pt, vars={'phi':phi})
        ax.plot(phi, line, linestyles[pn], linewidth=1, label=labels[pn]) 
        pn += 1
    labels = ['$(\\sigma^{\\uparrow} + \\sigma^{\\downarrow})/2$ ' + t.name for t in fits]
    pn = 0
    for approach in fits:
        approach.__class__.to_conventions(pt)
        approach.__class__.prepare(pt)
        lineBSS = approach.BSS(pt, vars={'phi':phi})
        ax.plot(phi, lineBSS, linestyles[pn], linewidth=3, label=labels[pn]) 
        pn += 1
    #ax.set_ylim(0.0, 0.5)
    # axes labels
    ax.set_xlabel('$\\phi\\quad {\\rm [rad]}$', fontsize=20)
    ax.set_ylabel('$\\sigma\\quad {\\rm [nb]}$', fontsize=20)
    #ax.legend(loc=(0.3,0))
    ax.legend(loc=(0.2, 0)).draw_frame(0)
    ax.text(0.7, 90., "EIC", fontsize=18)
    ax.text(0.7, 84., "$E_e = %d \\,{\\rm GeV}$" % pt.in1energy, fontsize=18)
    ax.text(0.7, 80., "$E_p = %d \\,{\\rm GeV}$" % pt.in2energy, fontsize=18)
    ax.text(0.7, 75., "$x_B = %s$" % pt.xB, fontsize=18)
    ax.text(0.7, 71., "$Q^2 = %s \\,{\\rm GeV}^2$" % pt.Q2, fontsize=18)
    ax.text(0.7, 67., "$t = %s \\,{\\rm GeV}^2$" % pt.t, fontsize=18)
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
    fig.subplots_adjust(bottom=0.15)
    # Left panel
    pt = Data.DummyPoint()
    pt.Q2 = 2.
    ax = fig.add_subplot(1,1,1)
    ax.set_xscale('log')  # x-axis to be logarithmic
    xval = np.power(10., np.arange(-3.5, 0, 0.01)) 
    colors = ['red', 'green', 'blue', 'purple']
    styles = ['-', '--', '-.', ':']
    pn = 0
    for t in theories:
        #pt.t = -0.1
        # kludge alert!
        #ImHvals = []
        #for xi in xval:
        #    if t.model.__dict__.has_key('Gepard'): t.m.g.newcall = 1
        #    ImHvals.append(t.model.ImH(pt, xi))
        #line1 = xval * np.array(ImHvals) / np.pi
        pt.t = -0.28
        ImHvals = []
        for xi in xval:
            if t.model.__dict__.has_key('Gepard'): t.m.g.newcall = 1
            ImHvals.append(t.model.ImH(pt, xi))
        line2 = xval * ImHvals / np.pi
        #ax.plot(xval, line1, color=colors[pn], linestyle=styles[pn], linewidth=2)
        ax.plot(xval, line2, color=colors[pn], linestyle=styles[pn], linewidth=3, label=t.name) 
        pn += 1
    ax.set_ylim(0.0, 0.4)
    ax.set_xlim(0.0005, 1.0)
    #plt.ylim(0.0, 0.5)
    # axes labels
    ax.set_xlabel('$\\xi = x_{\\rm B}/(2-x_{\\rm B})$', fontsize=18)
    ax.set_ylabel('$\\xi\\, H(\\xi, \\xi, t)$', fontsize=22)
    ax.legend().draw_frame(0)
    #ax.text(0.001, 0.405, "t = 0")
    #ax.text(0.001, 0.27, "$t = -0.1\\, {\\rm GeV}^2$", fontsize=15)
    ax.text(0.01, 0.25, "$t = -0.28\\, {\\rm GeV}^2$", fontsize=15)
    ax.text(0.01, 0.22, "$Q^2 = 2\\, {\\rm GeV}^2$", fontsize=15)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def Hval(theories=[], path=None, fmt='png'):
    """Makes plot of Im(cffH) in valence region.
    
    
    """
    title = '' # 'Fig 15'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    fig.subplots_adjust(bottom=0.15)
    pt = Data.DummyPoint()
    pt.Q2 = 2.
    pt.t = -0.28
    ax = fig.add_subplot(1,1,1)
    #ax.set_xscale('log')  # x-axis to be logarithmic
    xval = np.linspace(0.07, 0.38, 20) 
    colors = ['red', 'green', 'blue', 'purple']
    styles = ['-', '--', '-.', ':']
    pn = 0
    for t in theories:
        # kludge alert!
        ImHvals = []
        for xB in xval:
            xi = xB/(2.-xB)
            if t.model.__dict__.has_key('Gepard'): t.m.g.newcall = 1
            ImHvals.append(t.model.ImH(pt, xi)/np.pi)
        ax.plot(xval, ImHvals, color=colors[pn], linestyle=styles[pn], linewidth=3, label=t.name) 
        pn += 1
    ax.set_ylim(0.0, 4.5)
    ax.set_xlim(0.075, 0.38)
    #plt.ylim(0.0, 0.5)
    # axes labels
    ax.set_xlabel('$x_B$', fontsize=18)
    ax.set_ylabel('$\\Im\\! m \\mathcal{H}(x_{\\rm B}, t, Q^2)/\\pi$', fontsize=22)
    ax.legend().draw_frame(0)
    #ax.text(0.001, 0.405, "t = 0")
    ax.text(0.2, 3.2, "$t = -0.28\\, {\\rm GeV}^2$", fontsize=15)
    ax.text(0.2, 2.8, "$Q^2 = 2\\, {\\rm GeV}^2$", fontsize=15)
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
    linestyles = ['r-', 'g--', 'b-.', 'p:'] 
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
    pt.in1polarization = -0.8
    pt.s = 2 * Mp * pt.in1energy + Mp2
    pt.xB = 0.05
    pt.t = -0.2
    pt.Q2 = 2.
    ax = fig.add_subplot(1,2,2)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    phi = np.arange(0., np.pi, 0.2)
    utils.fill_kinematics(pt)
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
    """Make a band of fun defined by ndarray-valued fun."""
    pt = Data.DummyPoint()
    pt.t = -tm
    pt.xi = xvals
    pt.xB = 2*pt.xi/(1.+pt.xi)
    up = []
    down = []
    #res = xvals * fun(pt)
    res = fun(pt)
    mean = res.mean(axis=0)
    std = res.std(axis=0)
    if avg:
        up = mean + std
        down = mean - std
    else:
        up = res.max(axis=0)
        down = res.min(axis=0)
    x = plt.concatenate( (xvals, xvals[::-1]) )
    y = plt.concatenate( (up, down[::-1]) )
    ax.fill(x, y, facecolor=color, alpha=0.5)

def _axline(ax, tm, xvals, fun, **kwargs):
    """Make a line of x*fun defined by ndarray-valued fun."""
    res = []
    for xi in xvals:
        pt = Data.DummyPoint()
        pt.t = -tm
        pt.xi = xi
        pt.xB = 2*pt.xi/(1.+pt.xi)
        #res.append(xi*fun(pt))
        res.append(fun(pt))
    ax.plot(xvals, res, **kwargs)

def CFF(t, th=None, cffs=None, path=None, fmt='png', average=True):
    """Makes plot of cffs given by neural network
    
    t     -- Neural network 'theory'
    th    -- non-net theory
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
    colors = ['blue', 'red', 'brown', 'purple']
    linestyles = ['solid', 'dashed']
    #tms = [0.1, 0.3]
    tms = [0.2]
    # Define abscissas
    logxvals = np.power(10., np.arange(-3.0, -0.01, 0.1))  # left panel
    xvals = np.linspace(0.025, 0.2, 20) # right panel
    # ordinatas for  left ...
    #allylims = [(-0.3, 1.0), (-0.3, 1.0), (-0.45, 0.05)]
    #allylims = [(-3.3, 35.0), (-3.3, 10.0), (-4.45, 0.05)]
    # ... and right panles
    #ylims = [(-0.3, 1.0), (-0.3, 1.0), (-0.45, 0.05)]
    #ylims = [(-3.3, 35.0), (-3.3, 10.0), (-4.45, 0.05)]
    # Plot panels
    for n in range(len(cffs)):
        cff = cffs[n]
        # all-x logarithmic
        ax = fig.add_subplot(len(cffs), 2, 2*n+1)
        ax.set_xscale('log')  # x-axis to be logarithmic
        for tm in tms :
            _axband(ax, tm, logxvals, getattr(t.model, cff), 
                    color=colors[tms.index(tm)], avg=average)
            if th:
                _axline(ax, tm, logxvals, getattr(th.model, cff), 
                        color=colors[tms.index(tm)], 
                        linestyle='dashed', label='t = %s' % str(tm))
        ax.set_xlabel('$\\xi$', fontsize=15)
        ax.set_ylabel('%s' % cff, fontsize=18)
        ax.axhspan(-0.0005, 0.0005, facecolor='g', alpha=0.6)  # horizontal bar
        ax.axvspan(0.03, 0.093, facecolor='g', alpha=0.1)  # vertical band
        ax.text(0.03, -0.27, "data region", fontsize=14)
        #apply(ax.set_ylim, allylims[n])
        ax.set_xlim(0.005, 1.0)
        # measured x linear
        ax = fig.add_subplot(len(cffs), 2, 2*n+2)
        for tm in tms:
            _axband(ax, tm, xvals, getattr(t.model, cff), 
                    color=colors[tms.index(tm)], avg=average)
            if th:
                _axline(ax, tm, xvals, getattr(th.model, cff), 
                        color=colors[tms.index(tm)], 
                        linestyle='dashed', label='t = %s' % str(tm))
        ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
        ax.set_xlabel('$\\xi$', fontsize=15)
        #apply(ax.set_ylim, ylims[n])
        ax.set_xlim(0.03, 0.093)
        ax.text(0.03, -0.27, "data region only", fontsize=14)
        #ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.02))  # tickmarks
    t.model.parameters['nnet'] = old
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig


def HvalNN(theories=[], band=None, path=None, fmt='png'):
    """Makes plot of Im(cffH) in valence region.
    
    
    """
    title = '' # 'Fig 15'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    fig.subplots_adjust(bottom=0.15)
    pt = Data.DummyPoint()
    pt.Q2 = 2.
    pt.t = -0.28
    ax = fig.add_subplot(1,1,1)
    #ax.set_xscale('log')  # x-axis to be logarithmic
    xval = np.linspace(0.07, 0.38, 20) 
    colors = ['red', 'green', 'blue', 'purple']
    styles = ['-', '--', '-.', ':']
    pn = 0
    for t in theories:
        # kludge alert!
        ImHvals = []
        for xB in xval:
            xi = xB/(2.-xB)
            if t.model.__dict__.has_key('Gepard'): t.m.g.newcall = 1
            ImHvals.append(t.model.ImH(pt, xi)/np.pi)
        ax.plot(xval, ImHvals, color=colors[pn], linestyle=styles[pn], linewidth=3, label=t.name) 
        pn += 1
    ax.set_ylim(0.0, 4.5)
    ax.set_xlim(0.075, 0.38)
    #plt.ylim(0.0, 0.5)
    # axes labels
    ax.set_xlabel('$x_B$', fontsize=18)
    ax.set_ylabel('$\\Im\\! m \\mathcal{H}(x_{\\rm B}, t, Q^2)/\\pi$', fontsize=22)
    ax.legend().draw_frame(0)
    #ax.text(0.001, 0.405, "t = 0")
    ax.text(0.2, 3.2, "$t = -0.28\\, {\\rm GeV}^2$", fontsize=15)
    ax.text(0.2, 2.8, "$Q^2 = 2\\, {\\rm GeV}^2$", fontsize=15)
    if band:
        old = band.model.parameters['nnet']
        band.model.parameters['nnet'] = 'ALL'
        _axband(ax, -pt.t, xval, band.m.ImH, color='green', avg=False)
        band.model.parameters['nnet'] = old
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig


def Ht(theories=[], xB=0.36, path=None, fmt='png'):
    """Makes plot of cffH for M. Guidal's proceedings.
    
    """
    title = ''
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    fig.subplots_adjust(bottom=0.15)
    pt = Data.DummyPoint()
    pt.Q2 = 2.
    pt.xB = xB
    pt.xi = xB/(2.-xB)
    ax = fig.add_subplot(1,1,1)
    #ax.set_xscale('log')  # x-axis to be logarithmic
    mts = np.linspace(0.0, 0.6, 20) 
    colors = ['red', 'green', 'blue', 'purple']
    styles = ['-', '--', '-.', ':']
    pn = 0
    for th in theories:
        # kludge alert!
        ImHvals = []
        ReHvals = []
        for mt in mts:
            pt.mt = mt
            pt.t = -mt
            if th.model.__dict__.has_key('Gepard'): th.m.g.newcall = 1
            ImHvals.append(th.model.ImH(pt)/np.pi)  # M. Guidal conventions!
            ReHvals.append(-th.model.ReH(pt))       # M. Guidal conventions!
        ax.plot(mts, ImHvals, color=colors[pn], linestyle=styles[pn], linewidth=3, label=th.name) 
        ax.plot(mts, ReHvals, color=colors[pn], linestyle=styles[pn], linewidth=1, label=th.name) 
        pn += 1
    ax.set_ylim(-2.0, 11.0)
    #ax.set_xlim(0.075, 0.38)
    #plt.ylim(0.0, 0.5)
    # axes labels
    ax.set_xlabel('$-t \\, [{\\rm GeV}^2]$', fontsize=18)
    ax.set_ylabel('$\\Im\\! m \\mathcal{H}(x_{\\rm B}, t, Q^2=2\\, {\\rm GeV}^2)$', fontsize=22)
    ax.legend().draw_frame(0)
    #ax.text(0.001, 0.405, "t = 0")
    ax.text(0.2, 3.2, "$xB = %s$" % xB, fontsize=15)
    #ax.text(0.2, 2.8, "$Q^2 = 2\\, {\\rm GeV}^2$", fontsize=15)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def xBt(sets, path=None, fmt='png'):
    """Plot xB-t distribution of data sets. """
    
    title = 'xB-t distribution of data'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    ax = fig.add_subplot(1, 1, 1)
    setshapes = ['o', 's', '^', 'd']  # first circles, then squares ...
    setcolors = ['blue', 'black', 'purple', 'green']  # circles are blue, squares are black, ...
    setn = 0
    for set in sets:
        xval = []; yval = []; err = []
        for pt in set:
            xval.append(pt.xB)
            yval.append(pt.t)
            err.append(abs(pt.err/pt.val))
        nperr = np.array(err)
        ax.plot(xval, yval, linestyle='None', markersize=math.log(5e3*nperr.mean()),
                marker=setshapes[setn], color=setcolors[setn])
        setn += 1
    ax.set_xlabel(toTeX['xB'], fontsize=15)
    ax.set_ylabel(toTeX['t'], fontsize=15)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def uncertCFF(tnn=None, th=None, thband=None, cffs=['ImH'], path=None, fmt='png'):
    """Makes plot of cffs given by theory th, with uncertainty band.
    
    tnn     -- neural net (NN) theory
    th      -- non-NN theory to be plotted as line
    thband  -- non-NN theory to be plotted as 1-sigma error band
    cffs    -- List of CFFs to be plotted. Each produces two panels.

    """
    if tnn:
        title = tnn.description
        old = tnn.model.parameters['nnet']
        tnn.model.parameters['nnet'] = 'ALL'
    elif thband:
        title = thband.description
    elif th:
        title = th.description
    else:
        title = 'uncertCFF'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    colors = ['red', 'brown']     # worm human colors :-)
    nncolors = ['blue', 'green']  # cold computer colors
    linestyles = ['solid', 'dashed']
    tms = [0.2]
    # Define abscissas
    logxvals = np.power(10., np.arange(-3.0, -0.01, 0.1))  # left panel
    xvals = np.linspace(0.025, 0.2, 20) # right panel
    # ordinatas for  left ...
    #allylims = [(-0.3, 1.0), (-0.3, 1.0), (-0.45, 0.05)]
    #allylims = [(-3.3, 35.0), (-3.3, 10.0), (-4.45, 0.05)]
    # ... and right panles
    #ylims = [(-0.3, 1.0), (-0.3, 1.0), (-0.45, 0.05)]
    #ylims = [(-3.3, 35.0), (-3.3, 10.0), (-4.45, 0.05)]
    # Plot panels
    for n in range(len(cffs)):
        cff = cffs[n]
        # all-x logarithmic
        ax = fig.add_subplot(len(cffs), 2, 2*n+1)
        ax.set_xscale('log')  # x-axis to be logarithmic
        for tm in tms :
            if tnn:
                _axband(ax, tm, logxvals, getattr(tnn.model, cff), 
                        color=nncolors[tms.index(tm)])
            if th:
                _axline(ax, tm, logxvals, getattr(th.model, cff), 
                        color=colors[tms.index(tm)], 
                        linestyle='dashed', label='t = %s' % str(tm))
            if thband:
                _axband(ax, tm, logxvals, lambda pt: thband.model.uncert(cff, pt), 
                        color=colors[tms.index(tm)])
        ax.set_xlabel('$\\xi$', fontsize=15)
        ax.set_ylabel('%s' % cff, fontsize=18)
        ax.axhspan(-0.0005, 0.0005, facecolor='g', alpha=0.6)  # horizontal bar
        ax.axvspan(0.03, 0.093, facecolor='g', alpha=0.1)  # vertical band
        ax.text(0.03, -0.27, "data region", fontsize=14)
        #apply(ax.set_ylim, allylims[n])
        ax.set_xlim(0.005, 1.0)
        # measured x linear
        ax = fig.add_subplot(len(cffs), 2, 2*n+2)
        for tm in tms:
            if tnn:
                _axband(ax, tm, xvals, getattr(tnn.model, cff), 
                        color=nncolors[tms.index(tm)])
            if th:
                _axline(ax, tm, xvals, getattr(th.model, cff), 
                        color=colors[tms.index(tm)], 
                        linestyle='dashed', label='t = %s' % str(tm))
            if thband:
                _axband(ax, tm, xvals, lambda pt: thband.model.uncert(cff, pt), 
                        color=colors[tms.index(tm)])
        ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
        ax.set_xlabel('$\\xi$', fontsize=15)
        #apply(ax.set_ylim, ylims[n])
        ax.set_xlim(0.03, 0.093)
        ax.text(0.03, -0.27, "data region only", fontsize=14)
        #ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.02))  # tickmarks
    if tnn:
        tnn.model.parameters['nnet'] = old
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

