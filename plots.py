""" 
Plotting functions written for ad-hoc plotting
of some specific datasets and CFFs.

"""

#from IPython.Debugger import Tracer; debug_here = Tracer()

import sys, os, math, copy, string, commands
import numpy as np

import matplotlib
from matplotlib.backends.backend_pdf import PdfPages # for PDF creation
if os.sys.platform == 'win32':
    matplotlib.use('WxAgg')
#else: #linux
#    matplotlib.use('TkAgg')
import pylab as plt


import Data, Approach, utils
from constants import toTeX, Mp2, Mp, OBStoTeX
from results import *

# load experimental data
data = utils.loaddata('/home/kkumer/pype/data/ep2epgamma', approach=Approach.BMK) 
data.update(utils.loaddata('/home/kkumer/pype/data/gammastarp2gammap', approach=Approach.BMK)) 

###  subplots_adjust options and their meanings:
 # left  = 0.125  # the left side of the subplots of the figure
 # right = 0.9    # the right side of the subplots of the figure
 # bottom = 0.0   # the bottom of the subplots of the figure
 # top = 0.9      # the top of the subplots of the figure
 # wspace = 0.2   # the amount of width reserved for blank space between subplots
 # hspace = 0.2   # the amount of height reserved for white space between subplots
 # Usage: fig.subplots_adjust(wspace=0.7)

def mkpdf(filename):
    """Create publishable PDF and EPS figs."""
    root = filename[:-4]  # assumed .pdf !
    pp = PdfPages(filename)
    pp.savefig()
    pp.close()
    # fix PDF's bounding box
    commands.getoutput('pdfcrop --margins "0" --clip %s' % filename)
    commands.getoutput('mv %s-crop.pdf %s' % (root, filename))
    # create EPS
    commands.getoutput('acroread -toPostScript %s' % filename)
    commands.getoutput('ps2eps -f %s.ps' % root)
    commands.getoutput('rm %s.ps' % root)
    # fix EPS's bounding box
    commands.getoutput('mv %s.eps tmpfile.eps' % root)
    commands.getoutput('epstool --copy --bbox tmpfile.eps %s.eps' % root)
    commands.getoutput('rm tmpfile.eps')


#################################################################
##                                                             ##
##  [1]  Universal stuff, used to build up panels of figures   ##
##                                                             ##
#################################################################

def _axpoints(ax, pts, xaxis, **kwargs):
    """Make an errorbar plot defined by a single set of data points.

    pts   -- single set of data points
    xaxis -- 'Q2', 't', 'xB', ...

    """
    kwargs.pop('justbars', False)  # Kludge
    pts = [Approach.BMK.from_conventions(pt.copy()) for pt in pts]
    xvals = [getattr(pt, xaxis) for pt in pts]
    yvals = [pt.val for pt in pts]
    yerrs = [pt.err for pt in pts]
    ax.errorbar(xvals, yvals, yerrs, **kwargs)

def _axline(ax, fun, points, xaxis, **kwargs):
    """Make a line(s) corresponding to points 

    fun -- real-valued function of DataPoint instance
    points -- list of sets of DataPoint instances
    xaxis -- 'Q2', 't', 'xB', ...

    """
    for pts in points:
        xvals = [getattr(pt, xaxis) for pt in pts]
        yvals = [fun(pt) for pt in pts]
        kwargs.pop('justbars', False)  # Kludge
        if kwargs.pop('xF', False):  # do we want  xF(x) plotted?
            xis = np.array([pt.xi for pt in pts])
            yvals = xis*np.array(yvals)
        ax.plot(xvals, yvals, **kwargs)

def _axband(ax, fun, pts, xaxis, **kwargs):
    """Make a band corresponding to points 

    fun -- function of DataPoint instance returning tuple
           (mean, err) or (mean, err+, err-) defining errorband
    pts -- set of DataPoint instances (just one set!)
    xaxis -- 'Q2', 't', 'xB', ...

    """
    if not (isinstance(pts[0], Data.DataPoint) or isinstance(pts[0], Data.DummyPoint)):
        raise ValueError, "%s is not single dataset" % str(pts)

    xvals = [getattr(pt, xaxis) for pt in pts]
    res = [fun(pt) for pt in pts]
    if len(res[0]) == 2:
        # we have symmetric error
        res = [(m, err, err) for m, err in res]
    up, down = np.array([(m+errp, m-errm) for m,errp,errm in res]).transpose()
    if kwargs.pop('justbars', False):  # do we want just errorbars for theory?
        res = zip(up, down)
        #FIXME: symmetric error assumed
        yvals, yerrs = np.array([((up+down)/2, (up-down)/2) for up,down in res]).transpose()
        xvals = [x+0.1 for x in xvals]
        ax.errorbar(xvals, yvals, yerrs, marker='d', color='red', linestyle='None', 
                elinewidth=3)
    else:
        x = plt.concatenate( (xvals, xvals[::-1]) )
        y = plt.concatenate( (up, down[::-1]) )
        if kwargs.pop('xF', False):  # do we want  xF(x) plotted?
            xis = np.array([pt.xi for pt in pts+pts[::-1]])
            y = xis*y
        ax.fill(x, y, alpha=0.5, **kwargs)


def panel(ax, points=None, lines=None, bands=None, xaxis=None, xs=None, 
        kins={}, kinlabels=None, CL=False, **kwargs):
    """Plot datapoints together with fit/theory line(s) and band(s).

    ax -- subplot of matplotlib's figure i.e. ax = figure.add_subplot(..)

    Keyword arguments:

    points --  (list of) `DataSet` instance(s) to be plotted.
    lines -- (list of) 'theorie(s)' describing curves for plotting.
    bands  -- (list of) neural network 'theorie(s)' or 'theorie(s)' with
              defined uncertianties, defining band by their mean and std.dev.
    xaxis -- abscissa variable; if None, last of sets.xaxes is taken;
             if 'points' all points are just put equidistantly along x-axis
    xs    -- list of xvalues, needed if not specified by points

    Additional keywoard arguments:

    kinlabels -- list of constant kinematic variables whose values will
                 be extracted from points and put on plot as annotation
    CL -- if 68% C.L. is required for bands instead of one sigma

    **kwargs -- passed to matplotlib
    TODO: legend

    """
    if points:
        # If not provided, take xaxis from last axis from last set
        if not xaxis: xaxis = points[-1].xaxes[-1]
        if isinstance(points[0], Data.DataPoint): points = [points]

        pointshapes = ['o', 's', '^', 'd', 'o', 's', '^', 'd']  # first circles, then squares ...
        pointcolors = ['blue', 'black', 'purple', 'green', 'red', 'brown', 'blue', 'black']
        setn = 0
        for pts in points:
            _axpoints(ax, pts, xaxis, linestyle='None', elinewidth=1, 
                    marker=pointshapes[setn], color=pointcolors[setn], **kwargs)
            setn += 1
    else:
        # We have to create set of DataPoints.
        # (It is up to caller to provide everything needed
        # via xaxis and kins)
        pts = []
        for x in xs:
            pt = Data.DummyPoint(init=kins)
            setattr(pt, xaxis, x)
            # workaround for imperfect fill_kinematics
            if hasattr(pt, 'xi'):
                pt.xB = 2*pt.xi/(1.+pt.xi)
            # end of workaorund
            utils.fill_kinematics(pt)
            pts.append(pt)
        points = [pts]

    if lines:
        if not isinstance(lines, list): lines = [lines]
        lineshapes = ['s', '^', 'd', 'h']  # first squares, then triangles, diamonds, hexagons
        linecolors = ['black', 'blue', 'green', 'purple']  
        linestyles = ['-', '--', '-.', ':']  # solid, dashed, dot-dashed, dotted
        linen = 0
        for line in lines:
            _axline(ax, lambda pt: line.predict(pt, orig_conventions=True), points, xaxis=xaxis,
                    color=linecolors[linen], linestyle=linestyles[linen], 
                    linewidth=2, label=line.name, **kwargs)
            linen += 1

    if bands:
        if not isinstance(bands, list): bands = [bands]
        bandcolors = ['red', 'green', 'blue', 'purple']
        hatches = ['//', '\\\\', '|', '.']
        bandn = 0
        for band in bands:
            for pts in points:
                _axband(ax, lambda pt: band.predict(pt, error=True, CL=CL, orig_conventions=True), pts, xaxis=xaxis,
                        hatch=hatches[bandn], 
                        #color=bandcolors[bandn],
                        #facecolor=bandcolors[bandn], # band  color
                        facecolor='none',
                        edgecolor=bandcolors[bandn], # band edge color
                        linewidth=1,
                        label=band.name, **kwargs)
            bandn += 1

    if kinlabels:
        # constant kinematic variables positioning
        xvals = []; yvals = []
        for pts in points:
            for pt in pts:
                xvals.append(getattr(pt, xaxis))
                yvals.append(pt.val)
        labx = min(0, min(xvals)) + (max(xvals) - min(0, min(xvals))) * 0.35
        laby = min(0, min(yvals)) + (max(yvals) - min(0, min(yvals))) * 0.02
        labtxt = ""
        for lab in kinlabels:
            try:
                labtxt += toTeX[lab] + ' = ' + str(getattr(points[0],lab)) + ', '
            except AttributeError:
                # If dataset doesn't have it, all points should have it 
                labtxt += toTeX[lab] + ' = ' + str(getattr(points[0][0],lab)) + ', '
        ax.text(labx, laby, labtxt[:-2])

#################################################################
##                                                             ##
##  [2]  Comparisons with various measurements                 ##
##                                                             ##
#################################################################

def HERMES09BCA(path=None, fmt='png', **kwargs):
    """Plot HERMES 09 BCA."""
    id = 32
    title = 'HERMES 09 BCA'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    xaxes = ['tm', 'xB', 'Q2']
    # we have 3x6 points
    for x in range(3):
        ax = fig.add_subplot(3, 1, x+1)
        ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))  # tickmarks
        panel(ax, points=data[id][x*6+18:x*6+6+18], xaxis=xaxes[x], **kwargs)
        #apply(ax.set_ylim, ylims[(-0.30, 0.05)])
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def HERMES12(path=None, fmt='png', **kwargs):
    """Plot HERMES combined BCA and BSA data with fit lines."""

    title = 'HERMES-12'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    xaxes = ['tm', 'xB', 'Q2']
    ylims = [(-0.05, 0.3), (-0.15, 0.15), (-0.45, 0.05)]
    # we have 3x18=54 points to be separated in nine panels six points each:
    for y, id, shift in zip(range(3), [67, 67, 68], [18, 0, 0]):
        for x in range(3):
            npanel = 3*y + x + 1  # 1, 2, ..., 9
            ax = fig.add_subplot(3,3,npanel)
            ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))  # tickmarks
            panel(ax, points=data[id][x*6+shift:x*6+6+shift], xaxis=xaxes[x], **kwargs)
            apply(ax.set_ylim, ylims[y])
            if (npanel % 3) != 1:
                # Leave labels only on leftmost panels
                ax.set_ylabel('')
            else:
                ylabels = ['$BCA\\; \\cos \\phi$', '$BCA\\; \\cos 0\\phi$', '$ALUI\\; \\sin\\phi$']
                ax.set_ylabel(ylabels[(npanel-1)/3], fontsize=18)
            if npanel < 7:
                # Leave labels only on lowest panels
                ax.set_xlabel('')
            else:
                xlabels = ['$-t\\; [{\\rm GeV}^2]$', '$x_B$', '$Q^2\\; [{\\rm GeV}^2]$']
                ax.set_xlabel(xlabels[npanel-7], fontsize=18)

            if (npanel % 3) == 2:
                # Adjust x-axis on middle column
                ax.set_xlim(0.04, 0.25)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def HERMES12BCA(path=None, fmt='png', **kwargs):
    """Plot HERMES combined BCA data with fit lines."""

    #ids = [2, 4, 5]
    title = 'HERMES-12-BCA'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    xaxes = ['tm', 'xB', 'Q2']
    ylims = [(-0.15, 0.15), (-0.05, 0.3), (-0.15, 0.15)]
    # we have 3x18=54 points to be separated in nine panels six points each:
    for y, id, shift in zip(range(3), [67, 67, 67], [0, 18, 36]):
        for x in range(3):
            npanel = 3*y + x + 1  # 1, 2, ..., 9
            ax = fig.add_subplot(3,3,npanel)
            ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))  # tickmarks
            panel(ax, points=data[id][x*6+shift:x*6+6+shift], xaxis=xaxes[x], **kwargs)
            apply(ax.set_ylim, ylims[y])
            if (npanel % 3) != 1:
                # Leave labels only on leftmost panels
                ax.set_ylabel('')
            else:
                ylabels = ['$BCA\\; \\cos 0\\phi$', '$BCA\\; \\cos \\phi$', '$BCA\\; \\cos 2\\phi$']
                ax.set_ylabel(ylabels[(npanel-1)/3], fontsize=18)
            if npanel < 7:
                # Leave labels only on lowest panels
                ax.set_xlabel('')
            else:
                xlabels = ['$-t\\; [{\\rm GeV}^2]$', '$x_B$', '$Q^2\\; [{\\rm GeV}^2]$']
                ax.set_xlabel(xlabels[npanel-7], fontsize=18)
            if (npanel % 3) == 2:
                # Adjust x-axis on middle column
                ax.set_xlim(0.04, 0.25)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def HERMES12BSA(path=None, fmt='png', **kwargs):
    """Plot HERMES combined BSA data with fit lines."""

    #ids = [2, 4, 5]
    title = 'HERMES-12-BSA'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    xaxes = ['tm', 'xB', 'Q2']
    ylims = [(-0.45, 0.1), (-0.25, 0.25), (-0.15, 0.15)]
    # we have 3x18=54 points to be separated in nine panels six points each:
    for y, id, shift in zip(range(3), [68, 69, 68], [0, 0, 18]):
        for x in range(3):
            npanel = 3*y + x + 1  # 1, 2, ..., 9
            ax = fig.add_subplot(3,3,npanel)
            ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))  # tickmarks
            panel(ax, points=data[id][x*6+shift:x*6+6+shift], xaxis=xaxes[x], **kwargs)
            apply(ax.set_ylim, ylims[y])
            if (npanel % 3) != 1:
                # Leave labels only on leftmost panels
                ax.set_ylabel('')
            else:
                ylabels = ['$ALUI\\; \\sin \\phi$', '$ALUDVCS\\; \\sin \\phi$', '$ALUI\\; \\sin 2\\phi$']
                ax.set_ylabel(ylabels[(npanel-1)/3], fontsize=18)
            if npanel < 7:
                # Leave labels only on lowest panels
                ax.set_xlabel('')
            else:
                xlabels = ['$-t\\; [{\\rm GeV}^2]$', '$x_B$', '$Q^2\\; [{\\rm GeV}^2]$']
                ax.set_xlabel(xlabels[npanel-7], fontsize=18)
            if (npanel % 3) == 2:
                # Adjust x-axis on middle column
                ax.set_xlim(0.04, 0.25)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def HERMES09(path=None, fmt='png', **kwargs):
    """Plot HERMES 0909.3587 BCA and BSA data with fit lines."""

    #ids = [2, 4, 5]
    title = ''#HERMES-09'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    xaxes = ['tm', 'xB', 'Q2']
    ylims = [(-0.05, 0.3), (-0.45, 0.05)]
    xlims = [(0.0, 0.48), (0.04, 0.27), (0.85, 6)]
    # we have 3x18=54 points to be separated in nine panels six points each:
    for y, id, shift in zip(range(2), [32, 5], [18, 0]):
        for x in range(3):
            npanel = 3*y + x + 1  # 1, 2, ..., 6
            ax = fig.add_subplot(2,3,npanel)
            ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))  # tickmarks
            panel(ax, points=data[id][x*6+shift:x*6+6+shift], xaxis=xaxes[x], **kwargs)
            apply(ax.set_ylim, ylims[y])
            apply(ax.set_xlim, xlims[x])
            if (npanel % 3) != 1:
                # Leave labels only on leftmost panels
                ax.set_ylabel('')
                ax.set_yticklabels([])
            else:
                ylabels = ['$A_{C}^{\\cos\\phi}$', '$A_{LU}^{\\sin\\phi}$']
                ax.set_ylabel(ylabels[(npanel-1)/3], fontsize=20)
            if npanel < 4:
                # Leave labels only on lowest panels
                ax.set_xlabel('')
                ax.set_xticklabels([])
            else:
                xlabels = ['$-t\\; [{\\rm GeV}^2]$', '$x_B$', '$Q^2\\; [{\\rm GeV}^2]$']
                ax.set_xlabel(xlabels[npanel-7], fontsize=18)
            if npanel == 3:
                ax.legend(bbox_to_anchor=(0., 1.), loc='upper right',
                        borderaxespad=0.).draw_frame(0)
                pass
        #ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.2))  # tickmarks
            for label in ax.get_xticklabels() + ax.get_yticklabels():
                label.set_fontsize(14)
    fig.subplots_adjust(bottom=0.45, wspace=0.0, hspace=0.0)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def HERMES10LP(obs='TSA', path=None, fmt='png', **kwargs):
    """Plot HERMES 1004.0177 TSA data with fit lines."""

    title = 'HERMES-10-'+obs
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    xaxes = ['tm', 'xB', 'Q2']
    ylims = [(-0.05, 0.3), (-0.15, 0.15), (-0.45, 0.05)]
    xticks = [0.2, 0.1, 2]
    if obs == 'TSA':
        id = 52
        harmonics = [-1, -2, -3]
        fun = 'sin'
    elif obs == 'BTSA':
        id = 53
        harmonics = [0, 1, 2]
        fun = 'cos'
    else:
        raise ValueError, 'Observable %s nonexistent.' % obs
    subsets = {}
    for k in range(3):
        subsets[k] = utils.select(data[id], criteria=['FTn == %i' % harmonics[k]])
    # we have 3x12=36 points to be separated in nine panels four points each:
    for y in range(3):
        for x in range(3):
            npanel = 3*y + x + 1  # 1, 2, ..., 9
            ax = fig.add_subplot(3,3,npanel)
            ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))  # tickmarks
            ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(xticks[x]))
            panel(ax, points=subsets[y][x*4:x*4+4], xaxis=xaxes[x], **kwargs)
            #apply(ax.set_ylim, ylims[y])
            if (npanel % 3) != 1:
                # Leave labels only on leftmost panels
                ax.set_ylabel('')
            else:
                ylabels = ['$'+obs+'\\; \\'+fun+'(%i\\phi)$' % harmonics[k] for k in range(3)]
                ax.set_ylabel(ylabels[(npanel-1)/3], fontsize=18)
            if npanel < 7:
                # Leave labels only on lowest panels
                ax.set_xlabel('')
            else:
                xlabels = ['$-t\\; [{\\rm GeV}^2]$', '$x_B$', '$Q^2\\; [{\\rm GeV}^2]$']
                ax.set_xlabel(xlabels[npanel-7], fontsize=18)
            if (npanel % 3) == 2:
                # Adjust x-axis on middle column
                ax.set_xlim(0.04, 0.25)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def HERMES08TP(path=None, fmt='png', **kwargs):
    """Plot HERMES 08 TTSA data with fit lines."""

    title = 'HERMES-08 TTSA'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    xaxes = ['tm', 'xB', 'Q2']
    ylims = [(-0.42, 0.15), (-0.42, 0.15)]
    xlims = [(0.0, 0.48), (0.04, 0.27), (0.85, 6)]
    # we have 2x12=24 points to be separated in six panels four points each:
    for y, id, shift in zip(range(2), [66, 65], [12, 0]):
        for x in range(3):
            npanel = 3*y + x + 1  # 1, 2, ..., 6
            ax = fig.add_subplot(2,3,npanel)
            ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))  # tickmarks
            panel(ax, points=data[id][x*4+shift:x*4+4+shift], xaxis=xaxes[x], **kwargs)
            apply(ax.set_ylim, ylims[y])
            apply(ax.set_xlim, xlims[x])
            if (npanel % 3) != 1:
                # Leave labels only on leftmost panels
                ax.set_ylabel('')
                ax.set_yticklabels([])
            else:
                ylabels = ['$A_{UT,I}^{\\sin(\\phi-\\phi_S)\\cos\\phi}$', 
                           '$A_{UT,DVCS}^{\\sin(\\phi-\\phi_S)}$']
                ax.set_ylabel(ylabels[(npanel-1)/3], fontsize=20)
            if npanel < 4:
                # Leave labels only on lowest panels
                ax.set_xlabel('')
                ax.set_xticklabels([])
            else:
                xlabels = ['$-t\\; [{\\rm GeV}^2]$', '$x_B$', '$Q^2\\; [{\\rm GeV}^2]$']
                ax.set_xlabel(xlabels[npanel-7], fontsize=18)
            if npanel == 3 and kwargs:
                ax.legend(bbox_to_anchor=(0., 1.), loc='upper right',
                        borderaxespad=0.).draw_frame(0)
                pass
            ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.2))  # tickmarks
            ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
            for label in ax.get_xticklabels() + ax.get_yticklabels():
                label.set_fontsize(14)
    fig.subplots_adjust(bottom=0.35, wspace=0.0, hspace=0.0)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def CLAS(path=None, fmt='png', **kwargs):
    """Makes plot of CLAS BSA data with fit lines and bands"""

    dataset = data[25]
    title = 'CLAS-07'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    nmax = len(dataset) - 1
    # Each different Q2 has its own panel
    npt = 0
    panelpoints = []
    for npanel in [7, 8, 9, 4, 5, 6, 1, 2, 3]:
        ax = fig.add_subplot(3,3,npanel)
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
        if npanel == 1 or npanel == 3:
            # fake doubling of points to get line visibility
            ptd = copy.deepcopy(panelset[0])
            ptd.tm = ptd.tm+0.02
            panelset.append(ptd)
            panel(ax, points=panelset, xaxis='tm', kinlabels=['Q2', 'xB'], **kwargs)
        else:
            panel(ax, points=panelset, xaxis='tm', kinlabels=['Q2', 'xB'], **kwargs)

        plt.xlim(0.0, 0.6)
        plt.ylim(0.0, 0.4)
        ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
        if (npanel % 3) != 1:
            # Leave labels only on leftmost panels
            ax.set_ylabel('')
        if npanel < 7:
            ax.set_xlabel('')
        else:
            ax.set_xlabel('$-t\\; [{\\rm GeV}^2]$', fontsize=18)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def CLASTSA(path=None, fmt='png', **kwargs):
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
        panel(ax, points=data[id][x*3:x*3+3], xaxis=xaxes[x], **kwargs)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def HallAFT(path=None, fmt='png', **kwargs):
    """Makes plot of harmonics of Hall-A data with fit lines"""

    subsets = {}
    subsets[1] = utils.select(data[50], criteria=['Q2 == 1.5', 'FTn == -1'])
    subsets[2] = utils.select(data[50], criteria=['Q2 == 1.9', 'FTn == -1'])
    subsets[3] = utils.select(data[50], criteria=['Q2 == 2.3', 'FTn == -1'])
    subsets[4] = utils.select(data[51], criteria=['FTn == 0'])
    subsets[5] = utils.select(data[51], criteria=['FTn == 1'])
    #subsets[6] = utils.select(data[51], criteria=['FTn == 2'])

    title = 'Hall-A'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    Qs = ['1.5', '1.9', '2.3', '2.3', '2.3']
    for npanel in range(1,6):
        ax = fig.add_subplot(2,3,npanel)
        panel(ax, points=subsets[npanel], xaxis='t', **kwargs)
        ax.set_ylabel('%s(FTn = %i)' % (subsets[npanel][0].y1name, 
            subsets[npanel][0].FTn), fontsize=16)
        if npanel<5:
            ax.text(-0.31, 0.002, '$Q^2\\!= %s\\,{\\rm GeV}^2$' % Qs[npanel-1], fontsize=12)
        else:
            ax.text(-0.31, -0.008, '$Q^2\\!= %s\\,{\\rm GeV}^2$' % Qs[npanel-1], fontsize=12)
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
    fig.subplots_adjust(wspace=0.7)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def HallAphi(path=None, fmt='png', **kwargs):
    """Makes plot of Hall-A data phi-dependence with fit lines"""

    ids = [9, 14, 20, 21, 23, 24]
    subsets = {}
    subsets[9] = utils.select(data[33], criteria=['Q2 == 1.5', 't == -0.17'])
    subsets[14] = utils.select(data[33], criteria=['Q2 == 1.9', 't == -0.23'])
    subsets[20] = utils.select(data[33], criteria=['Q2 == 2.3', 't == -0.33'])
    subsets[21] = utils.select(data[34], criteria=['t == -0.17'])
    subsets[23] = utils.select(data[34], criteria=['t == -0.28'])
    subsets[24] = utils.select(data[34], criteria=['t == -0.33'])
    title = 'HALL-A-06'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    fig.subplots_adjust(wspace=0.2, right=0.95)
    npanel = 1
    for id in ids:
        ax = fig.add_subplot(2,3,npanel)
        panel(ax, points=subsets[id], xaxis='phi', kinlabels=['Q2', 't'], **kwargs)
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(2.))
        ax.set_xlabel('$\\phi\\; {\\rm [rad]}$', fontsize=16)
        if not(npanel==1 or npanel==4):
            # Leave labels only on leftmost panels
            ax.set_ylabel('')
        elif npanel==1:
            ax.set_ylabel('$d^4\\Sigma$', fontsize=18)
            #ax.set_ylabel('$d^4\\Sigma/(dQ^2\\! dx_{\\rm B}dt d\\phi)\\quad \
            #        [{\\rm nb/GeV}^4]$', fontsize=18)
        else: # panel 4
            ax.set_ylabel('$d^4\\sigma$', fontsize=18)
            #ax.set_ylabel('$d^4\\sigma/(dQ^2\\! dx_{\\rm B}dt d\\phi)\\quad  \
            #        [{\\rm nb/GeV}^4]$', fontsize=18)
        npanel += 1
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def H1ZEUS(path=None, fmt='png', **kwargs):
    """Makes plot of H1 DVCS data with fit lines"""

    subsets = {}
    subsets[1] = [utils.select(data[63], criteria=['Q2 == 8.']),
            utils.select(data[63], criteria=['Q2 == 15.5']), 
            utils.select(data[63], criteria=['Q2 == 25.'])]
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
    for npanel in range(1,5):
        ax = fig.add_subplot(2,2,npanel)
        ax.set_yscale('log')  # y-axis to be logarithmic
        panel(ax, points=subsets[npanel], xaxis=xs[npanel-1], **kwargs)
        if npanel < 3:
            ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.2))
        # y labels
        if npanel==1:
            ax.set_ylabel('$d\\sigma/dt\\quad [{\\rm nb/GeV}^2]$', fontsize=18)
        elif npanel==3:
            ax.set_ylabel('$\\sigma\\quad [{\\rm nb}]$', fontsize=18)
        else: # npanel 4
            ax.set_ylabel('')
        # x labels
        if npanel==1 or npanel==2:
            ax.set_xlabel('$t\\quad [{\\rm GeV}^2]$', fontsize=18)
        elif npanel==3:
            ax.set_xlabel('$W\\quad [{\\rm GeV}]$', fontsize=18)
        else: # npanel 4
            ax.set_xlabel('$Q^2\\quad [{\\rm GeV}^2]$', fontsize=18)
        if npanel==1:
            ax.text(-0.8, 22, '${\\rm H1}$', fontsize=16)
            ax.text(-0.8, 10, '${\\rm W = 82}\\, {\\rm GeV}$', fontsize=12)
            ax.text(-0.8, 3, '$Q^2\\!= 8,\\, 15.5,\\, 25\\,{\\rm GeV}^2$', fontsize=12)
        if npanel==2:
            ax.text(-0.3, 6, '${\\rm ZEUS}$', fontsize=16)
            ax.text(-0.3, 3.5, '${\\rm W = 104}\\, {\\rm GeV}$', fontsize=12)
            ax.text(-0.3, 2.1, '$Q^2\\!= 3.2\\,{\\rm GeV}^2$', fontsize=12)
        if npanel==3:
            ax.text(50, 28, '${\\rm ZEUS}\\, (idem):$', fontsize=16)
            ax.text(50, 5, '${\\rm H1}\\, (idem):$', fontsize=16)
        else: # npanel==4
            ax.text(30, 3, '${\\rm ZEUS}\\, (idem)$', fontsize=16)
            #ax.set_xlim(0, 80)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def COMPASSt(path=None, fmt='png', **kwargs):
    """Makes plot of COMPASS BCS asymmetry in terms of t for various values of xB"""

    title = ''#Prediction for COMPASS BCSA'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    #kinpoints = [(0.007, 1.5), (0.014, 2.5), (0.024, 3.7), 
    #                       (0.045, 3.0), (0.1, 3.0), (0.2, 4.4)]
    kinpoints = [(0.01, 1.4), (0.015, 2.4), (0.03, 5.0), 
                 (0.04, 1.4), (0.07, 2.4), (0.13, 5.0)]
    pn = 1
    ks = {'exptype' : 'fixed target'}
    ks['in1particle'] = 'e'
    ks['in1charge'] = 1
    ks['in1energy'] = 160.
    ks['in1polarization'] = -0.8
    ks['s'] = 2 * Mp * ks['in1energy'] + Mp2
    ks['phi'] = 3.141
    ks['units'] = {'phi' : 'radian'}
    ks['frame'] = 'BMK'
    ks['yaxis'] = 'BCSA'
    for xB, Q2 in kinpoints:
        ax = fig.add_subplot(2,3,pn)
        ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
        #linestyles = ['r-', 'g--', 'b-.', 'p:']
        #labels = ['HERMES+CLAS', 'HERMES+CLAS+HALLA', '+HALLA(phi)']
        tms = np.arange(0.05, 0.6, 0.02)
        ks.update({'xB':xB, 'Q2':Q2})
        panel(ax, xaxis='tm', xs=tms, kins=ks, **kwargs)
        ax.set_xlim(0.0, 0.7)
        ax.set_ylim(-0.1, 0.45)
        # axes labels
        ax.set_xlabel('$-t$')
        ax.text(0.03, 0.35, "%s = %s" % (toTeX['xB'], xB), fontsize=14)
        ax.text(0.03, 0.3, "%s = %s" % (toTeX['Q2'], Q2), fontsize=14)
        if pn == 1:
            #ax.legend(bbox_to_anchor=(0.9, 0.7), loc='upper right',
            #        borderaxespad=0.).draw_frame(0)
            pass
        if (pn != 1) and (pn != 4):
            # no y tick labels on right panels
            ax.set_yticklabels([])
        else:
            # y label on left panels
            ax.set_ylabel(toTeX['BCSA'], fontsize=20)
        if pn <= 3:
            # no x tick labels on upper panels
            ax.set_xticklabels([])
        if pn >= 4:
            # x-label only on lower panels
            ax.set_xlabel(toTeX['tm'], fontsize=17)
        pn += 1
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.2))  # tickmarks
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontsize(14)
    fig.subplots_adjust(bottom=0.12, wspace=0.0, hspace=0.0)
    #ax.set_ylabel('BCS Asymmetry')
    #ax.legend()
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

### FIXME: Following plots are done without using panel(), so they cannot plot bands

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

def EIC(lines=[], path=None, fmt='png'):
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
    labels = [r'$\sigma^{\uparrow}$ ' + t.name for t in lines]
    pn = 0
    for approach in lines:
        approach.__class__.to_conventions(pt)
        approach.__class__.prepare(pt)
        line = approach.Xunp(pt, vars={'phi':phi})
        ax.plot(phi, line, linestyles[pn], linewidth=1, label=labels[pn]) 
        pn += 1
    labels = ['$(\\sigma^{\\uparrow} + \\sigma^{\\downarrow})/2$ ' + t.name for t in lines]
    pn = 0
    for approach in lines:
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

def EICt(path=None, fmt='png', **kwargs):
    """Plot simulated EIC DVCS data."""
    id = 1002
    title = 'EIC (HERA-like simulated)'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    # we have 2x3 points
    ax = fig.add_subplot(1, 1, 1)
    #ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))  # tickmarks
    ax.set_yscale('log')  # y-axis to be logarithmic
    panel(ax, points=data[id], xaxis='tm', **kwargs)
    ax.set_ylabel('$d\\sigma/dt\\quad [{\\rm nb/GeV}^2]$', fontsize=18)
    ax.set_xlabel('$-t\\quad [{\\rm GeV}^2]$', fontsize=18)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def EICTTSA(lines=[], path=None, fmt='png'):
    """Plot EIC TTSA.
    
    FIXME: kinematic completion, charge etc. must be explicit here
    """

    title = 'EIC transversal target spin assymetry'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    fig.subplots_adjust(bottom=0.35, right=0.6)
    # Asymmetry panel
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
    pt.s = 2 * pt.in1energy * (pt.in2energy + math.sqrt(
        pt.in2energy**2 - Mp2)) + Mp2
    #pt.xB = 5.145e-4
    pt.xB = 8.2e-4
    pt.t = -0.275
    pt.Q2 = 4.4
    pt.varphi = -np.pi/2
    pt.units = {'phi' : 'radian'}
    ax = fig.add_subplot(1,1,1)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    phi = np.linspace(0., 2*np.pi, 50)
    utils.fill_kinematics(pt)
    linestyles = ['r-', 'g--', 'b-.', 'p:']
    dasheslengths = [(None, None), (20,5), (20,5,5,5), (5,5)]
    #labels = [r'$\sigma^{\uparrow}$ ' + t.name for t in lines]
    pn = 0
    for approach in lines:
        approach.__class__.to_conventions(pt)
        approach.__class__.prepare(pt)
        mem = approach.m.parameters['KAPS']
        for kaps in [-1.5, 0, 1.5]:
            approach.m.parameters['KAPS'] = kaps
            # Must go to BKM explicitely here
            line = approach.TSA(pt, vars={'phi':np.pi-phi})
            ax.plot(phi, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2, label='$\\kappa_{S}$ = % 4.1f' % kaps) 
            pn += 1
        approach.m.parameters['KAPS'] = mem

    #labels = ['$(\\sigma^{\\uparrow} + \\sigma^{\\downarrow})/2$ ' + t.name for t in lines]
    ax.set_xlim(0.0, 2*np.pi)
    # axes labels
    ax.set_xlabel('$\\phi\\quad {\\rm [rad]}$', fontsize=20)
    ax.set_ylabel('$A_{\\rm UT}^{\\sin(\\phi - \\phi_S)}$', fontsize=20)
    # Legend
    leg = ax.legend(loc='upper center', handlelength=4.0, fancybox=True)
    frame  = leg.get_frame()
    frame.set_facecolor('0.90')    # set the frame face color to light gray
    for t in leg.get_texts():
        t.set_fontsize(18)    # the legend text fontsize
    for l in leg.get_lines():
        l.set_linewidth(2.0)  # the legend line width
    #ax.text(0.2, -0.38, "EIC", fontsize=18)
    ax.text(0.1, -0.37, "$E_e = %d \\,{\\rm GeV}$" % pt.in1energy, fontsize=18)
    ax.text(0.1, -0.46, "$E_p = %d \\,{\\rm GeV}$" % pt.in2energy, fontsize=18)
    ax.text(0.1, -0.55, "$x_B = %s$" % pt.xB, fontsize=18)
    ax.text(0.1, -0.64, "$Q^2 = %s \\,{\\rm GeV}^2$" % pt.Q2, fontsize=18)
    ax.text(0.1, -0.73, "$t = %s \\,{\\rm GeV}^2$" % pt.t, fontsize=18)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def H(theories=[], path=None, fmt='png'):
    """Makes plot of x*Im(cffH)/pi = x * H
    
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

def EICspinALU(th, path=None, fmt='png'):
    """Plot EIC BSA (A_LU).  """
    title = ''
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    fig.subplots_adjust(bottom=0.6, right=0.9, wspace=0.)
    # Asymmetry panel
    pt = Data.DummyPoint()
    pt.exptype = 'collider'
    pt.in1particle = 'e'
    pt.in1charge = -1
    pt.in1energy = 5.
    #pt.in1energy = 20.
    pt.in1polarizationvector = 'L'
    pt.in1polarization = 1
    pt.in2particle = 'p'
    pt.in2energy = 100.
    #pt.in2energy = 250.
    pt.s = 2 * pt.in1energy * (pt.in2energy + math.sqrt(
        pt.in2energy**2 - Mp2)) + Mp2
    pt.xB = 3.2e-3
    pt.t = -0.25
    pt.Q2 = 4.4
    pt.units = {'phi' : 'radian'}
    ########   TSA  ###################
    #
    ###  phi panel
    ax = fig.add_subplot(1,3,1)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    phi = np.linspace(0., 2*np.pi, 50)
    utils.fill_kinematics(pt)
    linestyles = ['r--', 'r--', 'g-', 'b-.', 'p:']
    dasheslengths = [(None, None), (20,5), (20,5,5,5), (5,5)]
    #labels = [r'$\sigma^{\uparrow}$ ' + t.name for t in lines]
    pn = 0
    th.__class__.to_conventions(pt)
    th.__class__.prepare(pt)
    # Must go to BKM explicitely here
    line = th._BSA(pt, vars={'phi':np.pi-phi})
    ax.plot(phi, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2, label=th.name) 
    ax.set_xlim(0.0, 2*np.pi)
    ax.set_ylim(-0.65, 0.65)
    # axes labels
    ax.set_xlabel('$\\phi\\quad {\\rm [rad]}$', fontsize=20)
    ax.set_ylabel('$A_{\\rm LU}$', fontsize=20)
    ax.text(1, 0.35, "$t = %s \\,{\\rm GeV}^2$" % pt.t, fontsize=14)
    ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2])
    ax.set_xticklabels(['$0$', '$\\pi/2$', '$\\pi$', '$3 \\pi/2$'])
    # Legend
    leg = ax.legend(loc='upper center', handlelength=4.0, fancybox=True)
    frame  = leg.get_frame()
    frame.set_facecolor('0.90')    # set the frame face color to light gray
    for t in leg.get_texts():
        t.set_fontsize(18)    # the legend text fontsize
    for l in leg.get_lines():
        l.set_linewidth(2.0)  # the legend line width
    #
    ###  t panel
    pt.FTn = -1
    ax = fig.add_subplot(1,3,2)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    tms = list(np.linspace(0.,0.1))+list(np.linspace(0.1, 0.85))
    line = []
    for tm in tms:
        pt.t = -tm
        line.append(th.BSA(pt))
    ax.plot(tms, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2) 
    ax.set_xlim(0., 0.9)
    ax.set_ylim(-0.65, 0.65)
    ax.set_xlabel('$-t \\,{\\rm [GeV}^2{\\rm ]}$', fontsize=20)
    ax.set_yticklabels([])
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.2))  # tickmarks
    ax.text(0.1, -0.17, "$E_e = %d \\,{\\rm GeV}$" % pt.in1energy, fontsize=14)
    ax.text(0.1, -0.26, "$E_p = %d \\,{\\rm GeV}$" % pt.in2energy, fontsize=14)
    ax.text(0.1, -0.35, "$x_B = %s$" % pt.xB, fontsize=14)
    ax.text(0.1, -0.44, "$Q^2 = %s \\,{\\rm GeV}^2$" % pt.Q2, fontsize=14)
    #
    ###  xB panel
    pt.t = -0.25
    pt.tm = 0.25
    ax = fig.add_subplot(1,3,3)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    xBs = np.logspace(-4., -1., 40)
    line = []
    pt.Q2 = 2.5
    for xB in xBs:
        del pt.W
        pt.xB = xB
        utils._complete_xBWQ2(pt)
        line.append(th.BSA(pt))
    ax.semilogx(xBs, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2,
            label="$Q^2= %s \\,{\\rm GeV}^2$" % pt.Q2) 
    pn += 1
    line = []
    pt.Q2 = 13.9
    for xB in xBs:
        del pt.W
        pt.xB = xB
        utils._complete_xBWQ2(pt)
        line.append(th.BSA(pt))
    ax.semilogx(xBs, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2,
            label="$Q^2= %s \\,{\\rm GeV}^2$" % pt.Q2) 
    ax.set_xlim(5.e-5, 0.2)
    ax.set_ylim(-0.65, 0.65)
    ax.set_xlabel('$x_{\\rm B}$', fontsize=20)
    ax.set_yticklabels([])
    ax.text(0.001, 0.35, "$t = %s \\,{\\rm GeV}^2$" % pt.t, fontsize=14)
    # Legend
    leg = ax.legend(loc='lower center', handlelength=4.0, fancybox=True)
    frame  = leg.get_frame()
    frame.set_facecolor('0.90')    # set the frame face color to light gray
    for t in leg.get_texts():
        t.set_fontsize(14)    # the legend text fontsize
    for l in leg.get_lines():
        l.set_linewidth(2.0)  # the legend line width
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def EICspinAUT(th, path=None, fmt='png'):
    """Plot EIC TTSA.  """
    title = ''
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    fig.subplots_adjust(bottom=0.6, right=0.9, wspace=0.)
    # Asymmetry panel
    pt = Data.DummyPoint()
    pt.exptype = 'collider'
    pt.in1particle = 'e'
    pt.in1charge = -1
    pt.in1energy = 5.
    #pt.in1energy = 20.
    pt.in1polarization = 0.0
    pt.in2particle = 'p'
    pt.in2energy = 100.
    #pt.in2energy = 250.
    pt.in2polarizationvector = 'T'
    pt.in2polarization = 1
    pt.s = 2 * pt.in1energy * (pt.in2energy + math.sqrt(
        pt.in2energy**2 - Mp2)) + Mp2
    pt.xB = 8.155e-3
    pt.t = -0.25
    pt.Q2 = 4.4
    pt.varphi = -np.pi/2
    pt.units = {'phi' : 'radian'}
    ########   TSA  ###################
    #
    ###  phi panel
    ax = fig.add_subplot(1,3,1)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    phi = np.linspace(0., 2*np.pi, 50)
    utils.fill_kinematics(pt)
    linestyles = ['r--', 'r--', 'g-', 'b-.', 'p:']
    dasheslengths = [(None, None), (20,5), (20,5,5,5), (5,5)]
    #labels = [r'$\sigma^{\uparrow}$ ' + t.name for t in lines]
    pn = 0
    th.__class__.to_conventions(pt)
    th.__class__.prepare(pt)
    # Must go to BKM explicitely here
    line = th.TSA(pt, vars={'phi':np.pi-phi})
    ax.plot(phi, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2, label=th.name) 
    ax.set_xlim(0.0, 2*np.pi)
    ax.set_ylim(-0.65, 0.65)
    # axes labels
    ax.set_xlabel('$\\phi\\quad {\\rm [rad]}$', fontsize=20)
    ax.set_ylabel('$A_{\\rm UT}^{\\sin(\\phi - \\phi_S)}$', fontsize=20)
    ax.text(1, 0.35, "$t = %s \\,{\\rm GeV}^2$" % pt.t, fontsize=14)
    ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2])
    ax.set_xticklabels(['$0$', '$\\pi/2$', '$\\pi$', '$3 \\pi/2$'])
    # Legend
    leg = ax.legend(loc='upper center', handlelength=4.0, fancybox=True)
    frame  = leg.get_frame()
    frame.set_facecolor('0.90')    # set the frame face color to light gray
    for t in leg.get_texts():
        t.set_fontsize(18)    # the legend text fontsize
    for l in leg.get_lines():
        l.set_linewidth(2.0)  # the legend line width
    #
    ###  t panel
    pt.FTn = 1
    ax = fig.add_subplot(1,3,2)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    tms = list(np.linspace(0.,0.1))+list(np.linspace(0.1, 0.85))
    line = []
    for tm in tms:
        pt.t = -tm
        # Must change sign to go to Trento:
        line.append(-th.TSA(pt))
    ax.plot(tms, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2) 
    ax.set_xlim(0., 0.9)
    ax.set_ylim(-0.65, 0.65)
    ax.set_xlabel('$-t \\,{\\rm [GeV}^2{\\rm ]}$', fontsize=20)
    ax.set_yticklabels([])
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.2))  # tickmarks
    ax.text(0.1, -0.17, "$E_e = %d \\,{\\rm GeV}$" % pt.in1energy, fontsize=14)
    ax.text(0.1, -0.26, "$E_p = %d \\,{\\rm GeV}$" % pt.in2energy, fontsize=14)
    ax.text(0.1, -0.35, "$x_B = %s$" % pt.xB, fontsize=14)
    ax.text(0.1, -0.44, "$Q^2 = %s \\,{\\rm GeV}^2$" % pt.Q2, fontsize=14)
    #
    ###  xB panel
    pt.t = -0.25
    pt.tm = 0.25
    ax = fig.add_subplot(1,3,3)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    xBs = np.logspace(-4., -1., 40)
    line = []
    pt.Q2 = 2.5
    for xB in xBs:
        del pt.W
        pt.xB = xB
        utils._complete_xBWQ2(pt)
        # Must change sign to go to Trento:
        line.append(-th.TSA(pt))
    ax.semilogx(xBs, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2,
            label="$Q^2 = %s \\,{\\rm GeV}^2$" % pt.Q2) 
    pn += 1
    line = []
    pt.Q2 = 13.9
    for xB in xBs:
        del pt.W
        pt.xB = xB
        utils._complete_xBWQ2(pt)
        line.append(-th.TSA(pt))
    ax.semilogx(xBs, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2,
            label="$Q^2 = %s \\,{\\rm GeV}^2$" % pt.Q2) 
    ax.set_xlim(5.e-5, 0.2)
    ax.set_ylim(-0.65, 0.65)
    ax.set_xlabel('$x_{\\rm B}$', fontsize=20)
    ax.set_yticklabels([])
    ax.text(0.001, 0.35, "$t = %s \\,{\\rm GeV}^2$" % pt.t, fontsize=14)
    # Legend
    leg = ax.legend(loc='lower center', handlelength=4.0, fancybox=True)
    frame  = leg.get_frame()
    frame.set_facecolor('0.90')    # set the frame face color to light gray
    for t in leg.get_texts():
        t.set_fontsize(14)    # the legend text fontsize
    for l in leg.get_lines():
        l.set_linewidth(2.0)  # the legend line width
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def EICspinAUL(th, path=None, fmt='png'):
    """Plot EIC LTSA.  """
    title = ''
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    fig.subplots_adjust(bottom=0.6, right=0.9, wspace=0.)
    # Asymmetry panel
    pt = Data.DummyPoint()
    pt.exptype = 'collider'
    pt.in1particle = 'e'
    pt.in1charge = -1
    pt.in1energy = 5.
    #pt.in1energy = 20.
    pt.in1polarization = 0.0
    pt.in2particle = 'p'
    pt.in2energy = 100.
    #pt.in2energy = 250.
    pt.in2polarizationvector = 'L'
    pt.in2polarization = 1
    pt.s = 2 * pt.in1energy * (pt.in2energy + math.sqrt(
        pt.in2energy**2 - Mp2)) + Mp2
    pt.xB = 1.3e-2
    #pt.xB = 8.155e-3
    pt.t = -0.25
    pt.Q2 = 4.4
    pt.units = {'phi' : 'radian'}
    ########   TSA  ###################
    #
    ###  phi panel
    ax = fig.add_subplot(1,3,1)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    phi = np.linspace(0., 2*np.pi, 50)
    utils.fill_kinematics(pt)
    linestyles = ['r--', 'r--', 'g-', 'b-.', 'p:']
    dasheslengths = [(None, None), (20,5), (20,5,5,5), (5,5)]
    #labels = [r'$\sigma^{\uparrow}$ ' + t.name for t in lines]
    pn = 0
    th.__class__.to_conventions(pt)
    th.__class__.prepare(pt)
    # Must go to BKM explicitely here
    line = th.TSA(pt, vars={'phi':np.pi-phi})
    ax.plot(phi, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2, label=th.name) 
    ax.set_xlim(0.0, 2*np.pi)
    ax.set_ylim(-0.65, 0.65)
    # axes labels
    ax.set_xlabel('$\\phi\\quad {\\rm [rad]}$', fontsize=20)
    ax.set_ylabel('$A_{\\rm UL}$', fontsize=20)
    ax.text(1, 0.35, "$t = %s \\,{\\rm GeV}^2$" % pt.t, fontsize=14)
    ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2])
    ax.set_xticklabels(['$0$', '$\\pi/2$', '$\\pi$', '$3 \\pi/2$'])
    # Legend
    leg = ax.legend(loc='upper center', handlelength=4.0, fancybox=True)
    frame  = leg.get_frame()
    frame.set_facecolor('0.90')    # set the frame face color to light gray
    for t in leg.get_texts():
        t.set_fontsize(18)    # the legend text fontsize
    for l in leg.get_lines():
        l.set_linewidth(2.0)  # the legend line width
    #
    ###  t panel
    pt.FTn = -1
    ax = fig.add_subplot(1,3,2)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    tms = list(np.linspace(0.,0.1))+list(np.linspace(0.1, 0.85))
    line = []
    for tm in tms:
        pt.t = -tm
        line.append(th.TSA(pt))
    ax.plot(tms, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2) 
    ax.set_xlim(0., 0.9)
    ax.set_ylim(-0.65, 0.65)
    ax.set_xlabel('$-t \\,{\\rm [GeV}^2{\\rm ]}$', fontsize=20)
    ax.set_yticklabels([])
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.2))  # tickmarks
    ax.text(0.1, -0.17, "$E_e = %d \\,{\\rm GeV}$" % pt.in1energy, fontsize=14)
    ax.text(0.1, -0.26, "$E_p = %d \\,{\\rm GeV}$" % pt.in2energy, fontsize=14)
    ax.text(0.1, -0.35, "$x_B = %s$" % pt.xB, fontsize=14)
    ax.text(0.1, -0.44, "$Q^2 = %s \\,{\\rm GeV}^2$" % pt.Q2, fontsize=14)
    #
    ###  xB panel
    pt.t = -0.25
    pt.tm = 0.25
    ax = fig.add_subplot(1,3,3)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    xBs = np.logspace(-4., -1., 40)
    line = []
    pt.Q2 = 2.5
    for xB in xBs:
        del pt.W
        pt.xB = xB
        utils._complete_xBWQ2(pt)
        line.append(th.TSA(pt))
    ax.semilogx(xBs, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2,
            label="$Q^2 = %s \\,{\\rm GeV}^2$" % pt.Q2) 
    pn += 1
    line = []
    pt.Q2 = 13.9
    for xB in xBs:
        del pt.W
        pt.xB = xB
        utils._complete_xBWQ2(pt)
        line.append(th.TSA(pt))
    ax.semilogx(xBs, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2,
            label="$Q^2 = %s \\,{\\rm GeV}^2$" % pt.Q2) 
    ax.set_xlim(5.e-5, 0.2)
    ax.set_ylim(-0.65, 0.65)
    ax.set_xlabel('$x_{\\rm B}$', fontsize=20)
    ax.set_yticklabels([])
    ax.text(0.001, 0.35, "$t = %s \\,{\\rm GeV}^2$" % pt.t, fontsize=14)
    # Legend
    leg = ax.legend(loc='lower center', handlelength=4.0, fancybox=True)
    frame  = leg.get_frame()
    frame.set_facecolor('0.90')    # set the frame face color to light gray
    for t in leg.get_texts():
        t.set_fontsize(14)    # the legend text fontsize
    for l in leg.get_lines():
        l.set_linewidth(2.0)  # the legend line width
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def EICcharge(th, path=None, fmt='png'):
    """Plot EIC BCA (A_C).  """
    title = ''
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    fig.subplots_adjust(bottom=0.6, right=0.9, wspace=0.)
    # Asymmetry panel
    pt = Data.DummyPoint()
    pt.exptype = 'collider'
    pt.in1particle = 'e'
    pt.in1charge = 1
    pt.in1energy = 5.
    #pt.in1energy = 20.
    pt.in1polarization = 0
    pt.in2particle = 'p'
    pt.in2energy = 100.
    #pt.in2energy = 250.
    pt.s = 2 * pt.in1energy * (pt.in2energy + math.sqrt(
        pt.in2energy**2 - Mp2)) + Mp2
    pt.xB = 5.1e-3
    #pt.xB = 5.1e-4
    pt.t = -0.25
    pt.Q2 = 4.4
    pt.units = {'phi' : 'radian'}
    ########   TSA  ###################
    #
    ###  phi panel
    ax = fig.add_subplot(1,3,1)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    phi = np.linspace(0., 2*np.pi, 50)
    utils.fill_kinematics(pt)
    linestyles = ['r--', 'r--', 'g-', 'b-.', 'p:']
    dasheslengths = [(None, None), (20,5), (20,5,5,5), (5,5)]
    #labels = [r'$\sigma^{\uparrow}$ ' + t.name for t in lines]
    pn = 0
    th.__class__.to_conventions(pt)
    th.__class__.prepare(pt)
    # Must go to BKM explicitely here
    line = th.BCA(pt, vars={'phi':np.pi-phi})
    ax.plot(phi, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2, label=th.name) 
    ax.set_xlim(0.0, 2*np.pi)
    ax.set_ylim(-0.45, 0.45)
    # axes labels
    ax.set_xlabel('$\\phi\\quad {\\rm [rad]}$', fontsize=20)
    ax.set_ylabel('$A_{\\rm C}$', fontsize=20)
    ax.text(1, 0.35, "$t = %s \\,{\\rm GeV}^2$" % pt.t, fontsize=14)
    ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2])
    ax.set_xticklabels(['$0$', '$\\pi/2$', '$\\pi$', '$3 \\pi/2$'])
    # Legend
    leg = ax.legend(loc='upper center', handlelength=4.0, fancybox=True)
    frame  = leg.get_frame()
    frame.set_facecolor('0.90')    # set the frame face color to light gray
    for t in leg.get_texts():
        t.set_fontsize(18)    # the legend text fontsize
    for l in leg.get_lines():
        l.set_linewidth(2.0)  # the legend line width
    #
    ###  t panel
    pt.FTn = 1
    ax = fig.add_subplot(1,3,2)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    tms = list(np.linspace(0.,0.1))+list(np.linspace(0.1, 0.85))
    line = []
    for tm in tms:
        pt.t = -tm
        # Change sign to go to Trento
        line.append(-th.BCA(pt))
    ax.plot(tms, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2) 
    ax.set_xlim(0., 0.9)
    ax.set_ylim(-0.45, 0.45)
    ax.set_xlabel('$-t \\,{\\rm [GeV}^2{\\rm ]}$', fontsize=20)
    ax.set_yticklabels([])
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.2))  # tickmarks
    ax.text(0.1, -0.17, "$E_e = %d \\,{\\rm GeV}$" % pt.in1energy, fontsize=14)
    ax.text(0.1, -0.26, "$E_p = %d \\,{\\rm GeV}$" % pt.in2energy, fontsize=14)
    ax.text(0.1, -0.35, "$x_B = %s$" % pt.xB, fontsize=14)
    ax.text(0.1, -0.44, "$Q^2 = %s \\,{\\rm GeV}^2$" % pt.Q2, fontsize=14)
    #
    ###  xB panel
    pt.t = -0.25
    pt.tm = 0.25
    ax = fig.add_subplot(1,3,3)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    xBs = np.logspace(-4., -0.7, 40)
    line = []
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
    ax.set_xlabel('$x_{\\rm B}$', fontsize=20)
    ax.set_yticklabels([])
    ax.text(0.001, 0.35, "$t = %s \\,{\\rm GeV}^2$" % pt.t, fontsize=14)
    # Legend
    leg = ax.legend(loc='lower center', handlelength=4.0, fancybox=True)
    frame  = leg.get_frame()
    frame.set_facecolor('0.90')    # set the frame face color to light gray
    for t in leg.get_texts():
        t.set_fontsize(14)    # the legend text fontsize
    for l in leg.get_lines():
        l.set_linewidth(2.0)  # the legend line width
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def EICphase1(th, path=None, fmt='png'):
    """Plot various EIC 5x100 asymmetries.  """
    title = ''
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    fig.subplots_adjust(wspace=0.,hspace=0.)
    # Asymmetry panel
    pt = Data.DummyPoint()
    pt.exptype = 'collider'
    pt.in1particle = 'e'
    pt.in1charge = -1
    pt.in1energy = 5.
    pt.in1polarizationvector = 'L'
    pt.in1polarization = 1
    pt.in2particle = 'p'
    pt.in2energy = 100.
    pt.s = 2 * pt.in1energy * (pt.in2energy + math.sqrt(
        pt.in2energy**2 - Mp2)) + Mp2
    pt.xB = 3.2e-3
    pt.t = -0.25
    pt.Q2 = 4.4
    pt.units = {'phi' : 'radian'}
    ########   TSA  ###################
    #
    ###  phi panel
    ax = fig.add_subplot(4,3,1)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    phi = np.linspace(0., 2*np.pi, 50)
    utils.fill_kinematics(pt)
    linestyles = ['r--', 'r--', 'g-', 'b-.', 'p:']
    dasheslengths = [(None, None), (20,5), (20,5,5,5), (5,5)]
    #labels = [r'$\sigma^{\uparrow}$ ' + t.name for t in lines]
    pn = 0
    th.__class__.to_conventions(pt)
    th.__class__.prepare(pt)
    # Must go to BKM explicitely here
    line = th._BSA(pt, vars={'phi':np.pi-phi})
    ax.plot(phi, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2, label=th.name) 
    ax.set_xlim(0.0, 2*np.pi)
    ax.set_ylim(-0.65, 0.65)
    # axes labels
    ax.set_ylabel('$A_{\\rm LU}$', fontsize=20)
    ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2])
    ax.set_xticklabels(['$0$', '$\\pi/2$', '$\\pi$', '$3 \\pi/2$'])
    ax.text(0.1, -0.17, "$E_e = %d \\,{\\rm GeV}$" % pt.in1energy, fontsize=14)
    ax.text(0.1, -0.35, "$E_p = %d \\,{\\rm GeV}$" % pt.in2energy, fontsize=14)
    ax.text(0.1, -0.55, "$Q^2 = %s \\,{\\rm GeV}^2$" % pt.Q2, fontsize=14)
    # Legend
    leg = ax.legend(loc='upper center', handlelength=4.0, fancybox=True)
    frame  = leg.get_frame()
    frame.set_facecolor('0.90')    # set the frame face color to light gray
    for t in leg.get_texts():
        t.set_fontsize(18)    # the legend text fontsize
    for l in leg.get_lines():
        l.set_linewidth(2.0)  # the legend line width
    #
    ###  t panel
    pt.FTn = -1
    ax = fig.add_subplot(4,3,2)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    tms = list(np.linspace(0.,0.1))+list(np.linspace(0.1, 0.85))
    line = []
    for tm in tms:
        pt.t = -tm
        line.append(th.BSA(pt))
    ax.plot(tms, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2) 
    ax.set_xlim(0., 0.9)
    ax.set_ylim(-0.65, 0.65)
    ax.set_yticklabels([])
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.2))  # tickmarks
    ax.text(0.1, -0.35, "$x_B = %s$" % pt.xB, fontsize=14)
    #
    ###  xB panel
    pt.t = -0.25
    pt.tm = 0.25
    ax = fig.add_subplot(4,3,3)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    xBs = np.logspace(-4., -1., 40)
    line = []
    pt.Q2 = 2.5
    for xB in xBs:
        del pt.W
        pt.xB = xB
        utils._complete_xBWQ2(pt)
        line.append(th.BSA(pt))
    ax.semilogx(xBs, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2,
            label="$Q^2= %s \\,{\\rm GeV}^2$" % pt.Q2) 
    pn += 1
    line = []
    pt.Q2 = 13.9
    for xB in xBs:
        del pt.W
        pt.xB = xB
        utils._complete_xBWQ2(pt)
        line.append(th.BSA(pt))
    ax.semilogx(xBs, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2,
            label="$Q^2= %s \\,{\\rm GeV}^2$" % pt.Q2) 
    ax.set_xlim(5.e-5, 0.2)
    ax.set_ylim(-0.65, 0.65)
    ax.set_yticklabels([])
    # Legend
    leg = ax.legend(loc='lower center', handlelength=4.0, fancybox=True)
    frame  = leg.get_frame()
    frame.set_facecolor('0.90')    # set the frame face color to light gray
    for t in leg.get_texts():
        t.set_fontsize(10)    # the legend text fontsize
    for l in leg.get_lines():
        l.set_linewidth(2.0)  # the legend line width
    #
    ####    AUT
    #
    # Asymmetry panel
    pt = Data.DummyPoint()
    pt.exptype = 'collider'
    pt.in1particle = 'e'
    pt.in1charge = -1
    pt.in1energy = 5.
    pt.in1polarization = 0.0
    pt.in2particle = 'p'
    pt.in2energy = 100.
    pt.in2polarizationvector = 'T'
    pt.in2polarization = 1
    pt.s = 2 * pt.in1energy * (pt.in2energy + math.sqrt(
        pt.in2energy**2 - Mp2)) + Mp2
    pt.xB = 8.155e-3
    pt.t = -0.25
    pt.Q2 = 4.4
    pt.varphi = -np.pi/2
    pt.units = {'phi' : 'radian'}
    ########   TSA  ###################
    #
    ###  phi panel
    ax = fig.add_subplot(4,3,4)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    phi = np.linspace(0., 2*np.pi, 50)
    utils.fill_kinematics(pt)
    linestyles = ['r--', 'r--', 'g-', 'b-.', 'p:']
    dasheslengths = [(None, None), (20,5), (20,5,5,5), (5,5)]
    #labels = [r'$\sigma^{\uparrow}$ ' + t.name for t in lines]
    pn = 0
    th.__class__.to_conventions(pt)
    th.__class__.prepare(pt)
    # Must go to BKM explicitely here
    line = th.TSA(pt, vars={'phi':np.pi-phi})
    ax.plot(phi, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2, label=th.name) 
    ax.set_xlim(0.0, 2*np.pi)
    ax.set_ylim(-0.65, 0.65)
    # axes labels
    ax.set_ylabel('$A_{\\rm UT}^{\\sin(\\phi - \\phi_S)}$', fontsize=20)
    ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2])
    ax.set_xticklabels(['$0$', '$\\pi/2$', '$\\pi$', '$3 \\pi/2$'])
    #
    ###  t panel
    pt.FTn = 1
    ax = fig.add_subplot(4,3,5)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    tms = list(np.linspace(0.,0.1))+list(np.linspace(0.1, 0.85))
    line = []
    for tm in tms:
        pt.t = -tm
        # Must change sign to go to Trento:
        line.append(-th.TSA(pt))
    ax.plot(tms, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2) 
    ax.set_xlim(0., 0.9)
    ax.set_ylim(-0.65, 0.65)
    ax.set_yticklabels([])
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.2))  # tickmarks
    ax.text(0.1, -0.35, "$x_B = %s$" % pt.xB, fontsize=14)
    #
    ###  xB panel
    pt.t = -0.25
    pt.tm = 0.25
    ax = fig.add_subplot(4,3,6)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    xBs = np.logspace(-4., -1., 40)
    line = []
    pt.Q2 = 2.5
    for xB in xBs:
        del pt.W
        pt.xB = xB
        utils._complete_xBWQ2(pt)
        # Must change sign to go to Trento:
        line.append(-th.TSA(pt))
    ax.semilogx(xBs, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2,
            label="$Q^2 = %s \\,{\\rm GeV}^2$" % pt.Q2) 
    pn += 1
    line = []
    pt.Q2 = 13.9
    for xB in xBs:
        del pt.W
        pt.xB = xB
        utils._complete_xBWQ2(pt)
        line.append(-th.TSA(pt))
    ax.semilogx(xBs, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2,
            label="$Q^2 = %s \\,{\\rm GeV}^2$" % pt.Q2) 
    ax.set_xlim(5.e-5, 0.2)
    ax.set_ylim(-0.65, 0.65)
    ax.set_yticklabels([])
    #
    #   A_UL
    # 
    pt = Data.DummyPoint()
    pt.exptype = 'collider'
    pt.in1particle = 'e'
    pt.in1charge = -1
    pt.in1energy = 5.
    pt.in1polarization = 0.0
    pt.in2particle = 'p'
    pt.in2energy = 100.
    pt.in2polarizationvector = 'L'
    pt.in2polarization = 1
    pt.s = 2 * pt.in1energy * (pt.in2energy + math.sqrt(
        pt.in2energy**2 - Mp2)) + Mp2
    pt.xB = 1.3e-2
    pt.t = -0.25
    pt.Q2 = 4.4
    pt.units = {'phi' : 'radian'}
    ########   TSA  ###################
    #
    ###  phi panel
    ax = fig.add_subplot(4,3,7)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    phi = np.linspace(0., 2*np.pi, 50)
    utils.fill_kinematics(pt)
    linestyles = ['r--', 'r--', 'g-', 'b-.', 'p:']
    dasheslengths = [(None, None), (20,5), (20,5,5,5), (5,5)]
    #labels = [r'$\sigma^{\uparrow}$ ' + t.name for t in lines]
    pn = 0
    th.__class__.to_conventions(pt)
    th.__class__.prepare(pt)
    # Must go to BKM explicitely here
    line = th.TSA(pt, vars={'phi':np.pi-phi})
    ax.plot(phi, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2, label=th.name) 
    ax.set_xlim(0.0, 2*np.pi)
    ax.set_ylim(-0.65, 0.65)
    # axes labels
    ax.set_ylabel('$A_{\\rm UL}$', fontsize=20)
    ax.text(1, 0.35, "$t = %s \\,{\\rm GeV}^2$" % pt.t, fontsize=14)
    ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2])
    ax.set_xticklabels(['$0$', '$\\pi/2$', '$\\pi$', '$3 \\pi/2$'])
    ###  t panel
    pt.FTn = -1
    ax = fig.add_subplot(4,3,8)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    tms = list(np.linspace(0.,0.1))+list(np.linspace(0.1, 0.85))
    line = []
    for tm in tms:
        pt.t = -tm
        line.append(th.TSA(pt))
    ax.plot(tms, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2) 
    ax.set_xlim(0., 0.9)
    ax.set_ylim(-0.65, 0.65)
    ax.set_yticklabels([])
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.2))  # tickmarks
    ax.text(0.1, -0.35, "$x_B = %s$" % pt.xB, fontsize=14)
    #
    ###  xB panel
    pt.t = -0.25
    pt.tm = 0.25
    ax = fig.add_subplot(4,3,9)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    xBs = np.logspace(-4., -1., 40)
    line = []
    pt.Q2 = 2.5
    for xB in xBs:
        del pt.W
        pt.xB = xB
        utils._complete_xBWQ2(pt)
        line.append(th.TSA(pt))
    ax.semilogx(xBs, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2,
            label="$Q^2 = %s \\,{\\rm GeV}^2$" % pt.Q2) 
    pn += 1
    line = []
    pt.Q2 = 13.9
    for xB in xBs:
        del pt.W
        pt.xB = xB
        utils._complete_xBWQ2(pt)
        line.append(th.TSA(pt))
    ax.semilogx(xBs, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2,
            label="$Q^2 = %s \\,{\\rm GeV}^2$" % pt.Q2) 
    ax.set_xlim(5.e-5, 0.2)
    ax.set_ylim(-0.65, 0.65)
    ax.set_yticklabels([])
    #
    ###    BCA
    #
    pt = Data.DummyPoint()
    pt.exptype = 'collider'
    pt.in1particle = 'e'
    pt.in1charge = 1
    pt.in1energy = 5.
    pt.in1polarization = 0
    pt.in2particle = 'p'
    pt.in2energy = 100.
    pt.s = 2 * pt.in1energy * (pt.in2energy + math.sqrt(
        pt.in2energy**2 - Mp2)) + Mp2
    pt.xB = 5.1e-3
    pt.t = -0.25
    pt.Q2 = 4.4
    pt.units = {'phi' : 'radian'}
    ###  phi panel
    ax = fig.add_subplot(4,3,10)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    phi = np.linspace(0., 2*np.pi, 50)
    utils.fill_kinematics(pt)
    linestyles = ['r--', 'r--', 'g-', 'b-.', 'p:']
    dasheslengths = [(None, None), (20,5), (20,5,5,5), (5,5)]
    #labels = [r'$\sigma^{\uparrow}$ ' + t.name for t in lines]
    pn = 0
    th.__class__.to_conventions(pt)
    th.__class__.prepare(pt)
    # Must go to BKM explicitely here
    line = th.BCA(pt, vars={'phi':np.pi-phi})
    ax.plot(phi, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2, label=th.name) 
    ax.set_xlim(0.0, 2*np.pi)
    ax.set_ylim(-0.45, 0.45)
    # axes labels
    ax.set_xlabel('$\\phi\\quad {\\rm [rad]}$', fontsize=20)
    ax.set_ylabel('$A_{\\rm C}$', fontsize=20)
    ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2])
    ax.set_xticklabels(['$0$', '$\\pi/2$', '$\\pi$', '$3 \\pi/2$'])
    ###  t panel
    pt.FTn = 1
    ax = fig.add_subplot(4,3,11)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    tms = list(np.linspace(0.,0.1))+list(np.linspace(0.1, 0.85))
    line = []
    for tm in tms:
        pt.t = -tm
        # Change sign to go to Trento
        line.append(-th.BCA(pt))
    ax.plot(tms, line, linestyles[pn], dashes=dasheslengths[pn], linewidth=2) 
    ax.set_xlim(0., 0.9)
    ax.set_ylim(-0.45, 0.45)
    ax.set_xlabel('$-t \\,{\\rm [GeV}^2{\\rm ]}$', fontsize=20)
    ax.set_yticklabels([])
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.2))  # tickmarks
    ax.text(0.1, -0.35, "$x_B = %s$" % pt.xB, fontsize=14)
    #
    ###  xB panel
    pt.t = -0.25
    pt.tm = 0.25
    ax = fig.add_subplot(4,3,12)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    xBs = np.logspace(-4., -0.7, 40)
    line = []
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
    ax.set_xlabel('$x_{\\rm B}$', fontsize=20)
    ax.set_yticklabels([])
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

#################################################################
##                                                             ##
##  [3]  Other plots                                           ##
##                                                             ##
#################################################################

def xBt(points, path=None, fmt='png'):
    """Plot xB-t distribution of data sets. """
    
    title = 'xB-t distribution of data'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    ax = fig.add_subplot(1, 1, 1)
    setshapes = ['o', 's', '^', 'd']  # first circles, then squares ...
    setcolors = ['blue', 'black', 'purple', 'green']  # circles are blue, squares are black, ...
    setn = 0
    if isinstance(points[0], Data.DataPoint): points = [points]
    for set in points:
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

def CFF(cffs=['ImH', 'ReH'], path=None, fmt='png', **kwargs):
    """Makes plots of cffs given by various theories/models
    
    cffs    -- List of CFFs to be plotted. Each produces two panels.

    """
    title = ''
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    colors = ['red', 'brown']     # worm human colors :-)
    nncolors = ['blue', 'green']  # cold computer colors
    linestyles = ['solid', 'dashed']
    # Define abscissas
    #logxvals = np.power(10., np.arange(-3.0, -0.001, 0.1))  # left panel
    logxvals = np.logspace(-5.0, -0.01, 20)  # left panel
    xvals = np.linspace(0.02, 0.15, 20) # right panel
    # ordinates 
    #ylims = {'ImH': (-4.3, 35), 'ReH': (-6.5, 8),
    ylims = {'ImH': (-4.3, 35), 'ReH': (-15, 8),
             'ImE': (-40, 35), 'ReE': (-15, 30),
             'ImEt': (-200, 300), 'ReEt': (-150, 150),
             'ImHt': (-40, 50), 'ReHt': (-15, 30)}
    # Plot panels
    for n in range(len(cffs)):
        cff = cffs[n]
        # all-x logarithmic
        ax = fig.add_subplot(len(cffs), 2, 2*n+1)
        ax.set_xscale('log')  # x-axis to be logarithmic
        panel(ax, xaxis='xi', xs=logxvals, kins={'yaxis':cff, 't':-0.2, 'Q2':4.}, **kwargs)
        ax.set_xlabel(toTeX['xixB'], fontsize=15)
        ax.set_ylabel(toTeX['%s' % cff], fontsize=18)
        ax.axhspan(-0.0005, 0.0005, facecolor='g', alpha=0.6)  # horizontal bar
        #ax.axvspan(0.025, 0.136, facecolor='g', alpha=0.1)  # vertical band
        #ax.text(0.35, 0.95, "data region", transform=ax.transAxes, 
        #        fontsize=14, fontweight='bold', va='top')
        #apply(ax.set_ylim, ylims[cff])
        #ax.set_xlim(0.005, 1.0)
        if n == 0:
            ax.legend(loc='upper right')
            ax.legend().draw_frame(0)
            ax.text(0.006, 0.5, "$t = -0.2\\, {\\rm GeV}^2$",# transform=ax.transAxes, 
                    fontsize=15)
        # measured x linear
        ax = fig.add_subplot(len(cffs), 2, 2*n+2)
        panel(ax, xaxis='xi', xs=xvals, kins={'yaxis':cff, 't':-0.2,  'Q2':4.}, **kwargs)
        ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
        ax.set_xlabel(toTeX['xixB'], fontsize=15)
        #apply(ax.set_ylim, ylims[cff])
        #ax.axvspan(0.025, 0.136, facecolor='g', alpha=0.1)  # vertical band
        #ax.set_xlim(0.02, 0.142)
        #ax.text(0.20, 0.95, "d   a   t   a       r   e   g   i   o   n", 
        #        transform=ax.transAxes, fontsize=14, fontweight='bold', va='top')
        #ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.02))  # tickmarks
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def CFF2(cffs=['ImH', 'ReH'], path=None, fmt='png', **kwargs):
    """Makes plots of cffs at two t given by various theories/models
    
    cffs    -- List of CFFs to be plotted. Each produces two panels.

    """
    title = ''
    fig = plt.figure()
    #fig.canvas.set_window_title(title)
    #fig.suptitle(title)
    colors = ['red', 'brown']     # worm human colors :-)
    nncolors = ['blue', 'green']  # cold computer colors
    linestyles = ['solid', 'dashed']
    # Define abscissas
    logxvals = np.logspace(-3.0, -0.01, 40)
    # ordinates 
    #ylims = {'ImH': (-4.3, 35), 'ReH': (-6.5, 8)}
    ylims = {'ImH': (-4.3, 35), 'ReH': (-6, 8),
             'ImE': (-40, 35), 'ReE': (-15, 30),
             'ImEt': (-50, 100), 'ReEt': (-50, 100),
             'ImHt': (-10, 20), 'ReHt': (-10, 20)}
    # Plot panels
    ts = [-0.3, -0.12]
    for n in range(len(cffs)):
        for nt in range(len(ts)):
            cff = cffs[n]
            # all-x logarithmic
            ax = fig.add_subplot(len(cffs), 2, 2*n+nt+1)
            ax.set_xscale('log')  # x-axis to be logarithmic
            panel(ax, xaxis='xi', xs=logxvals, kins={
                    'yaxis':cff, 't':ts[nt], 'Q2':4.}, **kwargs)
            ax.axhspan(-0.0005, 0.0005, facecolor='g', alpha=0.6)  # horizontal bar
            ax.text(0.31, 0.02, "$t = %s\\, {\\rm GeV}^2$" % str(ts[nt]), 
                    transform=ax.transAxes, fontsize=18)
            ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%s'))
            if nt == 0:
                # plot data-region vertical band for left panels
                #ax.axvspan(0.025, 0.136, facecolor='g', edgecolor='black', alpha=0.2)
                ax.set_ylabel(toTeX['%s' % cff], fontsize=20)
            if nt == 1:
                # no y tick labels on right panels
                ax.set_yticklabels([])
            if n == 0:
                # no x tick labels on upper panels
                ax.set_xticklabels([])
            if n == 1:
                # x-label only on lower panels
                ax.set_xlabel(toTeX['xixB'], fontsize=18)
            if n == 0 and nt == 0:
                # put legend in first panel
                ax.legend(loc='upper right')
                ax.legend(bbox_to_anchor=(0.85, 0.85), loc='upper right',
                        borderaxespad=0.).draw_frame(0)
                pass
                #ax.legend().draw_frame(0)
                # 
                #ax.text(0.33, 0.95, "data region", transform=ax.transAxes, 
                #        fontsize=14, fontweight='bold', va='top')
            apply(ax.set_ylim, ylims[cff])
            ax.set_xlim(0.005, 1.0)
            #ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.02))  # tickmarks
            for label in ax.get_xticklabels() + ax.get_yticklabels():
                label.set_fontsize(14)
    fig.subplots_adjust(bottom=0.5, wspace=0.0, hspace=0.0)
    if path:
        #fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
        fig.set_size_inches((14, 16))
        mkpdf(path)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def CFFt(cffs=['ImH', 'ReH'], path=None, fmt='png', **kwargs):
    """Makes plots of cffs as function of t given by various theories/models
    
    cffs    -- List of CFFs to be plotted. Each produces two panels: ImF, ReF.

    """
    title = ''
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    colors = ['red', 'brown']     # worm human colors :-)
    nncolors = ['blue', 'green']  # cold computer colors
    linestyles = ['solid', 'dashed']
    # Define abscissas
    tmvals = np.linspace(0.02, 0.6, 10) # right panel
    # ordinates 
    #ylims = {'ImH': (-4.3, 35), 'ReH': (-6.5, 8),
    # Plot panels
    for n in range(len(cffs)):
        cff = cffs[n]
        # small x
        ax = fig.add_subplot(len(cffs), 2, 2*n+1)
        panel(ax, xaxis='tm', xs=tmvals, kins={'yaxis':cff, 'xB':0.001, 'Q2':4.}, **kwargs)
        ax.set_xlabel(toTeX['tm'], fontsize=15)
        ax.set_ylabel(toTeX['%s' % cff], fontsize=18)
        ax.axhspan(-0.0005, 0.0005, facecolor='g', alpha=0.6)  # horizontal bar
        if n == 0:
            ax.legend(loc='upper right')
            ax.legend().draw_frame(0)
            ax.text(0.03, 600, "$x_B = 0.001$",# transform=ax.transAxes, 
                    fontsize=15)
            ax.text(0.03, 40, "$Q^2 = 4\\, {\\rm GeV}^2$",# transform=ax.transAxes, 
                    fontsize=15)
        # moderate x
        ax = fig.add_subplot(len(cffs), 2, 2*n+2)
        panel(ax, xaxis='tm', xs=tmvals, kins={'yaxis':cff, 'xB':0.05, 'Q2':4.}, **kwargs)
        ax.set_xlabel(toTeX['tm'], fontsize=15)
        ax.set_ylabel(toTeX['%s' % cff], fontsize=18)
        ax.axhspan(-0.0005, 0.0005, facecolor='g', alpha=0.6)  # horizontal bar
        if n == 0:
            ax.legend(loc='upper right')
            ax.legend().draw_frame(0)
            ax.text(0.03, 4, "$x_B = 0.05$",# transform=ax.transAxes, 
                    fontsize=15)
            ax.text(0.03, 0.4, "$Q^2 = 4\\, {\\rm GeV}^2$",# transform=ax.transAxes, 
                    fontsize=15)
    fig.subplots_adjust(bottom=0.1)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def CFF3(path=None, fmt='png', **kwargs):
    """Makes plots of cffs given by various theories/models
    
    cffs    -- List of CFFs to be plotted. Each produces two panels.

    """
    title = ''
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    colors = ['red', 'brown']     # worm human colors :-)
    nncolors = ['blue', 'green']  # cold computer colors
    linestyles = ['solid', 'dashed']
    # Define abscissas
    logxvals = np.logspace(-5.0, -2, 20) 
    # ordinates 
    leftcffs = ['ImH', 'ImE']
    rightcffs = ['ReH', 'ReE']
    for n in range(2):
        # left panel
        leftcff = leftcffs[n]
        ax = fig.add_subplot(2, 2, 2*n+1)
        ax.set_xscale('log')  # x-axis to be logarithmic
        panel(ax, xaxis='xi', xs=logxvals, kins={'yaxis':leftcff, 't':-0.2, 'Q2':4.}, **kwargs)
        ax.set_xlabel(toTeX['xixB'], fontsize=15)
        ax.set_ylabel(toTeX['%s' % leftcff], fontsize=18)
        ax.axhspan(-0.0005, 0.0005, facecolor='g', alpha=0.6)  # horizontal bar
        # right panel
        rightcff = rightcffs[n]
        ax = fig.add_subplot(2, 2, 2*n+2)
        ax.set_xscale('log')  # x-axis to be logarithmic
        panel(ax, xaxis='xi', xs=logxvals, kins={'yaxis':rightcff, 't':-0.2, 'Q2':4.}, **kwargs)
        ax.set_xlabel(toTeX['xixB'], fontsize=15)
        ax.set_ylabel(toTeX['%s' % rightcff], fontsize=18)
        ax.axhspan(-0.0005, 0.0005, facecolor='g', alpha=0.6)  # horizontal bar
    fig.subplots_adjust(bottom=0.10)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def GPDt(path=None, fmt='png', **kwargs):
    """Makes plots of GPDs as function of t given by various theories/models

    """
    title = ''
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    colors = ['red', 'blue', 'brown', 'violet']
    linestyles = ['solid', 'dashed', 'solid', 'dashed']
    # Define abscissas
    tmvals = np.linspace(0.02, 0.8, 20) # right panel
    # ordinates 
    #ylims = {'ImH': (-4.3, 35), 'ReH': (-6.5, 8),
    # Plot panels
    obses = ['gpdHtrajQ', 'gpdHtrajG', 'gpdEtrajQ', 'gpdEtrajG']
    for n in range(len(obses)):
        obs = obses[n]
        ax = fig.add_subplot(2, 2, n+1)
        panel(ax, xaxis='tm', xs=tmvals, kins={'yaxis':obs, 'xB':0.001998, 'Q2':4.},
                xF=(obs[-1]=='Q'), **kwargs)
        ax.set_xlabel(toTeX['tm'], fontsize=15)
        #ax.set_ylabel(toTeX['%s' % cff], fontsize=18)
        ax.set_ylabel(obs, fontsize=18)
        ax.axhspan(-0.0005, 0.0005, facecolor='g', alpha=0.6)  # horizontal bar
        if n == 0:
            ax.legend(loc='upper right')
            ax.legend().draw_frame(0)
            ax.text(0.3, 60, "$x_B = 0.001$",# transform=ax.transAxes, 
                    fontsize=15)
            ax.text(0.3, 3000, "$Q^2 = 4\\, {\\rm GeV}^2$",# transform=ax.transAxes, 
                    fontsize=15)
    fig.subplots_adjust(bottom=0.1)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig


def jbod(path=None, fmt='png', **kwargs):
    """Plot "just a bunch of data" with lines
    
    """
    title = 'just-a-bunch-of-data'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    ax = fig.add_subplot(1, 1, 1)
    n = 1
    if isinstance(kwargs['points'][0], Data.DataPoint):
        kwargs['points'] = [kwargs['points']]
    for pts in kwargs['points']:
        label = string.join(set('%s-%s-%s' % (pt.collaboration, 
            str(pt.year)[-2:], pt.yaxis) for pt in pts), ' + ')
        ax.text(n, 0, label)
        for pt in pts:
            pt.npt = n
            n += 1
    panel(ax, xaxis='npt', **kwargs)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def bspace(th, parsets=False, path=None, fmt='png', error=False, **kwargs):
    """Makes b-space distribution plot of theory th. and pars ...

    """
    title = ''
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    colors = ['red', 'blue', 'black', 'violet']
    linestyles = ['solid', 'dashed', 'solid', 'dashed']
    pt = Data.DummyPoint()
    pt.xi = 0.001
    pt.Q2 = 4.0
    pt.t = -0.2 # not used
    cov = True
    if not parsets:
        parsets = ['{}']
        cov = False
    # Define abscissas
    bvals = np.linspace(-1.8, 1.8, 50)
    ###    FIRST panel
    k = 0
    ax = fig.add_subplot(1, 3, 1)
    for pars in parsets:
        th.m.parameters.update(eval(pars))
        if (error and cov): th.m.covariance = eval(pars+'cov')
        # ordinates 
        ys = []
        for b in bvals:
            pt.by = b
            pt.bx = 0
            ys.append(th.predict(pt, observable='gpdHQb', error=error))
        ys = np.array(ys)
        if error:
            yup, ydown = np.array([(m+err, m-err) for m,err in ys]).transpose()
            x = plt.concatenate( (bvals, bvals[::-1]) )
            y = pt.xi*plt.concatenate( (yup, ydown[::-1]) )
            ax.fill(x, y, alpha=0.5, color=colors[k+1], **kwargs)
        else:
            ax.plot(bvals, pt.xi*ys, color=colors[k])
        k += 1
    ax.set_xlabel('$b_y \\quad {\\rm [fm]}$', fontsize=24)
    ax.set_ylabel('$x q(x, b) \\quad {\\rm [fm}^{-2}{\\rm ]}$', fontsize=24)
    ax.axhline(y=0, color="black", linestyle="--", linewidth=1)    
    ax.axvline(x=0, color="black", linestyle="--", linewidth=1)    
    ax.set_ylim(0, 1.8)
    ###    SECOND panel
    k = 0
    ax = fig.add_subplot(1, 3, 2)
    for pars in parsets:
        th.m.parameters.update(eval(pars))
        if (error and cov): th.m.covariance = eval(pars+'cov')
        # ordinates 
        ys = []
        for b in bvals:
            pt.by = b
            pt.bx = 0
            ys.append(th.predict(pt, observable='gpdHQbpol', error=error))
        ys = np.array(ys)
        if error:
            yup, ydown = np.array([(m+err, m-err) for m,err in ys]).transpose()
            x = plt.concatenate( (bvals, bvals[::-1]) )
            y = pt.xi*plt.concatenate( (yup, ydown[::-1]) )
            ax.fill(x, y, alpha=0.5, color=colors[k+1], **kwargs)
        else:
            ax.plot(bvals, pt.xi*ys, color=colors[k])
        k += 1
    ax.set_xlabel('$b_y \\quad {\\rm [fm]}$', fontsize=24)
    ax.set_ylabel('$x q^{\\uparrow}(x, b) \\quad {\\rm [fm}^{-2}{\\rm ]}$', fontsize=24)
    ax.axhline(y=0, color="black", linestyle="--", linewidth=1)    
    ax.axvline(x=0, color="black", linestyle="--", linewidth=1)    
    ax.set_ylim(0, 1.8)
    ###    THIRD panel
    k = 0
    ax = fig.add_subplot(1, 3, 3)
    for pars in parsets:
        th.m.parameters.update(eval(pars))
        if (error and cov): th.m.covariance = eval(pars+'cov')
        # ordinates 
        ys = []
        for b in bvals:
            pt.by = b
            pt.bx = 0
            ys.append(th.predict(pt, observable='gpdHGb', error=error))
        ys = np.array(ys)
        if error:
            yup, ydown = np.array([(m+err, m-err) for m,err in ys]).transpose()
            x = plt.concatenate( (bvals, bvals[::-1]) )
            y = plt.concatenate( (yup, ydown[::-1]) )
            ax.fill(x, y, alpha=0.5, color=colors[k+1], **kwargs)
        else:
            ax.plot(bvals, ys, color=colors[k])
        k += 1
    ax.set_xlabel('$b_y \\quad {\\rm [fm]}$', fontsize=24)
    ax.set_ylabel('$x g(x, b) \\quad {\\rm [fm}^{-2}{\\rm ]}$', fontsize=24)
    ax.axhline(y=0, color="black", linestyle="--", linewidth=1)    
    ax.axvline(x=0, color="black", linestyle="--", linewidth=1)    
    ax.set_ylim(0, 10)
    fig.subplots_adjust(bottom=0.5)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def bspace2D(th, flavor='Q', path=None, fmt='png', **kwargs):
    """Makes b-space distribution plot of theory th.

    """
    title = ''
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    pt = Data.DummyPoint()
    pt.xi = 0.001
    pt.Q2 = 4.0
    pt.t = -0.2
    # ordinates 
    bxs = np.linspace(-1.7, 1.7, 20)
    bys = np.linspace(-1.7, 1.7, 20)
    X, Y = np.meshgrid(bxs, bys)
    Z = []
    Zpol = []
    for by in bys:
        pt.by = by
        aux = []
        auxpol = []
        for bx in bxs:
            pt.bx = bx
            if flavor=='Q':
                aux.append(pt.xi*th.m.gpdHQb(pt))
                auxpol.append(pt.xi*th.m.gpdHQbpol(pt))
            elif flavor=='G':
                aux.append(th.m.gpdHGb(pt))
                auxpol.append(th.m.gpdHGbpol(pt))
            else:
                raise ValueError('Flavor %s is unknown' % flavor)
        Z.append(aux)
        Zpol.append(auxpol)
    ax = fig.add_subplot(1, 2, 1)
    CS = ax.contourf(X, Y, Z, cmap=plt.cm.jet)
    plt.colorbar(CS)
    ax.set_xlabel('bx', fontsize=15)
    ax.set_ylabel('by', fontsize=18)
    props = dict(color="black", linestyle="--", linewidth=1)
    ax.axvline(x=0, **props)
    ax.axhline(y=0, **props)    
    ax = fig.add_subplot(1, 2, 2)
    ax.contourf(X, Y, Zpol, cmap=plt.cm.jet)
    ax.set_xlabel('bx', fontsize=15)
    ax.set_ylabel('by', fontsize=18)
    props = dict(color="black", linestyle="--", linewidth=1)
    ax.axvline(x=0, **props)
    ax.axhline(y=0, **props)    
    #fig.subplots_adjust(bottom=0.1)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def markus(th, error=False, path=None, fmt='png', **kwargs):
    """Makes b-space distribution plot of theory th.

    """
    title = ''
    fig = plt.figure(figsize=(9,10), edgecolor='b')
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    pt = Data.DummyPoint()
    pt.xi = 0.001
    pt.Q2 = 4.0
    pt.t = -0.2
    #
    resolution = 20
    bxs = np.linspace(-1.7, 1.7, resolution)
    bys = np.linspace(-1.7, 1.7, resolution)
    ### 1D panel
    ax = fig.add_subplot(2, 1, 1)
    ys = []
    for b in bxs:
        pt.by = 0
        pt.bx = b
        ys.append(th.predict(pt, observable='gpdHQbpol', error=error))
    ys = np.array(ys)
    if error:
        yup, ydown = np.array([(m+err, m-err) for m,err in ys]).transpose()
        x = plt.concatenate( (bxs, bxs[::-1]) )
        y = pt.xi*plt.concatenate( (yup, ydown[::-1]) )
        ax.fill(x, y, alpha=0.7, color='blue', **kwargs)
    else:
        ax.plot(bxs, pt.xi*ys, color='blue')
    ax.set_xlim(-1.7, 1.7)
    ax.set_ylim(0, 2.1)
    ax.set_xlabel('$b_x \\; {\\rm [fm]}$', fontsize=24)
    ax.set_ylabel('$x q^{\\Uparrow}(x, b) \\; {\\rm [fm}^{-2}{\\rm ]}$', fontsize=24)
    props = dict(color="green", linestyle="-.", linewidth=1)
    ax.axvline(x=0, **props)
    for label in ax.get_yticklabels():
        label.set_fontsize(16)
    ### 2D panel
    X, Y = np.meshgrid(bxs, bys)
    Z = []
    Zpol = []
    for by in bys:
        pt.by = by
        auxpol = []
        for bx in bxs:
            pt.bx = bx
            auxpol.append(pt.xi*max(th.m.gpdHQbpol(pt),0.01))
        Zpol.append(auxpol)
    ax = fig.add_subplot(2, 1, 2)
    #actual contour plot
    im = ax.contourf(X, Y, Zpol, np.linspace(0, 2, 40), cmap=plt.cm.Blues)
    # Just the single border line - contour style
    #cs = ax.contour(X, Y, Zpol, np.array([0.15]), linestyle='--')
    #for c in cs.collections:
    #    c.set_linestyle('dashed')
    #    c.set_linewidth(0.5)
    # Just the single border line - plot style
    phis = np.linspace(0, 2*np.pi, 20)
    bees = []
    for phi in phis:
        bees.append(th.m.rth(pt, phi))
    bees = np.array(bees)
    ax.plot(bees*np.cos(phis), bees*np.sin(phis), 'r--')
    # colorbar
    cb = plt.colorbar(im)
    cb.ax.set_yticks([0, 0.5, 1])
    cb.ax.set_yticklabels(['0.', '1.', '2.'])
    l,b,w,h = fig.gca().get_position().bounds
    ll,bb,ww,hh = cb.ax.get_position().bounds
    cb.ax.set_position([ll-0.07, b+0.0*h, ww, h*1.1])
    ax.set_ylim(-1.7, 1.7)
    ax.set_xlim(-1.7, 1.7)
    ax.set_xlabel('$b_x \\; {\\rm [fm]}$', fontsize=24)
    ax.set_ylabel('$b_y \\; {\\rm [fm]}$', fontsize=24)
    ax.axvline(x=0, **props)
    ax.axhline(y=0, **props)    
    #ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.02))  # tickmarks
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(16)
    fig.subplots_adjust(left=0.2, right=0.7, hspace=0)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def beff(m, pars1, pars2, path=None, fmt='png', **kwargs):
    """Makes beff slope plot

    """
    title = ''
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    colors = ['red', 'blue', 'brown', 'violet']
    linestyles = ['solid', 'dashed', 'solid', 'dashed']
    pt = Data.DummyPoint()
    pt.t = -0.2 # not used
    ax = fig.add_subplot(1, 2, 1)
    ax.set_xscale('log')  # x-axis to be logarithmic
    # Define abscissas
    logxBvals = np.logspace(-4, -2)
    # ordinates 
    n = 0
    m.parameters.update(pars1)
    for Q2 in [4., 8., 16.]:
        pt.Q2 = Q2
        beffs = []
        for xB in logxBvals:
            pt.xB = xB
            pt.xi = xB/(2.-xB)
            beffs.append(m.beff(pt))
        ax.plot(logxBvals, beffs, linestyle=linestyles[n], color=colors[n])
        n += 1
    ax.set_xlabel('xB', fontsize=15)
    ax.set_ylabel('b_eff(ImH)', fontsize=18)
    ax.axhspan(-0.0005, 0.0005, facecolor='g', alpha=0.6)  # horizontal bar
    ax.set_xlim(1e-4, 1e-2)
    ax.set_ylim(0, 7)
    ax = fig.add_subplot(1, 2, 2)
    ax.set_xscale('log')  # x-axis to be logarithmic
    # Define abscissas
    logxBvals = np.logspace(-4, -2)
    # ordinates 
    n = 0
    m.parameters.update(pars2)
    for Q2 in [4., 8., 16.]:
        pt.Q2 = Q2
        beffs = []
        for xB in logxBvals:
            pt.xB = xB
            pt.xi = xB/(2.-xB)
            beffs.append(m.beff(pt))
        ax.plot(logxBvals, beffs, linestyle=linestyles[n], color=colors[n])
        n += 1
    ax.set_xlabel('xB', fontsize=15)
    ax.set_ylabel('b_eff(ImH)', fontsize=18)
    ax.axhspan(-0.0005, 0.0005, facecolor='g', alpha=0.6)  # horizontal bar
    ax.set_xlim(1e-4, 1e-2)
    ax.set_ylim(0, 7)
    fig.subplots_adjust(bottom=0.1)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig

def binplot(path=None, fmt='png', **kwargs):
    """Plot local fit of a single kinematic bin."""
    title = kwargs.pop('title', 'N/A')
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    ax = fig.add_subplot(1, 1, 1)
    n = 1
    pts = kwargs['points']
    for pt in pts:
        pt.npt = n
        n += 1
    panel(ax, xaxis='npt', **kwargs)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    ax.set_xticks(range(1,n))
    ax.set_xticklabels([OBStoTeX[(pt.y1name, int(pt.FTn))] for pt in pts])
    for pt, label in zip(pts, ax.get_xticklabels()):
        if (pt.y1name == 'BTSA' and int(pt.FTn) == 1):
            label.set_fontsize(12)
        elif (pt.y1name == 'AUTI' and int(pt.FTn) == 0):
            label.set_fontsize(12)
        elif (pt.y1name == 'AUTI' and int(pt.FTn) == -1):
            label.set_fontsize(12)
        else:
            label.set_fontsize(16)
    for label in ax.get_xticklabels()[6:]:
            label.set_rotation(60)
    #xtickNames = plt.setp(ax)
    #plt.setp(xtickNames, rotation=45)
    fig.subplots_adjust(bottom=0.2)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig
