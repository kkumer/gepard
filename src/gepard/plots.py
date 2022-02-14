"""Plotting scripts. Not everything works."""

import itertools
import os

import matplotlib  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import numpy as np
import pandas as pd  # type: ignore
from scipy.interpolate import InterpolatedUnivariateSpline  # type: ignore

from . import constants, data

# if os.sys.platform == 'win32':
#     matplotlib.use('WxAgg')
# #else: #linux
# #    matplotlib.use('TkAgg')
# #import pylab as plt
# # Following works better in jupyter:
# matplotlib.rcParams['figure.figsize'] = (7.0, 5.0)

#  subplots_adjust options and their meanings:
# left  = 0.125  # the left side of the subplots of the figure
# right = 0.9    # the right side of the subplots of the figure
# bottom = 0.0   # the bottom of the subplots of the figure
# top = 0.9      # the top of the subplots of the figure
# wspace = 0.2   # the amount of width reserved for blank space between subplots
# hspace = 0.2   # the amount of height reserved for white space between subplots
# Usage: fig.subplots_adjust(wspace=0.7)

DISpoints = data.dset[201].copy()
for id in range(202, 213):
    DISpoints += data.dset[id]

#################################################################
#                                                               #
#   [1]  Universal stuff, used to build up panels of figures    #
#                                                               #
#################################################################

def _axpoints(ax, pts, xaxis, **kwargs):
    """Make an errorbar plot defined by a single set of data points.

    Args:
        ax: matplotlib axis object
        pts: single set of data points
        xaxis: 'Q2', 't', 'xB', ...
        **kwargs: keyword arguments

    """
    kwargs.pop('justbars', False)  # Kludge, FIXME
    pts = [pt.copy() for pt in pts]
    for pt in pts:
        pt.from_conventions()
    xvals = [getattr(pt, xaxis) for pt in pts]
    yvals = [pt.val for pt in pts]
    yerrs = [pt.err for pt in pts]
    ax.errorbar(xvals, yvals, yerrs, **kwargs)


def _axline(ax, fun, points, xaxis, **kwargs):
    """Make a line(s) corresponding to points.

    Args:
        ax: matplotlib axis object
        fun: real-valued function of DataPoint instance
        points: list of sets of DataPoint instances
        xaxis: 'Q2', 't', 'xB', ...
        **kwargs: keyword arguments

    """
    for pts in points:
        if 'phi' in pts[0] and pts[0].frame == 'Trento':
            # go back to degrees
            oldpts = [pt.copy() for pt in pts]
            for pt in oldpts:
                pt.from_conventions()
            xvals = [getattr(pt, xaxis) for pt in oldpts]
        else:
            xvals = [getattr(pt, xaxis) for pt in pts]
        yvals = [fun(pt) for pt in pts]
        kwargs.pop('justbars', False)  # Kludge
        if kwargs.pop('xF', False):  # do we want  xF(x) plotted?
            xis = np.array([pt.xi for pt in pts])
            yvals = xis*np.array(yvals)
        ax.plot(xvals, yvals, **kwargs)


def _axband(ax, fun, pts, xaxis, **kwargs):
    """Make a band corresponding to points.

    Args:
        ax: matplotlib axis object
        fun: function of DataPoint instance returning tuple
             (mean, err) or (mean, err+, err-) defining errorband
        pts: set of DataPoint instances (just one set!)
        xaxis: 'Q2', 't', 'xB', ...
        **kwargs: keyword arguments

    """
    if 'phi' in pts[0] and pts[0].frame == 'Trento':
        # go back to degrees
        oldpts = [pt.copy() for pt in pts]
        for pt in oldpts:
            pt.from_conventions()
        xvals = [getattr(pt, xaxis) for pt in oldpts]
    else:
        xvals = [getattr(pt, xaxis) for pt in pts]
    res = [fun(pt) for pt in pts]
    if len(res[0]) == 2:
        # we have symmetric error
        res = [(m, err, err) for m, err in res]
    up, down = np.array([(m+errp, m-errm) for m, errp, errm in res]).transpose()

    if kwargs.pop('justbars', False):  # do we want just errorbars for theory?
        res = list(zip(up, down))
        # FIXME: symmetric error assumed
        yvals, yerrs = np.array([((up+down)/2, (up-down)/2) for
                                 up, down in res]).transpose()
        xvals = [x+0.1 for x in xvals]
        ax.errorbar(xvals, yvals, yerrs, marker='d', color='red', linestyle='None',
                    elinewidth=3)
    else:
        x = np.concatenate((xvals, xvals[::-1]))
        y = np.concatenate((up, down[::-1]))
        if kwargs.pop('xF', False):  # do we want  xF(x) plotted?
            xis = np.array([pt.xi for pt in pts+pts[::-1]])
            y = xis*y
        ax.fill(x, y, alpha=0.5, **kwargs)


def panel(ax, points=None, lines=None, bands=None, xaxis=None, xs=None,
          kins={}, kinlabels=None, **kwargs):
    """Plot datapoints together with fit/theory line(s) and band(s).

    Args:
        ax: matplotlib axis object
        points: (list of) `DataSet` instance(s) to be plotted.
        lines: (list of) 'theorie(s)' describing curves for plotting.
        bands: (list of) neural network 'theorie(s)' or 'theorie(s)' with
               defined uncertianties, defining band by their mean and std.dev.
        xaxis: abscissa variable; if None, last of sets.xaxes is taken;
               if 'points' all points are just put equidistantly along x-axis
        xs: list of xvalues, needed if not specified by points
        kins: additional kinematic data
        kinlabels: list of constant kinematic variables whose values will
                   be extracted from points and put on plot as annotation
        **kwargs: keyword arguments passed to matplotlib

    """
    if points:
        # If not provided, take xaxis from last axis from last set
        if not xaxis:
            xaxis = points[-1].xaxes[-1]
        if isinstance(points[0], data.DataPoint):
            points = [points]

        pointshapes = ['o', 's', '^', 'd', 'o', 's', '^', 'd']
        pointcolors = ['blue', 'red', 'green', 'red', 'brown', 'blue', 'black']
        pointcolors = ['black', 'red', 'green', 'red', 'brown', 'blue', 'black']
        setn = 0
        for pts in points:
            _axpoints(ax, pts, xaxis, linestyle='None', elinewidth=1,
                      label='{} {}'.format(pts[0].collaboration, pts[0].year),
                      marker=pointshapes[setn], color=pointcolors[setn], **kwargs)
            setn += 1
    else:
        # We have to create set of DataPoints.
        # (It is up to caller to provide everything needed
        # via xaxis and kins)
        pts = []
        for x in xs:
            pt = data.DataPoint(init=kins)
            setattr(pt, xaxis, x)
            # workaround for imperfect fill_kinematics
            if hasattr(pt, 'xi'):
                pt.xB = 2*pt.xi/(1.+pt.xi)
            # end of workaorund
            data._fill_kinematics(pt)
            pts.append(pt)
        points = [pts]

    if lines:
        if not isinstance(lines, list):
            lines = [lines]
        linecolors = ['red', 'black', 'blue', 'green', 'darkorchid', 'olive',
                      'darkcyan', 'indianred', 'red', 'black', 'blue', 'green',
                      'darkorchid', 'olive', 'darkcyan', 'indianred']
        linestyles = ['-', '--', '-.', ':', '-', '--', '-.', ':',
                      '-', '--', '-.', ':', '-', '--', '-.', ':']
        linen = 0
        for line in lines:
            _axline(ax, lambda pt: line.predict(pt, orig_conventions=True), points,
                    xaxis=xaxis, color=linecolors[linen], linestyle=linestyles[linen],
                    linewidth=2, alpha=1.0, label=r'\texttt{{{}}}'.format(line.name),
                    **kwargs)
            linen += 1

    if bands:
        if not isinstance(bands, list):
            bands = [bands]
        bandcolors = ['xkcd:teal', 'red', 'darkcyan', 'blue', 'green', 'purple']
        hatches = ['\\\\', '////', '|', '.']
        bandn = 0
        for band in bands:
            for pts in points:
                _axband(ax, lambda pt: band.predict(pt, error=True,
                        orig_conventions=True), pts, xaxis=xaxis,
                        hatch=hatches[bandn],
                        # color=bandcolors[bandn],
                        # facecolor=bandcolors[bandn], # band  color
                        facecolor='none',
                        edgecolor=bandcolors[bandn],  # band edge color
                        linewidth=2,
                        label=r'\texttt{{{}}}'.format(band.name), **kwargs)
            bandn += 1

    if kinlabels:
        # constant kinematic variables positioning
        xvals = []
        yvals = []
        for pts in points:
            for pt in pts:
                if xaxis == 'phi' and pt.frame == 'Trento':
                    xvals.append((np.pi - pt.phi)/np.pi * 180.)
                else:
                    xvals.append(getattr(pt, xaxis))
                yvals.append(pt.val)
        labx = min(0, min(xvals)) + (max(xvals) - min(0, min(xvals))) * 0.15
        laby = min(0, min(yvals)) + (max(yvals) - min(0, min(yvals))) * 0.04
        labtxt = ""
        for lab in kinlabels:
            try:
                labtxt += constants.toTeX[lab] + ' = ' + str(getattr(points[0], lab)) + ', '
            except AttributeError:
                # If dataset doesn't have it, all points should have it
                labtxt += constants.toTeX[lab] + ' = ' + str(getattr(points[0][0],lab)) + ', '
        ax.text(labx, laby, labtxt[:-2], fontsize=14)


def thline(th, pts, ex_pt, npts=16):
    """Interpolated theory lines.

    For kinematics in pandas dataframe pts and DataPoint ex_pt.

    """
    tms = np.linspace(pts.tm.min()-0.02, pts.tm.max()+0.03, npts)
    try:
        spl_xB = InterpolatedUnivariateSpline(pts.tm.values, pts.xB.values, k=3)
        spl_Q2 = InterpolatedUnivariateSpline(pts.tm.values, pts.Q2.values, k=3)
    except:  # too few points in a bin
        spl_xB = InterpolatedUnivariateSpline(pts.tm.values, pts.xB.values, k=1)
        spl_Q2 = InterpolatedUnivariateSpline(pts.tm.values, pts.Q2.values, k=1)
    ys = []
    pt = ex_pt.copy()
    for tm in tms:
        pt.tm = tm
        pt.xB = spl_xB(tm)[()]
        pt.Q2 = spl_Q2(tm)[()]
        del pt.t
        del pt.W
        data._fill_kinematics(pt)
        ys.append(th.predict(pt))
    return tms, np.array(ys)

#################################################################
#                                                              ##
#   [2]  Comparisons with various measurements                 ##
#                                                              ##
#################################################################


def HERMES09AC(path=None, fmt='png', **kwargs):
    """Plot HERMES 09 AC."""
    id = 32
    title = 'HERMES 09 AC'
    fig = plt.figure()
    fig.suptitle(title)
    xaxes = ['tm', 'xB', 'Q2']
    # we have 3x6 points
    for x in range(3):
        ax = fig.add_subplot(3, 1, x+1)
        ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))  # tickmarks
        panel(ax, points=data.dset[id][x*6+18:x*6+6+18], xaxis=xaxes[x], **kwargs)
        #apply(ax.set_ylim, ylims[(-0.30, 0.05)])
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig

def HERMES12(path=None, fmt='png', **kwargs):
    """Plot HERMES combined AC and ALU data with fit lines."""

    title = ''
    fig = plt.figure(figsize=(14, 14))
    fig.suptitle(title)
    xaxes = ['tm', 'xB', 'Q2']
    ylims = [(-0.05, 0.3), (-0.15, 0.15), (-0.45, 0.05)]
    # we have 3x18=54 points to be separated in nine panels six points each:
    for y, id, shift in zip(list(range(3)), [67, 67, 68], [18, 0, 0]):
        for x in range(3):
            npanel = 3*y + x + 1  # 1, 2, ..., 9
            ax = fig.add_subplot(3, 3, npanel)
            ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))  # tickmarks
            panel(ax, points=data.dset[id][x*6+shift:x*6+6+shift], xaxis=xaxes[x], **kwargs)
            ax.set_ylim(*ylims[y])
            if (npanel % 3) != 1:
                # Leave labels only on leftmost panels
                ax.set_ylabel('')
            else:
                ylabels = ['$A_{C}^{\\cos \\phi}$', '$A_{C}^{\\cos 0\\phi}$', '$A_{LU,I}^{\\sin\\phi}$']
                ax.set_ylabel(ylabels[int((npanel-1)/3)], fontsize=16)
            if npanel < 7:
                # Leave labels only on lowest panels
                ax.set_xlabel('')
            else:
                xlabels = ['$-t\\; [{\\rm GeV}^2]$', '$x_B$', '$Q^2\\; [{\\rm GeV}^2]$']
                ax.set_xlabel(xlabels[npanel-7], fontsize=16)

            if (npanel % 3) == 2:
                # Adjust x-axis on middle column
                ax.set_xlim(0.04, 0.25)
                ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))  # tickmarks
            if npanel == 2:
                handles, labels = ax.get_legend_handles_labels()
                handles.pop()
                labels.pop()
                ax.legend(handles, labels, handlelength=2.0,
                 prop=matplotlib.font_manager.FontProperties(size="small"),
                    loc="upper center").draw_frame(0)
    fig.subplots_adjust(wspace=0.1, hspace=0.1)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig


def HERMES10t(lines=None, path=None, fmt='png'):
    """Plot AUL(t) and ALL(t) for HERMES 2010 data."""
    title = 'HERMES-10'
    NPTS = 30
    dfAUL = data.select(data.dset[52], criteria=['FTn == -1'])[:4].df()
    dfALL0 = data.select(data.dset[53], criteria=['FTn == 0'])[:4].df()
    dfALL1 = data.select(data.dset[53], criteria=['FTn == 1'])[:4].df()
    dfs = [dfAUL, dfALL1, dfALL0]
    ylabels = [r'$A_{UL,+}^{\sin\phi}$',
               r'$-A_{LL,+}^{\cos\phi}$',
               r'$A_{LL,+}^{\cos 0\phi}$']
    ylims = [(-0.32, 0.08), (-0.32, 0.22), (-0.28, 0.55)]
    yticks = [0.1, 0.1, 0.2]
    styles = ['b-', 'r--', 'g-.', 'k:']
    fig, axs = plt.subplots(3, 1, sharex=True, figsize=[7,9])
    for p, ax in enumerate(axs):
        pts = dfs[p]
        ax.errorbar(pts.tm.values, pts.val.values, pts.err.values,
                    linestyle='None', marker='s', color='black')
        if not isinstance(lines, list):
            lines = [lines]
        for l, th in enumerate(lines):
            xs, ys = thline(th, pts, pts.pt[0], npts=NPTS)
            ax.plot(xs, ys, styles[l], label=th.name)
        ax.set_xlim(0.0, 0.51)
        ax.set_ylim(*ylims[p])
        ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
        ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(yticks[p]))
        ax.set_ylabel(ylabels[p], fontsize=18)
        ax.set_xlabel('$-t\\; [{\\rm GeV}^2]$', fontsize=16)
    axs[1].legend(loc='lower center').draw_frame(0)
    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig


def HERMES12t(lines=None, path=None, fmt='png'):
    """Plot ALU(t) and AC(t) for HERMES data + recoil."""
    title = 'HERMES-ALU-AC-recoil'
    NPTS = 20
    dfAC = data.dset[67][18:18+6].df()
    dfALU = data.dset[68][:6].df()
    dfALUrec = data.dset[122][:4].df()
    dfs = [dfALU, dfAC]
    ylabels = [r'$A_{LU,I}^{\sin\phi}$', r'$-A_{C}^{\cos\phi}$']
    styles = ['b-', 'r--', 'g-.', 'k:']
    fig, axs = plt.subplots(2, 1, sharex=True, figsize=[7, 6])
    for p, ax in enumerate(axs):
        pts = dfs[p]
        ax.errorbar(pts.tm.values, pts.val.values, pts.err.values,
                    linestyle='None', marker='s', color='black')
        if p == 0:
            # add recoil data
            ax.errorbar(dfALUrec.tm.values, dfALUrec.val.values,
                        dfALUrec.err.values, fillstyle='none',
                        linestyle='None', marker='s', color='black')
        if not isinstance(lines, list):
            lines = [lines]
        for l, th in enumerate(lines):
            xs, ys = thline(th, pts, pts.pt[0], npts=NPTS)
            ax.plot(xs, ys, styles[l], label=th.name)
        ax.set_xlim(0.0, 0.51)
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
        ax.set_ylabel(ylabels[p], fontsize=18)
        ax.set_xlabel('$-t\\; [{\\rm GeV}^2]$', fontsize=16)
    axs[1].legend(loc='upper right').draw_frame(0)
    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig


def HERMES12AC(path=None, fmt='png', **kwargs):
    """Plot HERMES combined AC data with fit lines."""
    title = 'HERMES-12-AC'
    fig = plt.figure()
    fig.suptitle(title)
    xaxes = ['tm', 'xB', 'Q2']
    ylims = [(-0.15, 0.15), (-0.05, 0.3), (-0.15, 0.15)]
    # we have 3x18=54 points to be separated in nine panels six points each:
    for y, id, shift in zip(list(range(3)), [67, 67, 67], [0, 18, 36]):
        for x in range(3):
            npanel = 3*y + x + 1  # 1, 2, ..., 9
            ax = fig.add_subplot(3, 3, npanel)
            ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))  # tickmarks
            panel(ax, points=data.dset[id][x*6+shift:x*6+6+shift], xaxis=xaxes[x], **kwargs)
            ax.set_ylim(*ylims[y])
            if (npanel % 3) != 1:
                # Leave labels only on leftmost panels
                ax.set_ylabel('')
            else:
                ylabels = ['$AC\\; \\cos 0\\phi$', '$AC\\; \\cos \\phi$', '$AC\\; \\cos 2\\phi$']
                ax.set_ylabel(ylabels[int((npanel-1)/3)], fontsize=18)
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
        #fig.show()
    return fig


def HERMES12ALU(path=None, fmt='png', **kwargs):
    """Plot HERMES combined ALU data with fit lines."""

    title = 'HERMES-12-ALU'
    fig = plt.figure()
    fig.suptitle(title)
    xaxes = ['tm', 'xB', 'Q2']
    ylims = [(-0.45, 0.1), (-0.25, 0.25), (-0.15, 0.15)]
    # we have 3x18=54 points to be separated in nine panels six points each:
    for y, id, shift in zip(list(range(3)), [68, 69, 68], [0, 0, 18]):
        for x in range(3):
            npanel = 3*y + x + 1  # 1, 2, ..., 9
            ax = fig.add_subplot(3, 3, npanel)
            ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))  # tickmarks
            panel(ax, points=data.dset[id][x*6+shift:x*6+6+shift], xaxis=xaxes[x], **kwargs)
            ax.set_ylim(*ylims[y])
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
        #fig.show()
    return fig

def HERMES09(path=None, fmt='png', **kwargs):
    """Plot HERMES 0909.3587 AC and ALU data with fit lines."""

    #ids = [2, 4, 5]
    title = ''#HERMES-09'
    fig = plt.figure()
    fig.suptitle(title)
    xaxes = ['tm', 'xB', 'Q2']
    ylims = [(-0.05, 0.3), (-0.45, 0.05)]
    xlims = [(0.0, 0.48), (0.04, 0.27), (0.85, 6)]
    # we have 3x18=54 points to be separated in nine panels six points each:
    for y, id, shift in zip(list(range(2)), [32, 5], [18, 0]):
        for x in range(3):
            npanel = 3*y + x + 1  # 1, 2, ..., 6
            ax = fig.add_subplot(2,3,npanel)
            ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))  # tickmarks
            panel(ax, points=data.dset[id][x*6+shift:x*6+6+shift], xaxis=xaxes[x], **kwargs)
            ax.set_ylim(*ylims[y])
            ax.set_xlim(*xlims[x])
            if (npanel % 3) != 1:
                # Leave labels only on leftmost panels
                ax.set_ylabel('')
                ax.set_yticklabels([])
            else:
                ylabels = ['$A_{C}^{\\cos\\phi}$', '$A_{LU}^{\\sin\\phi}$']
                ax.set_ylabel(ylabels[(npanel-1)/3], fontsize=16)
            if npanel < 4:
                # Leave labels only on lowest panels
                ax.set_xlabel('')
                ax.set_xticklabels([])
            else:
                xlabels = ['$-t\\; [{\\rm GeV}^2]$', '$x_B$', '$Q^2\\; [{\\rm GeV}^2]$']
                ax.set_xlabel(xlabels[npanel-7], fontsize=16)
            if npanel == 3:
                ax.legend(bbox_to_anchor=(0., 1.), loc='upper right',
                        borderaxespad=0.).draw_frame(0)
                pass
        #ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.2))  # tickmarks
            for label in ax.get_xticklabels() + ax.get_yticklabels():
                label.set_fontsize(10)
    fig.subplots_adjust(bottom=0.45, wspace=0.0, hspace=0.0)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig

def HERMES10LP(observable='TSA', path=None, fmt='png', **kwargs):
    """Plot HERMES 1004.0177 TSA data with fit lines."""

    title = 'HERMES-10-'+observable
    fig = plt.figure()
    fig.suptitle(title)
    xaxes = ['tm', 'xB', 'Q2']
    ylims = [(-0.05, 0.3), (-0.15, 0.15), (-0.45, 0.05)]
    xticks = [0.2, 0.1, 2]
    if observable == 'TSA':
        id = 52
        harmonics = [-1, -2, -3]
        fun = 'sin'
        lbl = 'A_{UL}'
    elif observable == 'BTSA':
        id = 53
        harmonics = [0, 1, 2]
        fun = 'cos'
        lbl = 'A_{LL}'
    else:
        raise ValueError('Observable %s nonexistent.' % observable)
    subsets = {}
    for k in range(3):
        subsets[k] = data.select(data.dset[id], criteria=['FTn == %i' % harmonics[k]])
    # we have 3x12=36 points to be separated in nine panels four points each:
    for y in range(2):  # put this to 3 to plot also 3rd harmonic
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
                ylabels = ['$'+lbl+'^{\\'+fun+'%i\\phi}$' % harmonics[k] for k in range(3)]
                ax.set_ylabel(ylabels[(npanel-1)/3], fontsize=16)
            if npanel < 7:
                # Leave labels only on lowest panels
                ax.set_xlabel('')
            else:
                xlabels = ['$-t\\; [{\\rm GeV}^2]$', '$x_B$', '$Q^2\\; [{\\rm GeV}^2]$']
                ax.set_xlabel(xlabels[npanel-7], fontsize=16)
            if (npanel % 3) == 2:
                # Adjust x-axis on middle column
                ax.set_xlim(0.04, 0.25)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig

def HERMES08TP(path=None, fmt='png', **kwargs):
    """Plot HERMES 08 TTSA data with fit lines."""

    #title = 'HERMES-08 TTSA'
    title = ''
    fig = plt.figure()
    fig.suptitle(title)
    xaxes = ['tm', 'xB', 'Q2']
    ylims = [(-0.42, 0.15), (-0.42, 0.15)]
    xlims = [(0.0, 0.48), (0.04, 0.27), (0.85, 6)]
    # we have 2x12=24 points to be separated in six panels four points each:
    for y, id, shift in zip(list(range(2)), [66, 65], [12, 0]):
        for x in range(3):
            npanel = 3*y + x + 1  # 1, 2, ..., 6
            ax = fig.add_subplot(2,3,npanel)
            ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))  # tickmarks
            panel(ax, points=data.dset[id][x*4+shift:x*4+4+shift], xaxis=xaxes[x], **kwargs)
            ax.set_ylim(*ylims[y])
            ax.set_xlim(*xlims[x])
            if (npanel % 3) != 1:
                # Leave labels only on leftmost panels
                ax.set_ylabel('')
                ax.set_yticklabels([])
            else:
                ylabels = ['$A_{UT,I}^{\\sin(\\phi-\\phi_S)\\cos\\phi}$',
                           '$A_{UT,DVCS}^{\\sin(\\phi-\\phi_S)}$']
                ax.set_ylabel(ylabels[(npanel-1)/3], fontsize=16)
            if npanel < 4:
                # Leave labels only on lowest panels
                ax.set_xlabel('')
                ax.set_xticklabels([])
            else:
                xlabels = ['$-t\\; [{\\rm GeV}^2]$', '$x_B$', '$Q^2\\; [{\\rm GeV}^2]$']
                ax.set_xlabel(xlabels[npanel-7], fontsize=16)
            if npanel == 3 and kwargs:
                #ax.legend(bbox_to_anchor=(0., 1.), loc='upper right',
                #        borderaxespad=0.).draw_frame(0)
                pass
            ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.2))  # tickmarks
            ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
            for label in ax.get_xticklabels() + ax.get_yticklabels():
                label.set_fontsize(10)
    fig.subplots_adjust(bottom=0.135, wspace=0.0, hspace=0.0, left=0.15)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig

def CLAS(path=None, fmt='png', **kwargs):
    """Makes plot of CLAS ALU data with fit lines and bands"""

    dataset = data.select(data.dset[25], criteria=['Q2 > 1.57'])
    #dataset = data[25]
    title = ''
    title = 'CLAS (Girod:2007aa)'
    fig, axs = plt.subplots(2, 3, sharey=True, sharex=True)
    axs = axs.reshape(6)
    fig.suptitle(title)
    # Each different Q2 has its own panel
    panelorder = [4, 5, 6, 1, 2, 3]
    npanel = 0
    for Q2, grp in itertools.groupby(dataset, lambda pt: pt.Q2):
        np = panelorder[npanel]
        ax = axs[np-1]
        panelpoints = list(grp)
        # now transform list into DataSet instance ...
        panelset = data.DataSet(panelpoints)
        panelset.__dict__ = data.dset[25].__dict__.copy()
        # ... and plot
        panel(ax, points=panelset, xaxis='tm', kinlabels=['Q2', 'xB'], **kwargs)
        plt.xlim(0.0, 0.6)
        plt.ylim(0.0, 0.4)
        ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
        if (np % 3) != 1:
            # Leave labels only on leftmost panels
            ax.set_ylabel('')
        else:
            ax.set_ylabel('$A_{LU}^{\\sin\\phi}$', fontsize=16)
        if np < 4:
            ax.set_xlabel('')
        else:
            ax.set_xlabel('$-t\\; [{\\rm GeV}^2]$', fontsize=16)
        npanel += 1
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig

def CLAS14(observable='ALU', path=None, fmt='png', **kwargs):
    """Makes plot of CLAS-14 ALU data with fit lines and bands"""

    title = ''
    title = 'CLAS (prelim. 2014)'
    if observable == 'ALU':
        dataset = data.dset[85]
        lbl = '$A_{LU}^{\\sin\\phi}$'
        ymin, ymax = 0, 0.4
    elif observable == 'TSA1':
        dataset = data.DataSet(data.select(data.dset[86], criteria=['FTn == -1']))
        lbl = '$A_{UL}^{\\sin\\phi}$'
        ymin, ymax = 0, 0.4
    elif observable == 'TSA2':
        dataset = data.DataSet(data.select(data.dset[86], criteria=['FTn == -2']))
        lbl = '$A_{UL}^{\\sin2\\phi}$'
        ymin, ymax = -0.1, 0.32
    elif observable == 'BTSA0':
        dataset = data.DataSet(data.select(data.dset[87], criteria=['FTn == 0']))
        lbl = '$A_{LL}^{\\cos0\\phi}$'
        ymin, ymax = 0, 0.8
    elif observable == 'BTSA1':
        dataset = data.DataSet(data.select(data.dset[87], criteria=['FTn == 1']))
        lbl = '$A_{LL}^{\\cos\\phi}$'
        ymin, ymax = -0.3, 0.3
    else:
        raise ValueError('Observable %s unavailable.' % observable)
    fig, axs = plt.subplots(2, 2, sharey=True, sharex=True, figsize=[5,5])
    axs = axs.reshape(4)
    fig.suptitle(title + ' -- ' +  observable)
    # Each different Q2 has its own panel
    panelorder = [3, 4, 1, 2]
    npanel = 0
    for Q2, grp in itertools.groupby(dataset, lambda pt: pt.Q2):
        np = panelorder[npanel]
        ax = axs[np-1]
        panelpoints = list(grp)
        # now transform list into DataSet instance ...
        panelset = data.DataSet(panelpoints)
        panelset.__dict__ = dataset.__dict__.copy()
        # ... and plot
        if np == 2:
            # fake doubling of points to get line visibility
            ptd = panelset[0].copy()
            ptd.tm = ptd.tm+0.02
            panelset.append(ptd)
            panel(ax, points=panelset, xaxis='tm', kinlabels=['Q2', 'xB'], **kwargs)
        else:
            panel(ax, points=panelset, xaxis='tm', kinlabels=['Q2', 'xB'], **kwargs)
        plt.xlim(0.0, 0.6)
        plt.ylim(ymin, ymax)
        ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
        if observable == 'TSA2' or observable == 'BTSA1':
            ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
        if np == 2:
                ax.legend(loc='upper left', borderaxespad=0.).draw_frame(0)
        if (np % 2) != 1:
            # Leave labels only on leftmost panels
            ax.set_ylabel('')
        else:
            ax.set_ylabel(lbl, fontsize=16)
        if np < 3:
            ax.set_xlabel('')
        else:
            ax.set_xlabel('$-t\\; [{\\rm GeV}^2]$', fontsize=16)
        npanel += 1
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig

def CLAS15(observable='ALU', path=None, fmt='png', **kwargs):
    """Makes plot of CLAS 2015 data with fit lines and bands"""

    title = 'CLAS 2015 (Pisano:2015iqa)'
    if observable[:3] == 'ALU':
        dataset = data.dset[94]
        lbl = '$A_{LU}^{\\sin\\phi}$'
        ymin, ymax = 0, 0.49
        old = data.dset[25]
        oldpanels = [old[:1], old[3:6], old[6:9], old[11:14], None]
    elif observable == 'TSA':
        dataset = data.dset[95]
        lbl = '$A_{UL}^{\\sin\\phi}$'
        ymin, ymax = 0, 0.49
    elif observable == 'BTSA0':
        dataset = data.DataSet(data.select(data.dset[96], criteria=['FTn == 0']))
        lbl = '$A_{LL}^{\\cos0\\phi}$'
        ymin, ymax = 0, 0.85
    elif observable == 'BTSA1':
        dataset = data.DataSet(data.select(data.dset[96], criteria=['FTn == 1']))
        lbl = '$A_{LL}^{\\cos\\phi}$'
        ymin, ymax = -0.35, 0.35
    else:
        raise ValueError('Observable %s unavailable.' % observable)
    fig, axs = plt.subplots(3, 2, sharey=True, sharex=True)
    axs = axs.reshape(6)
    fig.suptitle(title + ' -- ' +  observable, fontsize=20)
    #
    panels = [ dataset[:1],
               dataset[1:4],
               dataset[4:7],
               dataset[7:9],
               dataset[9:]]
    for npanel, panelset in enumerate(panels):
        ax = axs[npanel]
        if (npanel == 0) or (npanel == 4):
            # fake doubling of points to get line visibility
            ptd = panelset[0].copy()
            ptd.tm = ptd.tm+0.02
            panelset.append(ptd)
            panel(ax, points=panelset, xaxis='tm', kinlabels=['Q2', 'xB'], **kwargs)
        else:
            if observable == 'ALUOLD':
                panel(ax, points=[panelset, oldpanels[npanel]], xaxis='tm', kinlabels=['Q2', 'xB'], **kwargs)
            else:
                panel(ax, points=panelset, xaxis='tm', kinlabels=['Q2', 'xB'], **kwargs)

        #panel(ax, points=panelset, xaxis='tm', kinlabels=['Q2', 'xB'], **kwargs)
        plt.xlim(0.0, 0.65)
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
        ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.02))
        plt.ylim(ymin, ymax)
        ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
        ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.02))
        if observable == 'BTSA1':
            ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
        if npanel == 3:
                # Create legend here, but draw it on last panel
                handles, labels = ax.get_legend_handles_labels()
        if (npanel % 2) == 1:
            # Remove y-axis label from right panels
            ax.set_ylabel('')
        else:
            ax.set_ylabel(lbl, fontsize=24)
        if npanel < 3:
            ax.set_xlabel('')
        else:
            ax.set_xlabel('$-t\\; [{\\rm GeV}^2]$', fontsize=20)
    axs[5].legend(handles, labels, loc="center", borderaxespad=0.).draw_frame(0)
    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig

def HERAF2Q2(path=None, fmt='png', **kwargs):
    """Makes plot of HERA F2(Q2) DIS data with fit lines and bands"""

    #title = 'H1 hep-ex/9603004 (Aid:1996au)'
    title = 'H1-F2'
    ymin, ymax = 0, 1.79
    fig, ax = plt.subplots(figsize=[8,6])
    panels = [ data.dset[208],
               data.dset[202] ]
    for npanel, panelset in enumerate(panels):
        panel(ax, points=panelset, xaxis='Q2', **kwargs)
        if npanel == 0:
            ax.text(20, 1.1, r'$x_B = {}$'.format(panelset[0].xB), fontsize=16)
            # Take legend info once
            handles, labels = ax.get_legend_handles_labels()
        else:
            ax.text(50, 0.63, r'$x_B = {}$'.format(panelset[0].xB), fontsize=16)

    ax.set_ylabel(r'$F_{2}^{p}(x_{\rm B}, Q^2)$', fontsize=16)
    ax.set_xlabel(r'$Q^2\; [{\rm GeV}^2]$', fontsize=16)
    plt.xlim(0.0, 65)
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(10))
    ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(5))
    ax.tick_params(axis='both', which='major', labelsize=14)
    plt.ylim(ymin, ymax)
    #ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
    #ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.02))
    ax.legend(handles, labels, loc="upper right", borderaxespad=0.).draw_frame(0)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig

def HERAF2xB(path=None, fmt='png', **kwargs):
    """Makes plot of HERA F2(xB) DIS data with fit lines and bands"""

    #title = 'H1 hep-ex/9603004 (Aid:1996au)'
    title = 'H1-F2'
    fig, ax = plt.subplots(figsize=[8,6])
    Q15 = data.select(DISpoints, criteria=['Q2 == 15'])
    Q8 = data.select(DISpoints, criteria=['Q2 == 8.5'])
    panels = [Q15, Q8 ]
    for npanel, panelset in enumerate(panels):
        panel(ax, points=panelset, xaxis='xB', **kwargs)
        ax.set_xscale('log')
        if npanel == 0:
            ax.text(0.0006, 1.3, r'$Q^2 = {:.0f}\, {{\rm GeV}}^2$'.format(panelset[0].Q2), fontsize=16)
            ax.text(0.00024, 0.82, r'${:.0f}\, {{\rm GeV}}^2$'.format(panelset[0].Q2), fontsize=16)
            # Take legend info once
            handles, labels = ax.get_legend_handles_labels()
        else:
            ax.text(0.00013, 1.05, r'${}\, {{\rm GeV}}^2$'.format(panelset[0].Q2), fontsize=16)
            ax.text(0.00018, 0.62, r'${}\, {{\rm GeV}}^2$'.format(panelset[0].Q2), fontsize=16)
    ax.set_ylabel(r'$F_{2}^{p}(x_{\rm B}, Q^2)$', fontsize=16)
    ax.set_xlabel(r'$x_B$', fontsize=16)
    plt.xlim(0.0001, 0.02)
    #ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(10))
    #ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(5))
    ax.tick_params(axis='both', which='major', labelsize=14)
    plt.ylim(0.3, 1.5)
    #ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
    #ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.02))
    ax.legend(handles, labels, loc="upper right", borderaxespad=0.).draw_frame(0)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig

def HERAF2(path=None, fmt='png', **kwargs):
    """Makes plot of HERA F2 DIS data with fit lines and bands"""

    #title = 'H1 hep-ex/9603004 (Aid:1996au)'
    title = 'H1-F2'
    lbl = r'$F_{2}^{p}(x_{\rm B}, Q^2)$'
    ymin, ymax = 0, 1.79
    fig, axs = plt.subplots(2, 1, sharey=True, sharex=True, figsize=[8,12])
    axs = axs.reshape(2)
    #fig.suptitle(title + ' -- ' +  observable)
    #
    panels = [ data.dset[208],
               data.dset[202] ]
    for npanel, panelset in enumerate(panels):
        ax = axs[npanel]
        panel(ax, points=panelset, xaxis='Q2', kinlabels=['xB'], **kwargs)
        plt.xlim(0.0, 65)
        #ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
        #ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.02))
        plt.ylim(ymin, ymax)
        #ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
        #ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.02))
        if npanel == 0:
                # Take legend info on this panel
                handles, labels = ax.get_legend_handles_labels()
        if (npanel % 2) == 13: #1:
            # Remove y-axis label from right panels
            ax.set_ylabel('')
        else:
            ax.set_ylabel(lbl, fontsize=16)
        if npanel < 0:
            # Remove x-axis label from upper panels
            ax.set_xlabel('')
        else:
            ax.set_xlabel(r'$Q^2\; [{\rm GeV}^2]$', fontsize=16)
    # Draw legend on this panel
    axs[0].legend(handles, labels, loc="upper right", borderaxespad=0.).draw_frame(0)
    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig

def CLASJ(path=None, fmt='png', **kwargs):
    """Makes plot of CLAS-14 ALU data with fit lines and bands.
       Some CLAS-07 data points are also plotted.
    """

    # old data
    aux = [(Q2, data.DataSet(grp)) for Q2, grp in
            itertools.groupby(data.dset[25], lambda pt: pt.Q2)]
    old = []
    for Q2, ds in aux:
        ds.__dict__ = data.dset[25].__dict__.copy()
        old.append((Q2, ds))
    # new data
    dataset = data.dset[85]
    title = ''
    title = 'CLAS (prelim. 2014)'
    fig, axs = plt.subplots(2, 2, sharey=True, sharex=True, figsize=[7,5])
    axs = axs.reshape(4)
    fig.suptitle(title)
    # Each different Q2 has its own panel
    panelorder = [3, 4, 1, 2]
    npanel = 0
    for Q2, grp in itertools.groupby(dataset, lambda pt: pt.Q2):
        npan = panelorder[npanel]
        ax = axs[npan-1]
        panelpoints = list(grp)
        # now transform list into DataSet instance ...
        panelset = data.DataSet(panelpoints)
        panelset.__dict__ = dataset.__dict__.copy()
        # ... and plot
        if npan == 2:
            # fake doubling of points to get line visibility
            ptd = panelset[0].copy()
            ptd.tm = ptd.tm+0.02
            panelset.append(ptd)
            panel(ax, points=panelset, xaxis='tm', kinlabels=['Q2', 'xB'], **kwargs)
        else:
            if npan == 1:
                panel(ax, points=[panelset, old[5][1]], xaxis='tm', kinlabels=['Q2', 'xB'], **kwargs)
            elif npan == 3:
                panel(ax, points=[panelset, old[2][1]], xaxis='tm', kinlabels=['Q2', 'xB'], **kwargs)
            else:
                panel(ax, points=panelset, xaxis='tm', kinlabels=['Q2', 'xB'], **kwargs)
        plt.xlim(0.0, 0.6)
        plt.ylim(0.0, 0.4)
        ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
        if npan == 1:
            ax.legend(loc='upper left', borderaxespad=0.).draw_frame(0)
        if (npan % 2) != 1:
            # Leave labels only on leftmost panels
            ax.set_ylabel('')
        else:
            ax.set_ylabel('$A_{LU}^{\\sin\\phi}$', fontsize=16)
        if npan < 3:
            ax.set_xlabel('')
        else:
            ax.set_xlabel('$-t\\; [{\\rm GeV}^2]$', fontsize=16)
        npanel += 1
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig


def CLAS08(path=None, fmt='png', **kwargs):
    """Makes plot of CLAS 2008 ALU data with fit lines and bands"""

    title = ''
    title = 'CLAS (Gavalian:2008aa)'
    alldata = data.select(data.dset[81], criteria=['FTn == -1'])
    datasets = [alldata[:3], alldata[-3:]]
    xaxes = ['Q2', 'tm']
    xlabels = ['$Q^2\\; [{\\rm GeV}^2]$', '$-t\\; [{\\rm GeV}^2]$']
    fig, axs = plt.subplots(1, 2, sharey=True, figsize=[6,3])
    fig.suptitle(title)
    # Each different Q2 has its own panel
    for npanel, ax in enumerate(axs):
        panel(ax, points=datasets[npanel], xaxis=xaxes[npanel], **kwargs)
        ax.set_xlabel(xlabels[npanel], fontsize=16)
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.2))
    axs[0].set_ylabel('$A_{LU}^{\\sin\\phi}$', fontsize=16)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig

def CLASTSA(path=None, fmt='png', **kwargs):
    """Plot CLAS TSA."""
    id = 54
    title = 'CLAS-TSA'
    fig = plt.figure()
    fig.suptitle(title)
    xaxes = ['tm', 'xB']
    # we have 2x3 points
    for x in range(2):
        ax = fig.add_subplot(2, 1, x+1)
        #ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))  # tickmarks
        panel(ax, points=data.dset[id][x*3:x*3+3], xaxis=xaxes[x], **kwargs)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig

def HallAFT(path=None, fmt='png', **kwargs):
    """Makes plot of harmonics of Hall-A data with fit lines"""

    subsets = {}
    subsets[1] = data.select(data.dset[50], criteria=['Q2 == 1.5', 'FTn == -1'])
    subsets[2] = data.select(data.dset[50], criteria=['Q2 == 1.9', 'FTn == -1'])
    subsets[3] = data.select(data.dset[50], criteria=['Q2 == 2.3', 'FTn == -1'])
    subsets[4] = data.select(data.dset[51], criteria=['FTn == 0'])
    subsets[5] = data.select(data.dset[51], criteria=['FTn == 1'])
    #subsets[6] = data.select(data.dset[51], criteria=['FTn == 2'])
    ylabels = 3*['$d^{\,4}\\Sigma^{\,\\sin\\phi,w}$'] + \
                ['$d^{\,4}\\sigma^{\,\\cos0\\phi,w}$'] + \
                ['$d^{\,4}\\sigma^{\,\\cos\\phi,w}$']
    title = 'MunozCamacho:2006hx'
    fig = plt.figure()
    fig.suptitle(title)
    Qs = ['1.5', '1.9', '2.3', '2.3', '2.3']
    for npanel in range(1,6):
        ax = fig.add_subplot(2,3,npanel)
        panel(ax, points=subsets[npanel], xaxis='tm', **kwargs)
        ax.set_ylabel('%s' % ylabels[npanel-1], fontsize=16)
        if npanel != 4:
            ax.text(-0.31, 0.017, '$Q^2\\!= %s\\,{\\rm GeV}^2$' % Qs[npanel-1], fontsize=12)
        else:
            ax.text(-0.31, 0.065, '$Q^2\\!= %s\\,{\\rm GeV}^2$' % Qs[npanel-1], fontsize=12)
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
    fig.subplots_adjust(wspace=0.7)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig


def HallA06(lines=None, path=None, fmt='png'):
    title = 'HallA06'
    NPTS = 12
    # tupls are (dataset id, ind_low, ind_high)
    xQbins = [1, 2, 3,
              'blank', 'blank', 4,
              'blank', 'blank', 5]
    subsets = {}
    subsets[1] = data.select(data.dset[50], criteria=['Q2 == 1.5', 'FTn == -1'])
    subsets[2] = data.select(data.dset[50], criteria=['Q2 == 1.9', 'FTn == -1'])
    subsets[3] = data.select(data.dset[50], criteria=['Q2 == 2.3', 'FTn == -1'])
    subsets[4] = data.select(data.dset[51], criteria=['FTn == 0'])
    subsets[5] = data.select(data.dset[51], criteria=['FTn == 1'])
    #fig, axs = plt.subplots(3, 3, sharey='row', sharex=True, figsize=[10,10])
    fig, axs = plt.subplots(3, 3, figsize=[7,9])
    axs = axs.reshape(9)
    for pn, iset in enumerate(xQbins):
        if iset == 'blank':
            axs[pn].axis('off')
            continue
        pts = pd.DataFrame([(pt.Q2, pt.xB, pt.tm, pt.val, pt.err) for pt in subsets[iset]], columns=('Q2', 'xB', 'tm', 'val', 'err'))
        axs[pn].errorbar(pts.tm.values, pts.val.values, pts.err.values, linestyle='None',
            capsize=2, color='black')
        if not isinstance(lines, list): lines = [lines]
        styles = ['r--', 'b-', 'g-.', 'b--', 'k:']
        for nl, th in enumerate(lines):
            xs, ys = thline(th, pts, subsets[iset][0], npts=NPTS)
            axs[pn].plot(xs, ys, styles[nl], label=th.name)
        axs[pn].set_xlim(0.12, 0.42)
        axs[pn].xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
        s1 = r'$Q^2 =\, {:.2f}$'.format(pts.Q2.mean())
        s1 += '\n'
        s1 += r'$x_B =\, {:.2f}$'.format(pts.xB.mean())
        if pn<3:
            axs[pn].text(0.25, 0.0175, s1, fontsize=14)
            axs[pn].set_ylim(0.0,0.024)
        if pn in [4,5]:
            #axs[pn].text(0.25, 0.06, s1, fontsize=16)
            axs[pn].set_ylim(0.0,0.088)
            axs[pn].yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.03))
        if pn in [7,8]:
            #axs[pn].text(0.25, 0.0005, s1, fontsize=16)
            axs[pn].set_ylim(-0.005, 0.017)
            axs[pn].axhline(y=0, linewidth=1, color='g')  # y=0 thin line
            axs[pn].yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.01))
        #else:
            #axs[pn].text(0.22, 0.05, s1, fontsize=16)
        if pn in [0,1,8]:
            axs[pn].set_xlabel(r'$-t\; [{\rm GeV}^2]$', fontsize=16)
        if pn == 5:
            axs[pn].set_ylabel(r'$d\sigma^{\cos 0\phi,w}$', fontsize=20)
        if pn == 8:
            axs[pn].set_ylabel(r'$d\sigma^{\cos \phi,w}$', fontsize=20)
            axs[pn].legend(loc=(-1.7,1.), fontsize=14, handlelength=3).draw_frame(0)
        if pn == 0:
            axs[pn].set_ylabel(r'$\Delta\sigma^{\sin\phi,w}$', fontsize=20)
        if pn in [1,2]:
            axs[pn].get_yaxis().set_visible(False)
    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig


def HallAphi(path=None, fmt='png', **kwargs):
    """Makes plot of Hall-A data phi-dependence with fit lines"""

    ids = [9, 14, 20, 21, 23, 24]
    subsets = {}
    subsets[9] = data.select(data.dset[33], criteria=['Q2 == 1.5', 't == -0.17'])
    subsets[14] = data.select(data.dset[33], criteria=['Q2 == 1.9', 't == -0.23'])
    subsets[20] = data.select(data.dset[33], criteria=['Q2 == 2.3', 't == -0.33'])
    subsets[21] = data.select(data.dset[34], criteria=['t == -0.17'])
    subsets[23] = data.select(data.dset[34], criteria=['t == -0.28'])
    subsets[24] = data.select(data.dset[34], criteria=['t == -0.33'])
    title = ''
    fig = plt.figure(figsize=[12,7])
    fig.suptitle(title)
    fig.subplots_adjust(wspace=0.2)
    for np in range(1,4):
        ax = fig.add_subplot(1,3,np)
        panel(ax, points=subsets[np], xaxis='phi', kinlabels=['Q2', 't'], **kwargs)
        #ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(2.))
        ax.set_xlabel('$\\phi\\; {\\rm [deg]}$', fontsize=16)
        ax.set_ylabel('{}'.format(subsets[np][0].observable), fontsize=18)
        if np == 1:
            ax.legend(loc='upper left', borderaxespad=0.).draw_frame(0)
        if np == 3:
            ax.set_ylim(0.0, 0.8)
        else:
            ax.set_ylim(-0.4, 0.4)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig

def CLAS15phi(path=None, fmt='png', **kwargs):
    """Makes plot of CLAS 15 data phi-dependence with fit lines"""

    title = 'CLAS-15-xs'
    fig, axs = plt.subplots(2, 3, sharey='row', sharex=True, figsize=[14,8])
    axs = axs.reshape(6)
    panels = [
       data.select(data.dset[98], criteria=['xB == 0.335', 'Q2 == 2.78', 't == -0.20']),
       data.select(data.dset[98], criteria=['xB == 0.335', 'Q2 == 2.78', 't == -0.26']),
       data.select(data.dset[98], criteria=['xB == 0.335', 'Q2 == 2.78', 't == -0.45']),
       data.select(data.dset[97], criteria=['xB == 0.335', 'Q2 == 2.78', 't == -0.20']),
       data.select(data.dset[97], criteria=['xB == 0.335', 'Q2 == 2.78', 't == -0.26']),
       data.select(data.dset[97], criteria=['xB == 0.335', 'Q2 == 2.78', 't == -0.45'])
               ]
    for np, panelset in enumerate(panels):
        ax = axs[np]
        panel(ax, points=panelset, xaxis='phi', kinlabels=['xB', 'Q2', 't'], **kwargs)
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(60))
        ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(6))
        #ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
        #ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.02))
        if np == 0:
            # Take legend info on this panel
            handles, labels = ax.get_legend_handles_labels()
            # Draw legend on this panel
            ax.legend(handles, labels, loc="upper center",
                borderpad=0.8, handlelength=3, fancybox=True, framealpha=1.0,
                #facecolor='navajowhite',
                prop=matplotlib.font_manager.FontProperties(size="x-large"))
            ax.set_ylabel(r'$d\sigma = d\sigma^{\leftarrow} + d\sigma^{\rightarrow}$', fontsize=20)
        elif np == 3:
            ax.set_ylabel(r'$\Delta\sigma = d\sigma^{\leftarrow} - d\sigma^{\rightarrow}$', fontsize=20)
        else:
            ax.set_ylabel('')
        if np < 3:
            # Remove x-axis label from upper panels
            ax.set_xlabel('')
        else:
            ax.set_xlabel(r'$\phi\; {\rm [deg]}$', fontsize=16)
        if np < 3:
            ax.set_ylim(-0.03, 0.3)
        else:
            ax.set_ylim(-0.1, 0.1)
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.tick_params(axis='both', which='minor', labelsize=14)
    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig

def CLAS15xs(lines=None, path=None, fmt='png'):
    title = 'CLAS15xs'
    NPTS = 12
    XLUw = pd.DataFrame([(pt.Q2, pt.xB, pt.tm, pt.val, pt.err) for pt in data.dataset[101]], columns=('Q2', 'xB', 'tm', 'val', 'err'))
    xQbins = [(0,4), (5,9), (11,15), (16,20), (21,26), (30,33), (35,38), (39,42), (43,45), (46,47)]
    fig, axs = plt.subplots(2,5, sharey='row', sharex=True, figsize=[12,9])
    axs = axs.reshape(10)
    for pn, (a,b) in enumerate(xQbins):
        pts = XLUw[a:b+1]
        axs[pn].errorbar(pts.tm.values, pts.val.values, pts.err.values, linestyle='None', color='black')
        if not isinstance(lines, list): lines = [lines]
        styles = ['r--', 'b-', 'g-.', 'b--', 'k:', 'r:', 'k--']
        for nl, th in enumerate(lines):
            xs, ys = thline(th, pts, data.dset[101][a], npts=NPTS)
            axs[pn].plot(xs, ys, styles[nl], label=th.name)
        axs[pn].set_xlim(0.05, 0.55)
        axs[pn].xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
        s1 = '$Q^2 = {:.2f}$'.format(pts.Q2.mean())
        s1 += '\n'
        s1 += '$x_B = {:.2f}$'.format(pts.xB.mean())
        if pn<5:
            axs[pn].text(0.2, 0.30, s1, fontsize=16)
        else:
            axs[pn].text(0.22, 0.045, s1, fontsize=16)
        if pn == 7:
            axs[pn].set_xlabel('$-t\\; [{\\rm GeV}^2]$', fontsize=16)
        if (pn % 5) == 0:
            axs[pn].set_ylabel(r'$\Delta\sigma^{\sin\phi,w}$', fontsize=24)
        if pn == 9:
            axs[pn].legend(loc=(-0.5,0.5),
                borderpad=0.8, handlelength=3, fancybox=True, framealpha=1.0,
                facecolor='navajowhite',
                prop=matplotlib.font_manager.FontProperties(size="x-large"))
        axs[pn].tick_params(axis='both', which='major', labelsize=12)
        axs[pn].tick_params(axis='both', which='minor', labelsize=12)
    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig


def HallA15(lines=None, enh=True, path=None, fmt='png'):
    title = 'HallA15'
    NPTS = 12
    # tupls are (dataset id, ind_low, ind_high)
    if enh:
        xQbins = [(117, 0, 4), (117, 5,9), (117, 10, 14),
            ('blank', 0, 0), (116, 0, 4), (116, 5, 9),
            ('blank', 0, 0), (116, 10, 14), (116, 15, 19)]
    else:
        xQbins = [(117, 0, 4), (117, 5,9), (117, 10, 14),
            ('blank', 0, 0), (123, 0, 4), (123, 5, 9),
            ('blank', 0, 0), (123, 10, 14), (123, 15, 19)]
    #fig, axs = plt.subplots(3, 3, sharey='row', sharex=True, figsize=[8,10])
    fig, axs = plt.subplots(3, 3, figsize=[7,9])
    axs = axs.reshape(9)
    for pn, (id, a,b) in enumerate(xQbins):
        if id == 'blank':
            axs[pn].axis('off')
            continue
        pts = pd.DataFrame([(pt.Q2, pt.xB, pt.tm, pt.val, pt.err) for pt in data.dset[id][a:b+1]], columns=('Q2', 'xB', 'tm', 'val', 'err'))
        axs[pn].errorbar(pts.tm.values, pts.val.values, pts.err.values, linestyle='None',
            capsize=2, color='black')
        if not isinstance(lines, list): lines = [lines]
        styles = ['r--', 'b-', 'g-.', 'b--', 'k:']
        for nl, th in enumerate(lines):
            xs, ys = thline(th, pts, data.dset[id][a], npts=NPTS)
            axs[pn].plot(xs, ys, styles[nl], label=th.name)
        axs[pn].set_xlim(0.12, 0.42)
        axs[pn].xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
        s1 = r'$Q^2 =\, {:.2f}$'.format(pts.Q2.mean())
        s1 += '\n'
        s1 += r'$x_B =\, {:.2f}$'.format(pts.xB.mean())
        if pn < 3:
            axs[pn].text(0.25, 0.015, s1, fontsize=14)
            axs[pn].set_ylim(0.0,0.02)
        if pn in [4,5]:
            #axs[pn].text(0.25, 0.06, s1, fontsize=16)
            axs[pn].set_ylim(0.0,0.088)
            axs[pn].yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.03))
        if pn in [7,8]:
            #axs[pn].text(0.25, 0.0005, s1, fontsize=16)
            axs[pn].set_ylim(-0.005, 0.017)
            axs[pn].axhline(y=0, linewidth=1, color='g')  # y=0 thin line
            axs[pn].yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.01))
        #else:
            #axs[pn].text(0.22, 0.05, s1, fontsize=16)
        if pn in [0,7,8]:
            axs[pn].set_xlabel(r'$-t\; [{\rm GeV}^2]$', fontsize=14)
        if pn == 4:
            axs[pn].set_ylabel(r'$d\sigma^{\cos 0\phi,w}$', fontsize=20)
        if pn == 7:
            axs[pn].set_ylabel(r'$d\sigma^{\cos \phi,w}$', fontsize=20)
            axs[pn].legend(loc=(-1.4,1.), fontsize=14, handlelength=3).draw_frame(0)
        if pn == 0:
            axs[pn].set_ylabel(r'$\Delta\sigma^{\sin\phi,w}$', fontsize=20)
        if pn in [1,2,5,8]:
            axs[pn].get_yaxis().set_visible(False)
    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig


def HallA17(observable='XUUcos0', lines=None, path=None, fmt='png'):
    title = 'HallA (2017) '+observable
    NPTS = 24
    #   E    Q2   ind_low  ind_high
    #  5.6  1.5   0  2
    #  3.4  1.5   3  5
    #  5.5  1.75  6  9
    #  4.5  1.75  10 13
    #  5.6  2     14 17
    #  4.5  2     18 21
    # tupls are (dataset id, ind_low, ind_high)
    xQbins = [(136, 14, 17), (136, 6, 9), (136, 0, 2),
          (136, 18, 21), (136, 10, 13), (136, 3, 5)]
    if observable == 'XUUcos0':
        ylab = r'$d\sigma^{\cos 0\phi,w}$'
        yminU, ymaxU = 0.0, 0.099   # upper row of panels
        yTU = 0.006
        ytickU = 0.02
        yminD, ymaxD = 0.0, 0.24   # lower row of panels
        yTD = 0.01
        ytickD = 0.05
        dshift = 0
        pshift = 0
    elif observable == 'XUUcos1':
        ylab = r'$d\sigma^{\cos \phi,w}$'
        yminU, ymaxU = -0.005, 0.06   # upper row of panels
        yTU = 0.0006
        ytickU = 0.02
        yminD, ymaxD = -0.005, 0.06   # lower row of panels
        yTD = -0.003
        ytickD = 0.02
        dshift = 0
        pshift = 22
    elif observable == 'XLUsin1':
        ylab = r'$\Delta\sigma^{\sin\phi,w}$'
        yminU, ymaxU = -0.00, 0.019   # upper row of panels
        yTU = 0.0006
        ytickU = 0.005
        yminD, ymaxD = -0.005, 0.047   # lower row of panels
        yTD = -0.003
        ytickD = 0.01
        dshift = -1
        pshift = 0
        xQbins[3] = ('blank', 0, 0)
    else:
        raise ValueError('Observable %s unavailable.' % observable)
    #fig, axs = plt.subplots(2, 3, sharey='row', sharex='all', figsize=[12,10])
    fig, axs = plt.subplots(2, 3, figsize=[12,10])
    axs = axs.reshape(6)
    for pn, (id, a, b) in enumerate(xQbins):
        if id == 'blank':
            axs[pn].axis('off')
            continue
        pts = pd.DataFrame([(pt.in1energy, pt.Q2, pt.xB, pt.tm, pt.val, pt.err) for
            pt in data.dset[id+dshift][a+pshift:b+pshift+1]], columns=('E', 'Q2', 'xB', 'tm', 'val', 'err'))
        axs[pn].errorbar(pts.tm.values, pts.val.values, pts.err.values, linestyle='None',
        capsize=2, color='black', label='Hall A 2017')
        if not isinstance(lines, list): lines = [lines]
        styles = ['r--', 'b-', 'g-.', 'b--', 'k:']
        for nl, th in enumerate(lines):
            xs, ys = thline(th, pts, data.dset[id+dshift][a+pshift], npts=NPTS)
            axs[pn].plot(xs, ys, styles[nl], label=th.name)
        axs[pn].set_xlim(0.12, 0.42)
        axs[pn].xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
        s1 = r'$Q^2 =\, {:.2f}$'.format(pts.Q2.mean())
        s1 += '\n'
        s1 += r'$E =\, {:.2f}$'.format(pts.E.mean())
        axs[pn].axhline(y=0, linewidth=0.5, linestyle='-.', color='g')  # y=0 thin line
        if pn < 3:
            # upper row of panels
            axs[pn].set_ylim(yminU,ymaxU)
            axs[pn].yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(ytickU))
            axs[pn].text(0.25, yTU, s1, fontsize=14)
        else:
            # lower row of panels
            axs[pn].set_xlabel(r'$-t\; [{\rm GeV}^2]$', fontsize=14)
            axs[pn].set_ylim(yminD,ymaxD)
            axs[pn].yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(ytickD))
            axs[pn].text(0.25, yTD, s1, fontsize=14)
        if pn in [0,3] or (observable == 'XLUsin1' and pn == 4):
            axs[pn].set_ylabel(ylab, fontsize=20)
        if observable[:6] == 'XUUcos' and pn in [1,2,4,5]:
            axs[pn].get_yaxis().set_visible(False)
        if observable == 'XLUsin1' and pn in [1,2,5]:
            axs[pn].get_yaxis().set_visible(False)
        if pn == 0:
            axs[pn].legend(loc=2, fontsize=14, handlelength=3).draw_frame(1)
    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig


def H1ZEUS(path=None, fmt='png', **kwargs):
    """Makes plot of H1 DVCS data with fit lines"""

    subsets = {}
    subsets[1] = [data.select(data.dset[63], criteria=['Q2 == 8.']),
            data.select(data.dset[63], criteria=['Q2 == 15.5']),
            data.select(data.dset[63], criteria=['Q2 == 25.'])]
    subsets[2] = [data.dset[46]] # ZEUS t-dep
    subsets[3] = [data.dset[48],
            data.select(data.dset[41], criteria=['Q2 == 8.']),
            data.select(data.dset[41], criteria=['Q2 == 15.5']),
            data.select(data.dset[41], criteria=['Q2 == 25.'])]
    subsets[4] = [data.dset[47]] # ZEUS Q2-dep
    xs = ['t', 't', 'W', 'Q2']
    title = 'H1_07_ZEUS_08_DVCS'
    fig, (axu, axl) = plt.subplots(2,2,figsize=[14,8])
    #fig.suptitle(title)
    #fig.subplots_adjust(bottom=0.1, hspace=0.3)
    #fig.subplots_adjust(hspace=0.4)
    for npanel in range(1,5):
        if npanel < 3:
            ax = axu[npanel-1]
        else:
            ax = axl[npanel-3]
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
            ax.text(-0.78, 22, '${\\rm H1}$', fontsize=16)
            ax.text(-0.78, 10, '${\\rm W = 82}\\, {\\rm GeV}$', fontsize=12)
            ax.text(-0.78, 3, '$Q^2\\!= 8,\\, 15.5,\\, 25\\,{\\rm GeV}^2$', fontsize=12)
        if npanel==2:
            ax.text(-0.3, 6, '${\\rm ZEUS}$', fontsize=16)
            ax.text(-0.3, 3.5, '${\\rm W = 104}\\, {\\rm GeV}$', fontsize=12)
            ax.text(-0.3, 2.1, '$Q^2\\!= 3.2\\,{\\rm GeV}^2$', fontsize=12)
            handles, labels = ax.get_legend_handles_labels()
            handles.pop()
            labels.pop()
            ax.legend(handles, labels, loc="upper left").draw_frame(0)
        if npanel==3:
            pass
            #ax.text(50, 28, '${\\rm ZEUS}\\, (idem):$', fontsize=16)
            #ax.text(50, 5, '${\\rm H1}\\, (idem):$', fontsize=16)
        else: # npanel==4
            pass
            #ax.text(30, 3, '${\\rm ZEUS}\\, (idem)$', fontsize=16)
            #ax.set_xlim(0, 80)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig

def DVMP(H109WdepXL, path=None, fmt='png', **kwargs):
    """Makes plot of H1 DVMP data with fit lines"""

    subsets = {}
    subsets[1] = [data.dset[76]]
    subsets[2] = [data.select(H109WdepXL, criteria=['Q2 == 6.6']),
                  data.select(H109WdepXL, criteria=['Q2 == 11.9']),
                  data.select(H109WdepXL, criteria=['Q2 == 19.5']),
                  data.select(H109WdepXL, criteria=['Q2 == 35.6'])]
    subsets[3] = []
    subsets[4] = []
    xs = ['Q2', 'W', 't', 'Q2']
    title = 'H1 DVMP'
    fig = plt.figure(figsize=(14,8))
    fig.suptitle(title)
    fig.subplots_adjust(bottom=0.1, hspace=0.3)
    for npanel in range(1,3):
        ax = fig.add_subplot(2,2,npanel)
        ax.set_yscale('log')  # y-axis to be logarithmic
        panel(ax, points=subsets[npanel], xaxis=xs[npanel-1], **kwargs)
        #if npanel < 3:
        #    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.2))
        # y labels
        if npanel==1:
            ax.set_ylabel('$\\sigma_L\\quad [{\\rm nb}]$', fontsize=18)
        elif npanel==3:
            ax.set_ylabel('$d\\sigma/dt\\quad [{\\rm nb/GeV}^2]$', fontsize=18)
        else: # npanel 4
            ax.set_ylabel('')
        # x labels
        if npanel==1:
            ax.set_xlabel('$Q^2\\quad [{\\rm GeV}^2]$', fontsize=18)
        elif npanel==2:
            ax.set_xlabel('$W\\quad [{\\rm GeV}]$', fontsize=18)
        #else: # npanel 4
        #    ax.set_xlabel('$t\\quad [{\\rm GeV}^2]$', fontsize=18)
        if npanel==1:
            ax.text(15, 100, '${\\rm W = 75}\\, {\\rm GeV}$', fontsize=16)
            ax.legend(loc='upper right', borderaxespad=0.).draw_frame(0)
        if npanel==2:
            ax.text(60, 300, '$Q^2 = 6.6,\\, 11.9,\\, 19.5,\\, 35.6\\, {\\rm GeV}^2$', fontsize=16)
        #if npanel==3:
        #    ax.text(50, 5, '${\\rm H1}\\, (idem):$', fontsize=16)
        #else: # npanel==4
        #    ax.text(30, 3, '${\\rm ZEUS}\\, (idem)$', fontsize=16)
        #    #ax.set_xlim(0, 80)
        ax.set_ylim(0.3, 1.e3)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig

def COMPASSt(path=None, fmt='png', **kwargs):
    """Makes plot of COMPASS BCS asymmetry in terms of t for various values of xB"""

    title = ''#Prediction for COMPASS BCSA'
    fig = plt.figure()
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
    ks['s'] = 2 * constants.Mp * ks['in1energy'] + constants.Mp2
    ks['phi'] = 3.141
    ks['units'] = {'phi' : 'radian'}
    ks['frame'] = 'BMK'
    ks['observable'] = 'BCSA'
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
        ax.text(0.03, 0.35, "%s = %s" % (constants.toTeX['xB'], xB), fontsize=14)
        ax.text(0.03, 0.3, "%s = %s" % (constants.toTeX['Q2'], Q2), fontsize=14)
        if pn == 1:
            #ax.legend(bbox_to_anchor=(0.9, 0.7), loc='upper right',
            #        borderaxespad=0.).draw_frame(0)
            pass
        if (pn != 1) and (pn != 4):
            # no y tick labels on right panels
            ax.set_yticklabels([])
        else:
            # y label on left panels
            ax.set_ylabel(constants.toTeX['BCSA'], fontsize=20)
        if pn <= 3:
            # no x tick labels on upper panels
            ax.set_xticklabels([])
        if pn >= 4:
            # x-label only on lower panels
            ax.set_xlabel(constants.toTeX['tm'], fontsize=17)
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
        #fig.show()
    return fig

#################################################################
##                                                             ##
##  [3]  Other plots                                           ##
##                                                             ##
#################################################################

def jbod(path=None, fmt='png', **kwargs):
    """Plot "just a bunch of data" with lines

    """
    title = 'just-a-bunch-of-data'
    fig = plt.figure()
    fig.suptitle(title)
    ax = fig.add_subplot(1, 1, 1)
    n = 1
    if isinstance(kwargs['points'], data.DataSet):
        kwargs['points'] = [kwargs['points']]
    for pts in kwargs['points']:
        label = ' + '.join(set('%s-%s-%s' % (pt.collaboration,
            str(pt.year)[-2:], pt.observable) for pt in pts))
        # x-coordinate is in data frame, while y-coordinate is in 0.0-1.0 frame
        ax.text(n, 0.9, label, transform=matplotlib.transforms.blended_transform_factory(
            ax.transData, ax.transAxes))
        for pt in pts:
            pt.npt = n
            n += 1
    panel(ax, xaxis='npt', **kwargs)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig

def CFF(cffs=['ImH', 'ReH'], path=None, fmt='png', **kwargs):
    """Makes plots of cffs given by various theories/models

    cffs    -- List of CFFs to be plotted. Each produces two panels.
        FIXME: gepard CFFs via DR are NOT plotted correctly!!!

    """
    title = 'CFF'
    fig, axs = plt.subplots(len(cffs), 2, sharex='col', sharey='row', figsize=(12,3.5*len(cffs)))
    #fig.suptitle(title)
    # Define abscissas
    xmin = 0.1
    xmax = 0.25
    xvals = np.linspace(xmin, xmax, 20) # left panel
    logxvals = np.logspace(-2.3, -0.01, 20)  # right panel
    # ranges of y axis
    ylims = {
              #'ImH': (-4.8, 35), 'ReH': (-3, 2),
              'ImH': (-5, 20), 'ReH': (-4, 2),
             'ImE': (-10, 0), 'ReE': (-40, 40),
             'ImEt': (-50, 50), 'ReEt': (-50, 50),
             #'ImHt': (-4, 8), 'ReHt': (-12, 12)}
             'ImHt': (-1, 3), 'ReHt': (-5, 10)}
    # Plot panels
    for n in range(len(cffs)):
        cff = cffs[n]
        # [LEFT:] data range
        if len(cffs) == 1:
            ax = axs[0]
        else:
            ax = axs[n, 0]
        panel(ax, xaxis='xi', xs=xvals, kins={'observable':cff, 't':-0.2, 'Q2':4.,
        'units':{'CFF': 1}, 'observable': 'CFF'}, **kwargs)
        ax.axhline(y=0, linewidth=0.5, color='k')  # y=0 thin line
        try:
            ax.set_ylim(*ylims[cff])
            ax.set_ylabel(constants.toTeX['%s' % cff], fontsize=20)
        except KeyError:
            pass
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.tick_params(axis='both', which='minor', labelsize=14)
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.05))  # tickmarks
        if n == 0:
            ax.legend(loc='upper right')
            ax.legend().draw_frame(0)
            ax.text(0.2, 0.9, r"$t = -0.2\,\, {\rm GeV}^2$", transform=ax.transAxes,
                    fontsize=12)
        if n == len(cffs)-1:
            ax.set_xlabel(constants.toTeX['xi'], fontsize=16) # [RIGHT:] extrapolation range (logarithmic)
        if len(cffs) == 1:
            ax = axs[1]
        else:
            ax = axs[n, 1]
        ax.set_xscale('log')  # x-axis to be logarithmic
        panel(ax, xaxis='xi', xs=logxvals, kins={'observable':cff, 't':-0.2, 'Q2':4.,
            'units':{'CFF': 1}, 'observable': 'CFF'}, **kwargs)
        ax.axvspan(0.1, 0.25, facecolor='g', alpha=0.1)  # vertical band
        ax.axhline(y=0, linewidth=0.5, color='k')  # y=0 thin line
        #ax.set_ylim(*ylims[cff])
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.tick_params(axis='both', which='minor', labelsize=14)
        if n == len(cffs)-1:
            ax.set_xlabel(constants.toTeX['xi'], fontsize=16)
    fig.subplots_adjust(wspace=0.04, hspace=0.06)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig

def CFF2(cffs=['ImH', 'ReH'], path=None, fmt='png', **kwargs):
    """Makes plots of cffs at two t given by various theories/models

    cffs    -- List of CFFs to be plotted. Each produces two panels.

    """
    title = ''
    fig = plt.figure()
    #fig.suptitle(title)
    colors = ['red', 'brown']     # worm human colors :-)
    nncolors = ['blue', 'green']  # cold computer colors
    linestyles = ['solid', 'dashed']
    # Define abscissas
    ximin = 0.01
    ximax = 0.3
    logxvals = np.logspace(np.log10(ximin), np.log10(ximax), 40)
    # ordinates
    #ylims = {'ImH': (-4.3, 35), 'ReH': (-6.5, 8)}
    ylims = {'ImH': (-4.3, 35), 'ReH': (-3.5, 2),
             'ImE': (-40, 35), 'ReE': (-25, 5),
             'ImEt': (-50, 100), 'ReEt': (-50, 100),
             'ImHt': (-6.9, 6.9), 'ReHt': (-10, 20)}
    # Plot panels
    ts = [-0.3, -0.12]
    for n in range(len(cffs)):
        for nt in range(len(ts)):
            cff = cffs[n]
            # all-x logarithmic
            ax = fig.add_subplot(len(cffs), 2, 2*n+nt+1)
            ax.set_xscale('log')  # x-axis to be logarithmic
            panel(ax, xaxis='xi', xs=logxvals, kins={
                    'observable':cff, 't':ts[nt], 'Q2':4., 'units':{'CFF': 1}, 'observable': 'CFF'}, **kwargs)
            ax.axhspan(-0.0005, 0.0005, facecolor='g', alpha=0.6)  # horizontal bar
            ax.text(0.31, 0.02, "$t = %s\\, {\\rm GeV}^2$" % str(ts[nt]),
                    transform=ax.transAxes, fontsize=14)
            ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%s'))
            if nt == 0:
                # plot data-region vertical band for left panels
                #ax.axvspan(0.025, 0.136, facecolor='g', edgecolor='black', alpha=0.2)
                ax.set_ylabel(constants.toTeX['%s' % cff], fontsize=16)
            if nt == 1:
                # no y tick labels on right panels
                ax.set_yticklabels([])
            if n < len(cffs)-1:
                # no x tick labels on upper panels
                ax.set_xticklabels([])
            if n == len(cffs)-1:
                # x-label only on lower panels
                ax.set_xlabel(constants.toTeX['xixB'], fontsize=14)
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
            ax.set_ylim(*ylims[cff])
            ax.set_xlim(ximin, ximax)
            #ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.02))  # tickmarks
            for label in ax.get_xticklabels() + ax.get_yticklabels():
                label.set_fontsize(12)
    fig.subplots_adjust(bottom=0.2, wspace=0.0, hspace=0.0)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig

def CFFt(cffs=['ImH', 'ReH'], path=None, fmt='png', **kwargs):
    """Makes plots of cffs as function of t given by various theories/models

    cffs    -- List of CFFs to be plotted.

    """
    title = ''
    fig = plt.figure(figsize=(12,4*len(cffs)))
    #fig.suptitle(title)
    # Define abscissas
    tmvals = np.linspace(0.05, 0.55, 10) # right panel
    # ordinates
    #ylims = {'ImH': (-4.3, 35), 'ReH': (-6.5, 8),
    # Plot panels
    for n in range(len(cffs)):
        cff = cffs[n]
        # smaller x
        ax = fig.add_subplot(len(cffs), 2, 2*n+1)
        panel(ax, xaxis='tm', xs=tmvals, kins={'observable':cff, 'xi':0.12, 'Q2':4.,
           'units':{'CFF': 1}, 'observable': 'CFF'}, **kwargs)
        ax.set_xlabel(constants.toTeX['tm'], fontsize=15)
        try:
            ax.set_ylabel(constants.toTeX['%s' % cff], fontsize=18)
        except KeyError:
            pass
        ax.axhline(y=0, linewidth=0.5, color='g')  # y=0 thin line
        if n == 0:
            ax.legend(loc='upper right')
            ax.legend().draw_frame(0)
            ax.text(0.1, 0.18, r"$\xi = 0.12$", transform=ax.transAxes,
                    fontsize=12)
            # ax.text(0.1, 0.12, "$Q^2 = 4\\, {\\rm GeV}^2$", transform=ax.transAxes,
                    # fontsize=12)
        # larger x
        ax = fig.add_subplot(len(cffs), 2, 2*n+2)
        panel(ax, xaxis='tm', xs=tmvals, kins={'observable':cff, 'xi':0.22, 'Q2':4.,
            'units':{'CFF': 1}, 'observable': 'CFF'}, **kwargs)
        ax.set_xlabel(constants.toTeX['tm'], fontsize=15)
        try:
            ax.set_ylabel(constants.toTeX['%s' % cff], fontsize=18)
        except KeyError:
            pass
        ax.axhline(y=0, linewidth=0.5, color='g')  # y=0 thin line
        if n == 0:
            ax.legend(loc='upper right')
            ax.legend().draw_frame(0)
            ax.text(0.1, 0.18, r"$\xi = 0.22$", transform=ax.transAxes,
                    fontsize=12)
            # ax.text(0.1, 0.12, "$Q^2 = 4\\, {\\rm GeV}^2$", transform=ax.transAxes,
                    # fontsize=12)
    fig.subplots_adjust(bottom=0.1)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig

def CFF3(cffs=['ImH', 'ReH', 'ImE', 'ReE'],
        path=None, fmt='png', **kwargs):
    """Return plot of cffs given by various theories/models

    cffs    -- List of even number of CFFs to be plotted.

    """
    title = 'CFF3'
    nrows = int(len(cffs)/2)
    fig, axs = plt.subplots(nrows, 2, figsize=[12, nrows*3], sharex='col')
    colors = ['red', 'brown']     # worm human colors :-)
    nncolors = ['blue', 'green']  # cold computer colors
    linestyles = ['solid', 'dashed']
    # Define abscissas
    xvals = np.linspace(0.05, 0.3, 100)
    # ordinates
    for pn, ax in enumerate(axs.flatten()):
        cff = cffs[pn]
        panel(ax, xaxis='xi', xs=xvals, kins={'observable':cff, 't':-0.2, 'Q2':4.,
            'units':{'CFF': 1}, 'observable': 'CFF'}, **kwargs)
        ax.set_ylabel(constants.toTeX['{}'.format(cff)], fontsize=18)
        ax.axhline(y=0, lw=0.5, color='k', ls=':')  # horizontal bar
        for ticklabel in ax.get_xticklabels() + ax.get_yticklabels():
                    ticklabel.set_fontsize(16)
        if pn >= 2*(nrows-1):
            ax.set_xlabel(constants.toTeX['xixB'], fontsize=16)
        if pn == 0:
            leg = ax.legend(loc='upper right', handlelength=4.0, fancybox=True)
            frame  = leg.get_frame()
            #frame.set_facecolor('0.90')    # set the frame face color to light gray
            for t in leg.get_texts():
                t.set_fontsize(16)    # the legend text fontsize
            for l in leg.get_lines():
                l.set_linewidth(2.0)  # the legend line width
        if cff == 'ImHt':
            ax.set_ylim(0,3.8)
            ax.yaxis.set_major_locator( matplotlib.ticker.MultipleLocator(1.)) 
            ax.yaxis.set_minor_locator( matplotlib.ticker.MultipleLocator(0.2)) 
        elif cff == 'ImH':
            ax.yaxis.set_major_locator( matplotlib.ticker.MultipleLocator(3.)) 
            ax.yaxis.set_minor_locator( matplotlib.ticker.MultipleLocator(1.)) 
        elif cff == 'ReE':
            ax.yaxis.set_major_locator( matplotlib.ticker.MultipleLocator(10.)) 
            ax.yaxis.set_minor_locator( matplotlib.ticker.MultipleLocator(2.)) 
        elif cff == 'ImEt':
            ax.yaxis.set_major_locator( matplotlib.ticker.MultipleLocator(20.)) 
            ax.yaxis.set_minor_locator( matplotlib.ticker.MultipleLocator(5.)) 
        elif cff == 'ImE':
            ax.set_ylabel(constants.toTeX['{}'.format(cff)], fontsize=18, labelpad=-8)

    fig.subplots_adjust(hspace=0., wspace=0.25)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig


def CFF3log(cffs=['ImH', 'ReH', 'ImE', 'ReE', 'ImHt', 'ImEt'], tval=-0.2,
        path=None, fmt='png', **kwargs):
    """Return plot of cffs given by various theories/models.

    cffs    -- List of even number of CFFs to be plotted.

    """
    title = 'CFF3log'
    nrows = int(len(cffs)/2)
    fig, axs = plt.subplots(nrows, 2, figsize=[12, nrows*3], sharex='col')
    colors = ['red', 'brown']     # worm human colors :-)
    nncolors = ['xkcd:teal', 'xkcd:tomato']  # cold computer colors
    linestyles = ['solid', 'dashed']
    # Define abscissas
    xvals = np.logspace(-2.3, -0.01, 50)
    # ordinates
    for pn, ax in enumerate(axs.flatten()):
        cff = cffs[pn]
        ax.set_xscale('log')  # x-axis to be logarithmic
        panel(ax, xaxis='xi', xs=xvals, kins={'observable':cff, 't':tval, 'Q2':4.,
            'units':{'CFF': 1}, 'observable': 'CFF'}, **kwargs)
        ax.set_ylabel(constants.toTeX['{}'.format(cff)], fontsize=18)
        ax.axhline(y=0, lw=0.5, color='k', ls=':')  # horizontal bar
        for ticklabel in ax.get_xticklabels() + ax.get_yticklabels():
                    ticklabel.set_fontsize(16)
        if pn >= 2*(nrows-1):
            ax.set_xlabel(constants.toTeX['xixB'], fontsize=16)
        if pn == 0:
            leg = ax.legend(loc='upper right', handlelength=4.0, fancybox=True)
            frame  = leg.get_frame()
            #frame.set_facecolor('0.90')    # set the frame face color to light gray
            for t in leg.get_texts():
                t.set_fontsize(16)    # the legend text fontsize
            for l in leg.get_lines():
                l.set_linewidth(2.0)  # the legend line width
        if cff == 'ImHt':
            # ax.set_ylim(0,3.8)
            ax.yaxis.set_major_locator( matplotlib.ticker.MultipleLocator(1.)) 
            ax.yaxis.set_minor_locator( matplotlib.ticker.MultipleLocator(0.2)) 
        elif cff == 'ImH':
            ax.yaxis.set_major_locator( matplotlib.ticker.MultipleLocator(3.)) 
            ax.yaxis.set_minor_locator( matplotlib.ticker.MultipleLocator(1.)) 
            ax.text(0.1, 0.1, r'$t = {}\, {{\rm GeV}}^2$'.format(tval),
                    transform=ax.transAxes, fontsize=16)
        elif cff == 'ReE':
            ax.yaxis.set_major_locator( matplotlib.ticker.MultipleLocator(10.)) 
            ax.yaxis.set_minor_locator( matplotlib.ticker.MultipleLocator(2.)) 
        elif cff == 'ImEt':
            ax.yaxis.set_major_locator( matplotlib.ticker.MultipleLocator(20.)) 
            ax.yaxis.set_minor_locator( matplotlib.ticker.MultipleLocator(5.)) 
        elif cff == 'ImE':
            ax.set_ylabel(constants.toTeX['{}'.format(cff)], fontsize=18, labelpad=-8)

    fig.subplots_adjust(hspace=0., wspace=0.25)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        #fig.show()
    return fig

