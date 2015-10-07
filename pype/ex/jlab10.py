#!/usr/bin/env python

""" 
Plotting NN CFFs for JLAB 2010 proceedings

Published in arXiv:1008.2762

"""

import os, shelve
import numpy as np

import matplotlib
if os.sys.platform == 'win32':
    matplotlib.use('WxAgg')
else: #linux
    matplotlib.use('TkAgg')
import pylab as plt


import Data, utils
from constants import toTeX, Mp2, Mp



def _axband(ax, tm, xvals, fun, color='g', avg=True):
    """Make a band of x*fun defined by ndarray-valued fun."""
    pt = Data.DummyPoint()
    pt.t = -tm
    pt.xB = xvals
    pt.xi = pt.xB/(2.-pt.xB)
    up = []
    down = []
    #res = xvals * fun(pt)
    res = fun(pt) / 3.14159
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

def CFF(t, ImH=True, path=None, fmt='png', average=True):
    """Makes plot of cffs given by neural network
    
    t     -- Neural network 'theory'
    ImH   -- plot 2-panel plot with ImH as well
    average -- plot 1-sigma band (othewise plot minimum-maximum band)

    """
    old = t.model.parameters['nnet']
    t.model.parameters['nnet'] = 'ALL'
    colors = ['red', 'blue', 'brown', 'purple']
    linestyles = ['solid', 'dashed']
    tms = [0.0, 0.3]
    # Define abscissas
    xvals = np.linspace(0.05, 0.4, 50) # right panel
    # ordinatas

    # Ploto ImH if requested
    if ImH:
      fig = plt.figure(figsize=[16./1.5, 5./1.3])
      ax = fig.add_subplot(1, 2, 1)
      for tm in tms:
          _axband(ax, tm, xvals, getattr(t.model, 'ImH'), 
                  color=colors[tms.index(tm)], avg=average)
      ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
      ax.set_xlabel('$x_{\\rm Bj}$', fontsize=15)
      #ax.set_ylabel('%s / $\\pi$' % cff, fontsize=18)
      ax.set_ylabel('$\\Im{\\rm m}\\, {\\cal H}(x_{\\rm Bj}, t) / \\pi$', fontsize=18)
      apply(ax.set_ylim, (-2.0, 4.0))
      ax.set_xlim(0.05, 0.4)
      ax.text(0.065, 0.5, "$t=0$", fontsize=14)
      ax.text(0.065, -1.5, "$t=-0.3\,{\\rm GeV}^2$", fontsize=14)
      #ax.text(0.03, -0.27, "data region only", fontsize=14)
      #ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.02))  # tickmarks
      ax = fig.add_subplot(1, 2, 2)
    else: 
      fig = plt.figure(figsize=[8./1.5, 5./1.3])
      ax = fig.add_subplot(1, 1, 1)

    # Always plot ReH
    for tm in tms:
        _axband(ax, tm, xvals, getattr(t.model, 'ReH'), 
                color=colors[tms.index(tm)], avg=average)
    ax.axhline(y=0, linewidth=1, color='g')  # y=0 thin line
    ax.set_xlabel('$x_{\\rm Bj}$', fontsize=15)
    #ax.set_ylabel('%s / $\\pi$' % cff, fontsize=18)
    ax.set_ylabel('$\\Re{\\rm e}\\, {\\cal H}(x_{\\rm Bj}, t) / \\pi$', fontsize=18)
    apply(ax.set_ylim, (-2.0, 4.0))
    ax.set_xlim(0.05, 0.4)
    ax.text(0.065, 0.5, "$t=0$", fontsize=14)
    ax.text(0.065, -1.5, "$t=-0.3\,{\\rm GeV}^2$", fontsize=14)
    #ax.text(0.03, -0.27, "data region only", fontsize=14)
    #ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.02))  # tickmarks

    # left  = 0.125  # the left side of the subplots of the figure
    # right = 0.9    # the right side of the subplots of the figure
    # bottom = 0.0   # the bottom of the subplots of the figure
    # top = 0.9      # the top of the subplots of the figure
    # wspace = 0.2   # the amount of width reserved for blank space between subplots
    # hspace = 0.2   # the amount of height reserved for white space between subplots
    fig.subplots_adjust(bottom=0.12, left=0.08, right=0.98, top=0.95)

    t.model.parameters['nnet'] = old
    if path:
        fig.savefig(os.path.join(path, 'jlab10-GLO.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig


if __name__ == '__main__':
    # Load fit results and choose one
    db = shelve.open('theories.db')
    tNN = db['EX-11-15-2']
    
    # plot a fig
    CFF(tNN, ImH=True, path='.', fmt='eps')

    # calculate some points
    tNN.model.parameters['nnet'] = 'ALL'
    pt = Data.DummyPoint()
    pt.t = -0.28

    for xB in [0.1, 0.2, 0.3]:
        pt.xB = xB
        pt.xi = pt.xB / (2. - pt.xB)
        res = tNN.m.ImH(pt)/np.pi
        print 'Im(CFF_H(xB = %s))/pi = %.2f +- %.2f' %(str(xB), res.mean(), res.std())

