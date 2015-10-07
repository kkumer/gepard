#!/usr/bin/env python

"""Plot of two NPB models for M. Guidal."""

# run it from ipy via 'run ex/HMG'


import shelve
import sys, os, math

import numpy as np

import Model, Approach, Fitter
import utils 

import matplotlib
if os.sys.platform == 'win32':
    matplotlib.use('WxAgg')
else: #linux
    matplotlib.use('TkAgg')
import pylab as plt


import Data, Approach, utils

from results import *


## Models
# Gepard only
mGepard = Model.ComptonGepard(cutq2=0.5)
tGepard = Approach.hotfixedBMK(mGepard)


# DR only
mDRonly = Model.ModelDR()
tDR = Approach.hotfixedBMK(mDRonly)
tDR.name = '(H1/ZEUS)+HERMES+CLAS'
tDR.m.parameters.update(DMepsGLO)

mDRonly1 = Model.ModelDR()
tDR1 = Approach.hotfixedBMK(mDRonly1)
tDR1.name = '(H1/ZEUS)+HERMES+CLAS+HallA'
tDR1.m.parameters.update(DMepsGLO1)

mDRPPsea = Model.ComptonModelDRPPsea()
m = Model.HybridDipole(mGepard, mDRPPsea)
t = Approach.BM10(m)
t.name = 'prelim. H1/ZEUS+HERMES+CLAS+HallA'
g = t.m.g
t.m.parameters.update(KKunp5)


def HMG(f, theories=[], path=None, fmt='png'):
    """Makes plot of cffH for M. Guidal's proceedings.

    Puts grids into file f.
    
    """
    title = 'HMG'
    fig = plt.figure()
    fig.canvas.set_window_title(title)
    fig.suptitle(title)
    fig.subplots_adjust(bottom=0.15)
    pt = Data.DummyPoint()
    pt.Q2 = 2.
    #ax.set_xscale('log')  # x-axis to be logarithmic
    mts = np.linspace(0.0, 0.6, 20) 
    colors = ['green', 'blue', 'purple']
    styles = ['--', '-', ':']
    lbl = ['without HallA', 'with Hall A']
    pn = 0
    for xB in [0.36, 0.25, 0.09]:
        pt.xB = xB
        pt.xi = xB/(2.-xB)
        # Upper pannels with Im(H)
        ln = 0
        ax = fig.add_subplot(2,3,pn+1)
        for th in theories:
            f.write('\n#  Im CFFH(xB = %s)/pi for fit %s\n' % (xB, lbl[ln]))
            f.write('# -t  Im(CFFH)/pi \n# ----------------\n')
            # kludge alert!
            ImHvals = []
            ReHvals = []
            for mt in mts:
                pt.mt = mt
                pt.t = -mt
                if th.model.__dict__.has_key('Gepard'): th.m.g.newcall = 1
                res = th.model.ImH(pt)/np.pi  # M. Guidal conventions!
                f.write('%5.3f   %5.3f\n' % (mt, res))
                ImHvals.append(res)
            ax.plot(mts, ImHvals, color=colors[ln], linestyle=styles[ln],
                         linewidth=3, label=th.name) 
            # printing to file
            ln += 1
        ax.set_ylim(-2.0, 11.0)
        ax.text(0.2, 8., "$x_B = %s$" % xB, fontsize=14)
        if pn == 0:
            ax.set_ylabel('$\\Im\\! m \\mathcal{H}(x_{\\rm B}, t, Q^2=2\\,\
                     {\\rm GeV}^2)/\\pi$', fontsize=10)
        # Lower pannels with Re(H)
        ln = 0
        ax = fig.add_subplot(2,3,pn+4)
        for th in theories:
            f.write('\n# -Re CFFH(xB = %s) for fit %s\n' % (xB, lbl[ln]))
            f.write('# -t  -Re(CFFH) \n# ----------------\n')
            # kludge alert!
            ImHvals = []
            ReHvals = []
            for mt in mts:
                pt.mt = mt
                pt.t = -mt
                if th.model.__dict__.has_key('Gepard'): th.m.g.newcall = 1
                res = -th.model.ReH(pt)       # M. Guidal conventions!
                ReHvals.append(res)
                f.write('%5.3f   %5.3f\n' % (mt, res))
            ax.plot(mts, ReHvals, color=colors[ln], linestyle=styles[ln], 
                            linewidth=3, label=lbl[ln]) 
            ln += 1
        ax.set_ylim(-2.0, 11.0)
        ax.set_xlabel('$-t \\, [{\\rm GeV}^2]$', fontsize=10)
        if pn == 0:
            ax.set_ylabel('$-\\Re\\! e \\mathcal{H}(x_{\\rm B}, t, Q^2=2\\,\
                     {\\rm GeV}^2)$', fontsize=10)
            ax.legend(prop=matplotlib.font_manager.FontProperties(size="smaller")).draw_frame(0)
        pn += 1
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig


f = open('grids.dat', 'w')
HMG(f, [tDR, tDR1], path='.', fmt='pdf')
f.close()
