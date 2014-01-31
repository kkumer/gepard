#!/usr/bin/env python

"""Plot comparison with Pegasus PDF evolution routine."""

# run it from ipy via 'run ex/LesHouches'


import shelve
import sys, os, math

import numpy as np

import Data, Model, Approach, utils

import matplotlib
#if os.sys.platform == 'win32':
    #matplotlib.use('WxAgg')
#else: #linux
    #matplotlib.use('TkAgg')
import pylab as plt


from results import *


ansaetze = ['HOUCHE', 'HOUCHE', 'NSPHOU', 'NSMHOU']
fftypes = ['SINGLET', 'SINGLET', 'NONSINGLET', 'NONSINGLET']
flavlbl = ['gluon', 'singlet Q', 'C=1 NS Q', 'C=-1 NS Q']


ss = np.loadtxt('ex/LesHouches-input.dat'); Q2 = 2.
#ss = np.loadtxt('ex/LesHouches-1e2.dat'); Q2 = 1.e2
#ss = np.loadtxt('ex/LesHouches.dat'); Q2 = 1.e4
ssr = ss.reshape((3,21,7))

# ... which is formatted as:
# 0     1     2       3               4             5    6
# xB   Q2   ORD   x(u_v-d_v)   x(u+ubar-d-dbar)   xSIG  xg

xBs = [x for x,p in ssr[0,:,0:3:2]] 

def LesHouches(path=None, fmt='png'):
    """Make Gepard-Pegasus comparison plot."""
    title = 'GepVsPeg'
    suptitle = 'GeParD vs. Pegausus @ %.0f GeV^2' % Q2
    fig = plt.figure()
    #fig.canvas.set_window_title(title)
    fig.suptitle(suptitle)
    fig.subplots_adjust(bottom=0.15)
    ax = fig.add_subplot(1,1,1)
    ax.set_xscale('log')
    ax.set_yscale('log')

    pt = Data.DummyPoint()
    pt.Q2 = Q2
    pt.t = 0
    colors = ['red', 'blue', 'green', 'black']
    styles = [':', '--', '-']
    lbl = ['LO', 'NLO', 'NNLO']
    for flav in range(4):
        m = Model.ComptonGepard(ansatz=ansaetze[flav], 
                fftype=fftypes[flav], q02=2.0)
        t = Approach.hotfixedBMK(m)
        t.m.g.parint.nf = 4
        asp0 = 0.35 / 2./ np.pi
        t.m.g.astrong.asp = np.array([asp0, asp0, asp0])
        t.m.g.astrong.mu02 = 2.0
        t.m.g.mbcont.phi = 1.9
        #t.m.g.mbcont.c = 0.25
        t.m.g.parint.pid = -2
        t.m.g.parint.acc = 4
        for ord in range(3):
            if ((flav >=2) and (ord == 2)): break
            t.m.g.parint.p = ord
            t.m.g.init()
            xF = []
            for xB in xBs:
                pt.xi = xB
                t.m.g.newcall = 1
                if flav == 0:
                    xF.append(m.gpdHzeroG(pt))
                else:
                    xF.append(pt.xi*m.gpdHzeroQ(pt))
            err = [abs((peg-gep)/peg) for peg,gep in zip(ssr[ord,:,6-flav], xF)]
            if np.any([np.isnan(x) for x in err]):
                ax.plot([], [], color=colors[flav], linestyle=styles[ord],
                        linewidth=2, label='NAN: %s %s' % (lbl[ord], flavlbl[flav])) 
            else:
                ax.plot(xBs, err, color=colors[flav], linestyle=styles[ord],
                             linewidth=2, label='%s %s' % (lbl[ord], flavlbl[flav])) 
    ax.set_xlim(1.0e-7, 1.0)
    #ax.text(0.2, 8., "$x_B = %s$" % xB, fontsize=14)
    #ax.set_ylabel('$\\Im\\! m \\mathcal{H}(x_{\\rm B}, t, Q^2=2\\,\
    #         {\\rm GeV}^2)/\\pi$', fontsize=10)
    # Lower pannels with Re(H)
    #ax.legend(prop=matplotlib.font_manager.FontProperties(size="smaller")).draw_frame(0)
    ax.legend().draw_frame(0)
    if path:
        fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    else:
        fig.canvas.draw()
        fig.show()
    return fig


if __name__ == '__main__':
    LesHouches(path='.', fmt='pdf')
