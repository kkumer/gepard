#!/usr/bin/env python

"""Plotting kinematic coverage of data."""

# run it from ipy via 'run data/polys'


#import shelve
import sys, string, os
import copy

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

#import Data, Model, Approach, utils
from abbrevs import *


######################################################################
# Convex Hull code by  Dinu C. Gherman
# http://code.activestate.com/recipes/66527-finding-the-convex-hull-of-a-set-of-2d-points/
######################################################################

def _myDet(p, q, r):
    """Calc. determinant of a special matrix with three 2D points.

    The sign, "-" or "+", determines the side, right or left,
    respectivly, on which the point r lies, when measured against
    a directed vector from p to q.
    """

    # We use Sarrus' Rule to calculate the determinant.
    # (could also use the Numeric package...)
    sum1 = q[0]*r[1] + p[0]*q[1] + r[0]*p[1]
    sum2 = q[0]*p[1] + r[0]*q[1] + p[0]*r[1]

    return sum1 - sum2


def _isRightTurn((p, q, r)):
    "Do the vectors pq:qr form a right turn, or not?"

    assert p != q and q != r and p != r
            
    if _myDet(p, q, r) < 0:
	return 1
    else:
        return 0


def _isPointInPolygon(r, P):
    "Is point r inside a given polygon P?"

    # We assume the polygon is a list of points, listed clockwise!
    for i in xrange(len(P[:-1])):
        p, q = P[i], P[i+1]
        if not _isRightTurn((p, q, r)):
            return 0 # Out!        

    return 1 # It's within!

def convexHull(P):
    "Calculate the convex hull of a set of points."

    # Get a local list copy of the unique points and sort them lexically.
    unique = {}
    for p in P:
        unique[p] = 1

    points = unique.keys()
    points.sort()

    # Build upper half of the hull.
    upper = [points[0], points[1]]
    for p in points[2:]:
	upper.append(p)
	while len(upper) > 2 and not _isRightTurn(upper[-3:]):
	    del upper[-2]

    # Build lower half of the hull.
    points.reverse()
    lower = [points[0], points[1]]
    for p in points[2:]:
	lower.append(p)
	while len(lower) > 2 and not _isRightTurn(lower[-3:]):
	    del lower[-2]

    # Remove duplicates.
    del lower[0]
    del lower[-1]

    # Concatenate both halfs and return.
    return tuple(upper + lower)


collider = [('H1', H1all), ('ZEUS', ZEUSall)]
shifted = []
for pt in HallAall:
    pta = pt.copy()
    pta.xB = pta.xB+0.005
    shifted.append(pta)
    
fixed = [('CLAS', CLASall), ('HERMES', HERMESall), ('Hall A', HallAall+shifted)]

fig = plt.figure()

colors = ['red', 'blue', 'green', 'purple', 'black']
markers = ['o', '+', 'v', '+', 's']

k = 0
# collider
ax = fig.add_subplot(121)
ax.set_xscale('log')
ax.set_yscale('log')
for label, set in collider:
    points = [(pt.xB, pt.Q2) for pt in set]
    hull = convexHull(points)
    ax.add_patch(Polygon(hull, closed=True, fill=True, 
        facecolor=colors[k], alpha=0.4, label=label))
    ax.scatter([x for x,y in points], [y for x,y in points], 
            facecolor=colors[k],marker=markers[k])
    k += 1
#ax3.set_xlim((0,6))
#ax3.set_ylim((0,2.5))
ax.legend(loc='upper left',
        prop=matplotlib.font_manager.FontProperties(size="larger")).draw_frame(0)
ax.set_xlabel('$x_B$', fontsize=20)
ax.set_ylabel('$Q^2\\,[{\\rm GeV}^2]$', fontsize=20)
ax.text(3e-4, 1.2, "collider", fontsize=24)
ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%3i'))
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontsize(14)
ax.autoscale_view()
# fixed target
ax = fig.add_subplot(122)
#ax.set_xscale('log')
#ax.set_yscale('log')
for label, set in fixed:
    points = [(pt.xB, pt.Q2) for pt in set]
    hull = convexHull(points)
    ax.add_patch(Polygon(hull, closed=True, fill=True, 
        facecolor=colors[k], alpha=0.5, label=label))
    ax.scatter([x for x,y in points], [y for x,y in points], 
            facecolor=colors[k],marker=markers[k])
    k += 1
ax.legend(prop=matplotlib.font_manager.FontProperties(size="larger")).draw_frame(0)
ax.set_xlabel('$x_B$', fontsize=20)
#ax.set_ylabel('$Q^2\\,[{\\rm GeV}^2]$', fontsize=20)
ax.text(0.15, 0.3, "fixed targets", fontsize=24)
ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%2i'))
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontsize(14)
ax.autoscale_view()
fig.subplots_adjust(bottom=0.5, hspace=0.1)
plt.show()
#fig.savefig('data.pdf', format='pdf')

#def LesHouches(path=None, fmt='png'):
    #"""Make Gepard-Pegasus comparison plot."""
    #title = 'GepVsPeg'
    #fig = plt.figure()
    #fig.canvas.set_window_title(title)
    #fig.suptitle(title)
    #fig.subplots_adjust(bottom=0.15)
    #ax = fig.add_subplot(1,1,1)
    #ax.set_xscale('log')
    #ax.set_yscale('log')

    #pt = Data.DummyPoint()
    #pt.Q2 = 1.0e4
    #pt.t = 0
    #colors = ['red', 'blue', 'green', 'black']
    #styles = [':', '--', '-']
    #lbl = ['LO', 'NLO', 'NNLO']
    #for flav in range(4):
        #m = Model.ComptonGepard(ansatz=ansaetze[flav], 
                #fftype=fftypes[flav], q02=2.0)
        #t = Approach.hotfixedBMK(m)
        #t.m.g.parint.nf = 4
        #asp0 = 0.35 / 2./ np.pi
        #t.m.g.astrong.asp = np.array([asp0, asp0, asp0])
        #t.m.g.astrong.mu02 = 2.0
        #t.m.g.mbcont.phi = 1.9
        ##t.m.g.mbcont.c = 0.25
        #t.m.g.parint.pid = -2
        #t.m.g.parint.acc = 4
        #for ord in range(3):
            #if ((flav >=2) and (ord == 2)): break
            #t.m.g.parint.p = ord
            #t.m.g.init()
            #xF = []
            #for xB in xBs:
                #pt.xi = xB
                #t.m.g.newcall = 1
                #if flav == 0:
                    #xF.append(m.gpdHzeroG(pt))
                #else:
                    #xF.append(pt.xi*m.gpdHzeroQ(pt))
            #err = [abs((peg-gep)/peg) for peg,gep in zip(ssr[ord,:,6-flav], xF)]
            #if np.any([np.isnan(x) for x in err]):
                #ax.plot([], [], color=colors[flav], linestyle=styles[ord],
                        #linewidth=2, label='NAN: %s %s' % (lbl[ord], flavlbl[flav])) 
            #else:
                #ax.plot(xBs, err, color=colors[flav], linestyle=styles[ord],
                             #linewidth=2, label='%s %s' % (lbl[ord], flavlbl[flav])) 
    #ax.set_xlim(1.0e-7, 1.0)
    ##ax.text(0.2, 8., "$x_B = %s$" % xB, fontsize=14)
    ##ax.set_ylabel('$\\Im\\! m \\mathcal{H}(x_{\\rm B}, t, Q^2=2\\,\
    ##         {\\rm GeV}^2)/\\pi$', fontsize=10)
    ## Lower pannels with Re(H)
    ##ax.legend(prop=matplotlib.font_manager.FontProperties(size="smaller")).draw_frame(0)
    #ax.legend().draw_frame(0)
    #if path:
        #fig.savefig(os.path.join(path, title+'.'+fmt), format=fmt)
    #else:
        #fig.canvas.draw()
        #fig.show()
    #return fig


#if __name__ == '__main__':
    #draw_poly(path='.', fmt='pdf')
