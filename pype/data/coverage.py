#!/usr/bin/env python

"""Plotting kinematic coverage of data."""

# run it from ipy via 'run data/coverage'


#import shelve
import sys, string, os
import copy

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

#import Data, Model, Approach, utils
from abbrevs import *

RASTER = False
rorder = -3000
if RASTER:
    rorder = 1

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
    
fixed = [
         ('HERMES', HERMESall), 
         ('CLAS', CLASall+data[90]+data[97]), 
         ('Hall A', data[109]+data[110]+data[111]+HallAall+shifted)]

# New temporary for NN fits
HA15pts = H_BSDwpts + H_BSSw0pts + H_BSSw1pts
HA17pts = H17_BSDwpts + H17_BSSw0pts + H17_BSSw1pts
fixed = [
         ('HERMES', BCApts+ALUIpts+LPpoints[::2][:6]+TPpoints), 
         ('CLAS', C_BSDwpts[::2]+ C_BSSw0pts[::2] + C_BSSw1pts[::2]), 
         ('Hall A', HA17pts[::2]+HA15pts)]

fig = plt.figure(figsize=(8,4))

colors = ['red', 'blue', 'green', 'purple', 'orange']
#markers = ['o', '+', 'v', '+', 's']
markers = ['.', '.', '.', '.', '.']
styles = ['solid', 'dashed', 'solid', 'dashed', 'dotted']

ordinate = sys.argv[1]

k = 0
# collider
ax = fig.add_subplot(121)
ax.set_rasterization_zorder(rorder)
ax.set_xscale('log')
if ordinate == 'Q2':
    ax.set_yscale('log')
for label, set in collider:
    points = []
    for pt in set:
        try: 
            points.append((pt.xB, getattr(pt, ordinate)))
        except AttributeError:
            pass
    hull = convexHull(points)
    ax.add_patch(Polygon(hull, closed=True, fill=True, 
        ls=styles[k], lw=2.2,
        facecolor=colors[k], alpha=0.5, label=label, zorder=0))
    ax.scatter([x for x,y in points], [y for x,y in points], 
            facecolor=colors[k],marker=markers[k], s=0.5, zorder=1)
    k += 1
ax.legend(loc=(0.05,0.7),
        prop=matplotlib.font_manager.FontProperties(size="larger")).draw_frame(0)
ax.set_xlabel('$x_B$', fontsize=14)
if ordinate == 't':
    ax.set_ylabel('$t\\,[{\\rm GeV}^2]$', fontsize=14)
    ax.text(1.5e-4, 0.05, "collider", fontsize=12)
else:
    ax.set_ylabel('$Q^{2}\\,[{\\rm GeV}^2]$', fontsize=14)
    ax.text(1.5e-5, 75., "collider", fontsize=12)
    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%3i'))
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontsize(12)
ax.autoscale_view()
# fixed target
ax = fig.add_subplot(122)
ax.set_rasterization_zorder(rorder)
for label, set in fixed:
    points = [(pt.xB, getattr(pt, ordinate)) for pt in set]
    hull = convexHull(points)
    ax.add_patch(Polygon(hull, closed=True, fill=True,
        ls=styles[k], lw=2.2,
        facecolor=colors[k], alpha=0.5, label=label, zorder=0))
    ax.scatter([x for x,y in points], [y for x,y in points], 
            facecolor=colors[k],marker=markers[k], s=0.5, zorder=1)
    k += 1
l = ax.legend(loc=(0.48,0.65), 
        #prop=matplotlib.font_manager.FontProperties(size="larger")
             )
l.draw_frame(0)
#l.set_zorder(rorder)  # to rasterize, but you don't want that
ax.set_xlabel('$x_B$', fontsize=14)
#ax.set_ylabel('$Q^2\\,[{\\rm GeV}^2]$', fontsize=12)
if ordinate == 't':
    ax.text(0.03, 0.1, "fixed targets", fontsize=12)
else:
    ax.text(0.03, 6.5, "fixed targets", fontsize=12)
    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%2i'))
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontsize(12)
ax.autoscale_view()
fig.tight_layout()
#fig.subplots_adjust(wspace=0.4)
#plt.show()
#fig.savefig('coverage.png', format='png')
#fig.savefig('coverage.eps', format='eps')   # and set RASTER=True above
fig.savefig('coverage.pdf', format='pdf')
# And then rasterize  with pdftops -eps -r 600 coverage.pdf coverage.eps


#if __name__ == '__main__':
    #draw_poly(path='.', fmt='pdf')
