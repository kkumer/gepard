#!/usr/bin/env python

# This just to try commiting ...

import sys, os
import matplotlib

if os.sys.platform == 'win32':
    from minuit import *
    matplotlib.use('WxAgg')
else: #linux
    from minuit2 import Minuit2 as Minuit
    matplotlib.use('TkAgg')

from scipy.special import gammainc


import utils 
import Approach
from constants import *
from plots import *

# [1] Load data

data = utils.loaddata()   # dictionary {1 : dataset, ...}

# choose datapoints for fitting

fitpoints = data[1] + data[8] + data[29]  # DM's GLO set
#fitpoints = data[1] + data[5] + data[25]  # KK's set
#fitpoints = data[1][:3] + data[8][:3]  # test set
#fitpoints = data[1] + data[8] + data[29] + data[30]  # DM's GLO1 set


# [2] Choose theoretical approach

bmk = Approach.BMK()
b = Approach.hotfixedBMK()

# Transform data into conventions of approach used, and 
# pre-calculate what can be pre-calculated

## Prepare just fitpoints ...
# [pt.prepare(b) for pt in fitpoints]
## ... or prepare ALL available datapoints
[[pt.prepare(b) for pt in set] for set in data.values()]

# [3] Define function to be minimized i.e. chi-square

def fcn(NS, alS, alpS, M2S, sS, bS, Nv, alv, alpv, M2v, sv, bv, C, M2C):
   pars = (NS, alS, alpS, M2S, sS, bS, Nv, alv, alpv, M2v, sv, bv, C, M2C)
   chisq = 0.
   for pt in fitpoints:
       chisq = chisq + (
               (getattr(b, pt.yaxis)(pt, pars) - pt.val)**2 / pt.err**2 )
   return chisq


# [3] Minimization

m = Minuit(fcn, NS=1.5,    fix_NS=True,
                alS=1.13,  fix_alS=True,
                alpS=0.15, fix_alpS=True, 
                M2S=0.5,   fix_M2S=True,
                sS=0.0,    fix_sS=True,
                bS=43.9, 
                Nv=1.35,   fix_Nv=True,
                alv=0.43,  fix_alv=True,
                alpv=0.85, fix_alpv=True, 
                M2v=0.64,  fix_M2v=True,
                sv=0.535, # fix_sv=True,
                bv=2.3,   # fix_bv=True,
                C=-5.3,   # fix_C=True,
                M2C=0.227**2#, fix_M2C=True
                )

#m.tol = 0.1
#m.printMode = 1

m.migrad()

print "ncalls = ", m.ncalls
fitpars = tuple([m.values[key] for key in m.parameters])
dof = len(fitpoints) - utils.npars(m)
fitprob = (1.-gammainc(dof/2., m.fval/2.)) # probability of this chi-sq
fitres = (m.fval, dof, fitprob)
print 'P(chi-square, d.o.f) = P(%1.2f, %2d) = %5.4f' % fitres

sigmas = [(getattr(b, pt.yaxis)(pt, fitpars) - pt.val) / pt.err for
            pt in fitpoints]
utils.prettyprint(sigmas)


# DM's GLO fit to BCAcos1 and ALUIsin1 (HERMES) and ALUa (CLAS)
dmpars= (1.5, 1.13, 0.15, 0.5, 0.0, 43.89782184615788,
         1.35, 0.43, 0.85, 0.64, 0.5351586374852629, 2.3013417610626226,
         -5.301310593853695, 0.22676583042331164**2)
