#!/usr/bin/env python

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
import models
import Approach
from constants import *
from plots import *
from fits import *

# [1] Load data

data = utils.loaddata()   # dictionary {1 : dataset, ...}

# choose datapoints for fitting

#fitpoints = data[1] + data[8] + data[29]  # DM's GLO set
#fitpoints = data[1] + data[5] + data[25]  # KK's set
fitpoints = data[2][:3] + data[8][:3]  # test set
#fitpoints = data[1] + data[8] + data[29] + data[30]  # DM's GLO1 set
#fitpoints = data[2] + data[5] + data[25] + data[26] + data[27] + data[28]  # KK's global set
#fitpoints = data[1] + data[2] + data[8] + data[22] + data[29] + data[30] 


# [2] Choose theoretical approach

ff = models.FormFactors()
#bmk = Approach.BMK()
b = Approach.hotfixedBMK(ff)

# Transform data into conventions of approach used, and 
# pre-calculate what can be pre-calculated

## Prepare just fitpoints ...
# [pt.prepare(b) for pt in fitpoints]
## ... or prepare ALL available datapoints
[[pt.prepare(b) for pt in set] for set in data.values()]

# [3] Define function to be minimized i.e. chi-square

def fcn(NS, alS, alpS, M2S, rS, bS, Nv, alv, alpv, M2v, rv, bv, C, M2C, tNv, tM2v, trv, tbv):
   pars = (NS, alS, alpS, M2S, rS, bS, Nv, alv, alpv, M2v, rv, bv, C, M2C, tNv, tM2v, trv, tbv)
   chisq = 0.
   for pt in fitpoints:
       chisq = chisq + (
               (getattr(b, pt.yaxis)(pt, pars) - pt.val)**2 / pt.err**2 )
   return chisq

def fcn2(*args):
    pars = {'NS': args[0], 'alS':args[1]}
    chisq = 0.
    for pt in fitpoints:
        chisq = chisq + (
                (getattr(b, pt.yaxis)(pt, pars) - pt.val)**2 / pt.err**2 )
    return chisq

# [3] Minimization

# Fitting parameters:
varpars = ['NS', 'alS']

# Constructing function with variable number of parameters by
# string manipulation.
sargs = ", ".join(varpars)
m = Minuit(eval("lambda %s: fcn2(%s)" % (sargs, sargs)))

# Parameters initialization
for par in varpars:
    m.values[par] = ff.H[par]

#m.tol = 20.0
#m.strategy = 2
m.printMode = 0
m.maxcalls = 4000

m.migrad()

print "ncalls = ", m.ncalls
#fitpars = tuple([m.values[key] for key in m.parameters])
dof = len(fitpoints) - utils.npars(m)
fitprob = (1.-gammainc(dof/2., m.fval/2.)) # probability of this chi-sq
fitres = (m.fval, dof, fitprob)
print 'P(chi-square, d.o.f) = P(%1.2f, %2d) = %5.4f' % fitres

sigmas = [(getattr(b, pt.yaxis)(pt, m.values) - pt.val) / pt.err for
            pt in fitpoints]
utils.prettyprint(sigmas)


