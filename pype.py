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
from fits import *

# [1] Load data

data = utils.loaddata()   # dictionary {1 : dataset, ...}

# choose datapoints for fitting

fitpoints = data[1] + data[8] + data[29]  # DM's GLO set
#fitpoints = data[1] + data[5] + data[25]  # KK's set
#fitpoints = data[2][:3] + data[8][:3]  # test set
#fitpoints = data[1] + data[8] + data[29] + data[30]  # DM's GLO1 set
#fitpoints = data[2] + data[5] + data[25] + data[26] + data[27] + data[28]  # KK's global set
#fitpoints = data[1] + data[2] + data[8] + data[22] + data[29] + data[30] 


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

def fcn(NS, alS, alpS, M2S, rS, bS, Nv, alv, alpv, M2v, rv, bv, C, M2C, tNv, tM2v, trv, tbv):
   pars = (NS, alS, alpS, M2S, rS, bS, Nv, alv, alpv, M2v, rv, bv, C, M2C, tNv, tM2v, trv, tbv)
   chisq = 0.
   for pt in fitpoints:
       chisq = chisq + (
               (getattr(b, pt.yaxis)(pt, pars) - pt.val)**2 / pt.err**2 )
   return chisq

# [3] Minimization

# GENERIC
m = Minuit(fcn,
  NS = 1.5,      fix_NS = True,
 alS = 1.13,    fix_alS = True,
alpS = 0.15,   fix_alpS = True, 
 M2S = 0.5,     fix_M2S = True,
  rS = 1.0,      fix_rS = True,
  bS = 3.0,      fix_bS = False,
  Nv = 1.35,     fix_Nv = True,
 alv = 0.43,    fix_alv = True,
alpv = 0.85,   fix_alpv = True, 
 M2v = 0.64,    fix_M2v = False,
  rv = 1.5,      fix_rv = False,
  bv = 2.0,      fix_bv = False,    limit_bv = (0., 200.),
   C = 5.0,       fix_C = False,
 M2C = 0.7,     fix_M2C = False,
 tNv = 0.6,     fix_tNv = True,
tM2v = 0.64,   fix_tM2v = False,
 trv = 1.0,     fix_trv = False,
 tbv = 3.0,     fix_tbv = False,    limit_tbv = (0., 200.)
)

# DMGLO1
m = Minuit(fcn,
  NS = 1.5,         fix_NS = True,
 alS = 1.13,       fix_alS = True,
alpS = 0.15,      fix_alpS = True, 
 M2S = 0.5,        fix_M2S = True,     
  rS = 1.0,         fix_rS = True,
  bS = 2.00203,     fix_bS = False,      limit_bS = (1.5, 2.5),
  Nv = 1.35,        fix_Nv = True,
 alv = 0.43,       fix_alv = True,
alpv = 0.85,      fix_alpv = True, 
 M2v = 1.01097**2,   fix_M2v = False,    limit_M2v = (0.2, 5.),
  rv = 0.496383,      fix_rv = False,    limit_rv = (0.3, 0.7),
  bv = 2.15682,       fix_bv = False,    limit_bv = (1.2, 2.8),
   C = 6.90484,        fix_C = False,    limit_C = (6., 8.),
 M2C = 1.33924**2,   fix_M2C = False,    limit_M2C = (0., 15.),
 tNv = 0.6,          fix_tNv = True,
tM2v = 2.69667**2,  fix_tM2v = False,    limit_tM2v = (4., 9.),
 trv = 5.97923,      fix_trv = False,    limit_trv = (5., 7.),
 tbv = 3.25607,      fix_tbv = False,    limit_tbv = (2., 4.)
)

# SCAN
m = Minuit(fcn,
  NS = 1.5,         fix_NS = True,
 alS = 1.13,       fix_alS = True,
alpS = 0.15,      fix_alpS = True, 
 M2S = 0.5,        fix_M2S = True,     
  rS = 1.0,         fix_rS = True,
  bS = 2.0,         fix_bS = False,      limit_bS = (0.5, 50.),
  Nv = 1.35,        fix_Nv = True,
 alv = 0.43,       fix_alv = True,
alpv = 0.85,      fix_alpv = True, 
 M2v = 0.8,          fix_M2v = False,    limit_M2v = (0.2, 5.),
  rv = 0.5,           fix_rv = False,    limit_rv = (0.0, 5.0),
  bv = 4.0,           fix_bv = False,    limit_bv = (0.5, 25.),
   C = 5.,             fix_C = False,    limit_C = (0.0, 15.),
 M2C = 3.0,          fix_M2C = False,    limit_M2C = (0.2, 5.),
 #tNv = 0.6,          fix_tNv = True,
 tNv = 0.0,          fix_tNv = True,
tM2v = 2.69667**2,  fix_tM2v = True,    limit_tM2v = (0.2, 15.),
 trv = 5.97923,      fix_trv = True,    limit_trv = (0., 5.),
 tbv = 3.25607,      fix_tbv = True,    limit_tbv = (0., 25.)
)

# DMGLOB
m = Minuit(fcn,
  NS = 1.5,         fix_NS = True,
 alS = 1.13,       fix_alS = True,
alpS = 0.15,      fix_alpS = True, 
 M2S = 0.5,        fix_M2S = True,     
  rS = 1.0,         fix_rS = True,
  bS = 3.19,        fix_bS = False,      limit_bS = (3.1, 3.3),
  Nv = 1.35,        fix_Nv = True,
 alv = 0.43,       fix_alv = True,
alpv = 0.85,      fix_alpv = True, 
 M2v = 0.99,         fix_M2v = False,    limit_M2v = (0.9, 1.1),
  rv = 0.71,          fix_rv = False,    limit_rv = (0.6, 0.8),
  bv = 0.5,           fix_bv = True,    limit_bv = (0.4, 0.6),
   C = 1.49,           fix_C = False,    limit_C = (1.4, 1.5),
 M2C = 1.192,       fix_M2C = False,    limit_M2C = (1.1, 1.3),
 tNv = 0.0,          fix_tNv = True,
tM2v = 2.69667**2,  fix_tM2v = True,    limit_tM2v = (4., 9.),
 trv = 5.97923,      fix_trv = True,    limit_trv = (5., 7.),
 tbv = 3.25607,      fix_tbv = True,    limit_tbv = (2., 4.)
)


    
# DMGLO
m = Minuit(fcn,
  NS = 1.5,         fix_NS = True,
 alS = 1.13,       fix_alS = True,
alpS = 0.15,      fix_alpS = True, 
 M2S = 0.5,        fix_M2S = True,     
  rS = 1.0,         fix_rS = True,
  bS = 2.25,        fix_bS = False,      limit_bS = (3.1, 3.3),
  Nv = 1.35,        fix_Nv = True,
 alv = 0.43,       fix_alv = True,
alpv = 0.85,      fix_alpv = True, 
 M2v = 0.683**2,         fix_M2v = False,    limit_M2v = (0.9, 1.1),
  rv = 0.684,          fix_rv = False,    limit_rv = (0.6, 0.8),
  bv = 0.5,           fix_bv = True,    limit_bv = (0.4, 0.6),
   C = 1.12,           fix_C = False,    limit_C = (1.4, 1.5),
 M2C = 1.22*2,       fix_M2C = False,    limit_M2C = (1.1, 1.3),
 tNv = 0.0,          fix_tNv = True,
tM2v = 2.69667**2,  fix_tM2v = True,    limit_tM2v = (4., 9.),
 trv = 5.97923,      fix_trv = True,    limit_trv = (5., 7.),
 tbv = 3.25607,      fix_tbv = True,    limit_tbv = (2., 4.)
)



#m.tol = 20.0
#m.strategy = 2
m.printMode = 0
m.maxcalls = 4000

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


