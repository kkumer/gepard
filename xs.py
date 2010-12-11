#!/usr/bin/env python


import sys
from math import sqrt

import Data, Model, Approach, utils

from constants import Mp, Mp2

from results import *


def theory(id):
    """ Return requested theory. """

    if id == 0:
        # for debugging, always returns 42
        class dummy(object):
            def to_conventions(self, pt):
                pass
            def prepare(self, pt):
                pass
            def Xunp(self, pt):
                return 42
        dummyo = dummy()
        return dummyo
    elif id == 1:
        # NPB fit withouth Hall A
        mDRonly = Model.ModelDR()
        tDR = Approach.hotfixedBMK(mDRonly)
        tDR.name = '(H1/ZEUS)+HERMES+CLAS'
        tDR.m.parameters.update(DMepsGLO)
        return tDR
    elif id == 2:
        # NPB fit with Hall A
        mDRonly1 = Model.ModelDR()
        tDR1 = Approach.hotfixedBMK(mDRonly1)
        tDR1.name = '(H1/ZEUS)+HERMES+CLAS+HallA'
        tDR1.m.parameters.update(DMepsGLO1)
        return tDR1
    elif id == 3:
        # Hybrid preliminary fit shown in Trento
        mGepard = Model.ComptonGepard(cutq2=0.5)
        mDRPPsea = Model.ComptonModelDRPPsea()
        m = Model.Hybrid(mGepard, mDRPPsea)
        t = Approach.BM10(m)
        t.name = 'prelim. H1/ZEUS+HERMES+CLAS+HallA'
        g = t.m.g
        t.m.parameters.update(KKunp5)
        return t
    else:
        sys.stdout.write('Unknown model: %d' % id)
        sys.exit(1)
    
def tmin(xB, Q2):
    """BMK Eq. (31)"""

    eps2 = (4. * xB**2 * Mp2) / Q2
    return (-Q2 * ( 2. * (1.-xB)*(1. - sqrt(1.+eps2)) + eps2 ) / (
            4. * xB * (1.-xB) + eps2 ))

def XSunp(id, Q, lam, Ee, Ep, xB, Q2, t, phi):
    """ Return cross section of scattering on unpolarized proton

    Here phi is assumed to be in Trento conventions.
    """

#    sys.stderr.write('called grid(%i, %i, %i, %.1f, %.1f, 
#       %.2f, %.2f, %.2f, %.2f)\n' % (
#       id, Q, lam, Ee, Ep, xB, Q2, t, phi))
    pt0 = Data.DummyPoint()
    pt0.in1energy = Ee
    if Ep > Mp:
        pt0.in2energy = Ep
        pt0.s = 2 * pt0.in1energy * (pt0.in2energy + sqrt(
                    pt0.in2energy**2 - Mp2)) + Mp2
    else:
        pt0.s = 2 * Mp * pt0.in1energy + Mp2
    if xB<0.000001 or xB>0.4:
        sys.stdout.write('xB = %.4g is not in reliable range 1e-6 < xB < 0.4\n' % xB)
        sys.stdout.flush()
        sys.exit(1)
    elif (Q2 >= (pt0.s-Mp2)*xB) or (t >= tmin(xB,Q2)):
        return 0
    else:
        pt0.in1charge = Q
        pt0.in1polarization = lam
        pt0.xB = xB
        pt0.Q2 = Q2
        pt0.t = t
        pt0.xi = pt0.xB/(2.-pt0.xB)
        pt0.phi = phi
        pt0.units = {'phi': 'radian'}
        pt0.frame = 'Trento'
        utils.fill_kinematics(pt0)
        t = theory(id)
        t.to_conventions(pt0)
        t.prepare(pt0)
        return t.Xunp(pt0)

if __name__ == '__main__':
    usage = """

  xs.exe  ModelID  Charge  Polarization  Ee  Ep  xB  Q2  t  phi

returns cross section (in nb) for scattering of lepton of energy Ee 
on unpolarized proton target of energy Ep. Charge=-1 is for electron. 

ModelID is one of  
   0 debug, always returns 42, 
   1 arXiv:0904.0458 fit without Hall A,
   2 arXiv:0904.0458 fit with Hall A, 
   3 preliminary hybrid fit with LO sea evolution, from Trento presentation.

xB Q2 t phi  -- usual kinematics (phi in Trento convention)

For convenience, if last argument (phi=n) is larger than 2pi, you get grid 
of n equidistant points with phi=0..2pi.

Example:
    ./xs.exe 1 -1 1 27.6 0.938  0.111 3. -0.3 20

"""

    try:
        args = [int(k) for k in 
                sys.argv[1:4]] + [float(k) for k in sys.argv[4:]]
    except:
        sys.stdout.write(usage)
        sys.exit(1)
    if len(args) != 9:
        sys.stdout.write(usage)
        sys.exit(1)
    if args[-1] <= 2*3.141592:
        res = apply(XSunp, args)
        sys.stdout.write("%.14f\n" % res)
    else:
        # we want phi-grid
        res = []
        npts = int(args[-1]) # number of points
        for k in range(npts):
            phi = 2*3.141*k / npts
            args[-1] = phi
            sys.stdout.write("%.3f  %.14f\n" % (phi, apply(XSunp, args)))
    sys.stdout.flush()
