#!/usr/bin/env python2

#from IPython.Debugger import Tracer; debug_here = Tracer()

import sys
from math import sqrt, pi

import Data, Model, Approach, utils

from constants import Mp, Mp2

from results import *


def theory(id):
    """ Return requested theory. """

    if id == 0:
        # for debugging, always returns 42
        class dummymodel(object):
            def is_within_model_kinematics(self, pt):
                return True
        class dummy(object):
            def __init__(self):
                self.model = dummymodel()
                self.m = self.model
            def to_conventions(self, pt):
                pass
            def prepare(self, pt):
                pass
            @staticmethod
            def is_within_phase_space(pt):
                return True
            def Xunp(self, pt):
                return 42
        dummyo = dummy()
        return dummyo
    elif id == 1:
        # NPB fit withouth Hall A
        mDRonly = Model.ModelDR()
        tDR = Approach.hotfixedBMK(mDRonly)
        tDR.name = 'KM09a'
        tDR.m.parameters.update(DMepsGLO)
        return tDR
    elif id == 2:
        # NPB fit with Hall A
        mDRonly1 = Model.ModelDR()
        tDR1 = Approach.hotfixedBMK(mDRonly1)
        tDR1.name = 'KM09b'
        tDR1.m.parameters.update(DMepsGLO1)
        return tDR1
    elif id == 3:
        # Hybrid preliminary fit shown in Trento
        # P(132.80, 160) = 0.943  (for UNP5 and H1ZEUS)
        mGepard = Model.ComptonGepard(cutq2=0.5)
        mDRPPsea = Model.ComptonModelDRPPsea()
        m = Model.HybridDipole(mGepard, mDRPPsea)
        th = Approach.BM10(m)
        th.name = 'KM10'
        g = th.m.g
        th.m.parameters.update(KM10)  # KKunp5
        return th
    elif id == 4:
        # Hybrid fit + 3rd PW to DVCSpoints + data[48] + ALTGLOpoints
        # P(133.06, 155) = 0.8983
        mGepard = Model.ComptonGepard(cutq2=0.5)
        mDRsea = Model.ComptonModelDRsea()
        m = Model.HybridDipole(mGepard, mDRsea)
        th = Approach.hotfixedBMK(m)
        th.name = 'KM10a'
        g = th.m.g
        th.m.parameters.update(KM10a)  # ALTGLO
        return th
    elif id == 5:
        # Hybrid fit + 3rd PW to DVCSpoints and GLO1points
        # P(42.51/40) = 0.3633
        mGepard = Model.ComptonGepard(cutq2=0.5)
        mDRPPsea = Model.ComptonModelDRPPsea()
        m = Model.HybridKelly(mGepard, mDRPPsea)
        th = Approach.BM10(m)
        th.name = 'KM10b'
        g = th.m.g
        th.m.parameters.update(KM10b)  # DM email only!
        return th
    elif id == 6:
        # Hybrid fit 
        # P(chi-square, d.o.f) = P(124.12, 80) = 0.0012
# 
# GLOnoBSS2 = H1ZEUS + ALUIpts + BCApts + CLASpts + BSDwpoints + AULpts + ALLpts + AUTIpts
# + BSSwpoints
# 
#    Mv =    0.951 +- 0.282  (p = 0.00114)          H1ZEUS: chi/npts =  28.57/35
#    rv =    1.121 +- 0.099  (p = 0)               ALUIpts: chi/npts =  12.47/6
#    bv =    0.400 +- 0.000  (p = 0)                BCApts: chi/npts =   9.35/12
#     C =    1.003 +- 0.565  (p = 0.0794)          CLASpts: chi/npts =   5.14/4
#    MC =    2.080 +- 3.754  (p = 0.581)        BSDwpoints: chi/npts =  12.74/12
#   tMv =    3.523 +- 13.175  (p = 0.79)        BSSwpoints: chi/npts =  18.80/8
#   trv =    1.302 +- 0.206  (p = 1.38e-08)         AULpts: chi/npts =  19.66/10
#   tbv =    0.400 +- 0.001  (p = 0)                ALLpts: chi/npts =  13.77/4
#   rpi =    3.837 +- 0.141  (p = 0)               AUTIpts: chi/npts =   3.61/4
#   Mpi =    4.000 +- 0.036  (p = 0)
#  M02S =    0.462 +- 0.032  (p = 0)
#  SECS =    0.313 +- 0.039  (p = 8.67e-12)
#  THIS =   -0.138 +- 0.012  (p = 0)
#  SECG =   -2.771 +- 0.228  (p = 0)
#  THIG =    0.945 +- 0.107  (p = 2.03e-13)
        mGepard = Model.ComptonGepard(cutq2=0.5)
        mDRPPsea = Model.ComptonModelDRPPsea()
        m = Model.HybridKelly(mGepard, mDRPPsea)
        th = Approach.BM10(m)
        th.name = 'KMM12'
        g = th.m.g
        th.m.parameters.update(KMM12)
        return th
    else:
        sys.stdout.write('Unknown model: %d\n' % id)
        sys.exit(1)
    
def tmin(xB, Q2):
    """BMK Eq. (31)"""

    eps2 = (4. * xB**2 * Mp2) / Q2
    return (-Q2 * ( 2. * (1.-xB)*(1. - sqrt(1.+eps2)) + eps2 ) / (
            4. * xB * (1.-xB) + eps2 ))

def XSall(id, Q, lam, Ee, Ep, xB, Q2, t, phi):
    """ Return cross section of scattering on unpolarized proton

    Here phi is assumed to be in Trento conventions, and results
    are also in Trento convention.
    """

#    sys.stderr.write('called grid(%i, %i, %i, %.1f, %.1f, 
#       %.2f, %.2f, %.2f, %.2f)\n' % (
#       id, Q, lam, Ee, Ep, xB, Q2, t, phi))
    pt0 = Data.DummyPoint()
    pt0.in1energy = Ee
    if Ep > Mp:
        # collider
        pt0.in2energy = Ep
        pt0.s = 2 * pt0.in1energy * (pt0.in2energy + sqrt(
                    pt0.in2energy**2 - Mp2)) + Mp2
    else:
        # fixed target
        pt0.s = 2 * Mp * pt0.in1energy + Mp2
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
    th = theory(id)
    th.to_conventions(pt0)
    th.prepare(pt0)
    if not th.is_within_phase_space(pt0):
        return 0
    elif not th.m.is_within_model_kinematics(pt0):
        sys.stdout.write('Error: (xB, Q2, t) = (%.4g, %.4g, %.4g) is outside of kinematic bounds of model %i\n' % (xB, Q2, t, id))
        sys.stdout.flush()
        sys.exit(1)
    else:
        xstot = [th.Xunp(pt0), 0, 0, 0]
        if id > 5:
            pt0.in2polarizationvector = 'L'
            pt0.in2polarization = 1
            xstot[3] = - th.XLP(pt0)  # minus to go to Trento
            pt0.in2polarizationvector = 'T'
            pt0.varphi = 0
            xstot[1] =  - th.XTP(pt0)
            pt0.varphi = pi/2
            xstot[2] = - th.XTP(pt0)
        return xstot

if __name__ == '__main__':
    usage = """

  xs.exe  ModelID  Charge  Polarization  Ee  Ep  xB  Q2  t  phi

returns cross section (in nb) for scattering of lepton of energy Ee on proton 
of energy Ep. xB, Q2 and t is usual kinematics. Charge=-1 is for electron. 
Polarization=+1 is for lepton polarization along the beam.  Output is:

    phi xs_unp  xs_TPcos  xs_TPsin  xs_LP

where total cross section is

 xs = xs_unp + sin(theta_S) cos(phi-phi_S) xs_TPcos
             + sin(theta_S) sin(phi-phi_S) xs_TPsin
             + cos(theta_S) xs_LP

and theta_S and phi_S are proton polarization polar and
azimuthal angles, while phi is angle between lepton
and reaction planes. All in Trento conventions.

ModelID is one of  
   0 debug, always returns 42, 
   1 KM09a - arXiv:0904.0458 fit without Hall A,
   2 KM09b - arXiv:0904.0458 fit with Hall A, 
   3 KM10  - preliminary hybrid fit with LO sea evolution, from Trento presentation,
   4 KM10a - preliminary hybrid fit with LO sea evolution, without Hall A data
   5 KM10b - preliminary hybrid fit with LO sea evolution, with Hall A data
   6 KMM12 - hybrid global fit to unpolarized and polarized DVCS data
where models 1-5 are for unpolarized target only.
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
        lst = apply(XSall, args)
        lst.insert(0, args[-1])
        if args[0] > 5:
            sys.stdout.write("%.3f  %.14f  % .14f  % .14f  % .14f\n" % tuple(lst))
        else:
            sys.stdout.write("%.3f  %.14f\n" % tuple(lst[:2]))
    else:
        # we want phi-grid
        res = []
        npts = int(args[-1]) # number of points
        for k in range(npts):
            phi = 2*pi*k / (npts-1)
            args[-1] = phi
            lst = apply(XSall, args)
            lst.insert(0, phi)
            if args[0] > 5:
                sys.stdout.write("%.3f  %.14f  % .14f  % .14f  % .14f\n" % tuple(lst))
            else:
                sys.stdout.write("%.3f  %.14f\n" % tuple(lst[:2]))
    sys.stdout.flush()