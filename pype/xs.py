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
        mGepard = Model.ComptonGepard(cutq2=0.5)
        mDRPPsea = Model.ComptonModelDRPPsea()
        m = Model.HybridDipole(mGepard, mDRPPsea)
        th = Approach.BM10(m)
        th.name = 'KMM12'
        g = th.m.g
        th.m.parameters.update(KMM12)
        return th
    elif id == 7:
        # Hybrid fit  KM15   chisq/dof = 240./275
        mGepard = Model.ComptonGepard(cutq2=0.5)
        mDRPPsea = Model.ComptonModelDRPPsea()
        m = Model.HybridKelly(mGepard, mDRPPsea)
        th = Approach.BM10tw2(m)
        th.name = 'KM15'
        g = th.m.g
        th.m.parameters.update(KM15)
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
    pt0.y1name = 'BSS'
    pt0.units = {'phi': 'radian', 'BSS': 'nb/GeV^4'}
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
        #print "BH2  = {}".format(th.TBH2unp(pt0))
        #print "PreFacBH  = {}".format(th.PreFacBH(pt0))
        #print "cBH0 = {}".format(th.cBH0unp(pt0))
        #print "cBH1 = {}".format(th.cBH1unp(pt0))
        #print "cBH2 = {}".format(th.cBH2unp(pt0))
        #print "phi = {}".format(pt0.phi)
        #print "F1 = {}".format(th.m.F1(pt0.t))
        #print "F2 = {}".format(th.m.F2(pt0.t))
        #print "t = {}".format(pt0.t)
        #print th.m.F1
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
and reaction planes. All in radians and Trento conventions.

ModelID is one of  
   0 debug, always returns 42, 
   1 KM09a - arXiv:0904.0458 fit without Hall A data,
   2 KM09b - arXiv:0904.0458 fit with Hall A harmonics ratio,
   3 KM10  - arXiv:1105.0899 fit with Hall A harmonics
   4 KM10a - arXiv:1105.0899 fit without Hall A data
   5 KM10b - arXiv:1105.0899 fit with Hall A harmonics ratio
   6 KMM12 - arXiv:1301.1230 fit with Hall A harmonics and polarized target
   7 KM15  - arXiv:1512.09014 fit now includes 2015 CLAS and Hall A data
where models 1-5 are for unpolarized target only.
For convenience, if last argument (phi=n) is larger than 2pi, you get grid 
of n equidistant points with phi=0..2pi.

Example:
    ./xs.exe  7 -1  1  5.75  0  0.36  2.3  -0.17  0.131

(Output:)
0.131  0.07584357734528  -0.03809893007524   0.00826740897951  -0.03819278821799

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
