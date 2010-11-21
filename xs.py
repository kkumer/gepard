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
            def prepare(self, pt):
                pass
            def BSA(self, pt):
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
        sys.stderr.write('Unknown model: %d' % id)
        sys.exit(1)
    

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

