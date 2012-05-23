#!/usr/bin/env python

"""Plot predictions for HERMES according to two KM10 models."""

# Results sent to Morgan and Caro produced with svn version 193


import sys
import numpy as np

import Model, Approach, Data
import utils 

from constants import Mp, Mp2
from results import *

## Define theories/models

if len(sys.argv) != 2:
    sys.stdout.write('Usage: morgan12 [ KM10a | KM10b ]')
    sys.exit(1)

if sys.argv[1] == 'KM10a':
    mGepard = Model.ComptonGepard(cutq2=0.5)
    mDRsea = Model.ComptonModelDRsea()
    m = Model.HybridDipole(mGepard, mDRsea)
    th = Approach.hotfixedBMK(m)
    th.name = 'KM10a'
    g = th.m.g
    th.m.parameters.update(KM10a)  # ALTGLO
    th.description = 'Model KM10a (without Hall A), arXiv: 1105.0899,1108.1713'
elif sys.argv[1] == 'KM10b':
    mGepard = Model.ComptonGepard(cutq2=0.5)
    mDRPPsea = Model.ComptonModelDRPPsea()
    m = Model.HybridKelly(mGepard, mDRPPsea)
    th = Approach.BM10(m)
    th.name = 'KM10b'
    g = th.m.g
    th.m.parameters.update(KM10b)  # DM email only!
    th.description ='Model KM10b (with Hall A), arXiv: 1105.0899,1108.1713'
else:
    sys.stdout.write('Unknown model: %d\n' % sys.argv[1])
    sys.exit(1)


print "Calculating predictions for HERMES with recoil ..."

kins = [
(0.102, 2.640, -0.130),
(0.084, 2.138, -0.037),
(0.100, 2.514, -0.095),
(0.115, 3.011, -0.201),
(0.131, 3.860, -0.408),
(0.055, 1.469, -0.114),
(0.084, 2.162, -0.111),
(0.122, 3.136, -0.132),
(0.197, 5.063, -0.198),
(0.058, 1.252, -0.098),
(0.080, 1.869, -0.111),
(0.109, 2.836, -0.131),
(0.170, 4.894, -0.188)
]

obs = [('BSA', -1), ('ALUI', -1), ('ALUDVCS', -1), ('ALUI', -2)]


sys.stdout.write('\n')
sys.stdout.write('# '+ th.description + '\n')
sys.stdout.write
sys.stdout.write('#%3s  %5s  %5s   %8s  %8s  %8s  %8s\n' % 
        ('xB', 'Q2', ' t', 'ALU1', 'ALUI1', 'ADVCS1', 'ALUI2') )
sys.stdout.write('# '+ 72*'-' + '\n')
for kin in kins:
    pt = Data.DummyPoint()
    pt.in1energy = 27.6
    pt.s = 2 * Mp * pt.in1energy + Mp2
    pt.in1charge = 1
    pt.in1polarization = 1
    pt.xB, pt.Q2, pt.t = kin
    pt.frame = 'Trento'
    pt.units = {'phi': 'radian'}
    utils.fill_kinematics(pt)
    th.to_conventions(pt)
    th.prepare(pt)
    res = []
    for o in obs:
        pt.FTn = o[1]
        res.append(getattr(th, o[0])(pt))
    sys.stdout.write('%.3f  %.3f  %.3f   % 8.5f % 8.5f % 8.5f % 8.5f\n' % tuple(list(kin)+res))







