#!/usr/bin/env python

"""Plot predictions for HERMES according to two NPB models."""

# Results sent to Morgan Murray produced with svn version 130


import sys
import numpy as np

import Model, Approach, Data
import utils 

from constants import Mp, Mp2
from results import *

## [3] Create theories/models

# Gepard only
mGepard = Model.ComptonGepard(cutq2=0.5)
tGepard = Approach.hotfixedBMK(mGepard)


# DR only
mDRonly = Model.ModelDR()
tDR = Approach.hotfixedBMK(mDRonly)
tDR.name = '(H1/ZEUS)+HERMES+CLAS'
tDR.m.parameters.update(DMepsGLO)

mDRonly1 = Model.ModelDR()
tDR1 = Approach.hotfixedBMK(mDRonly1)
tDR1.name = '(H1/ZEUS)+HERMES+CLAS+HallA'
tDR1.m.parameters.update(DMepsGLO1)


## Hybrid: Gepard+DR (can reuse above Gepard)
#mDRsea = Model.ComptonModelDRsea()
#m = Model.HybridDipole(mGepard, mDRsea)
#t = Approach.hotfixedBMK(m)
#t.name = 'DR + Gepard sea'
#g = t.m.g

mDRPPsea = Model.ComptonModelDRPPsea()
m = Model.HybridDipole(mGepard, mDRPPsea)
th = Approach.BM10(m)
th.name = 'prelim. H1/ZEUS+HERMES+CLAS+HallA'
g = th.m.g
th.m.parameters.update(KKunp5)




models = [th, tDR, tDR1]

print "Calculating predictions for HERMES ..."

kins = [
(0.118, 0.097,  2.505),
(0.019, 0.069, 1.722),
(0.044, 0.088, 2.245),
(0.079, 0.099, 2.493),
(0.143, 0.109, 2.759),
(0.261, 0.119, 3.232),
(0.463, 0.122, 3.727),
(0.099, 0.049, 1.343),
(0.093, 0.070, 1.794),
(0.106, 0.089, 2.303),
(0.122, 0.114, 2.936),
(0.160, 0.157, 4.062),
(0.233, 0.244, 6.126),
(0.078, 0.055, 1.199),
(0.092, 0.069, 1.590),
(0.106, 0.085, 2.079),
(0.127, 0.105, 2.768),
(0.152, 0.134, 3.766),
(0.220, 0.199, 5.789),
]

obs = [('BCA', 0), ('BCA', 1), ('BCA', 2), ('ALUI', -1), ('ALUDVCS', -1), ('ALUI', -2)]

models = [tDR, tDR1]
descriptions = ['Model KM09a (without Hall A)  NPB841 (2010) 1-58, arXiv:0904.0458',
                'Model KM09b (with Hall A)  NPB841 (2010) 1-58, arXiv:0904.0458']

for m, d in zip(models, descriptions):
    sys.stdout.write('\n')
    sys.stdout.write('# '+ d + '\n')
    sys.stdout.write
    sys.stdout.write('#%3s  %5s  %5s   %8s %8s %8s %8s %8s %8s\n' % 
            ('-t', 'xB', 'Q2', 'BCA0',  'BCA1', 'BCA2', 'ALUI1', 'ADVCS1', 'ALUI2') )
    sys.stdout.write('# '+ 72*'-' + '\n')
    for kin in kins:
        pt = Data.DummyPoint()
        pt.in1energy = 27.6
        pt.s = 2 * Mp * pt.in1energy + Mp2
        pt.in1charge = 1
        pt.in1polarization = 1
        pt.tm, pt.xB, pt.Q2 = kin
        pt.frame = 'Trento'
        pt.units = {'phi': 'radian'}
        utils.fill_kinematics(pt)
        m.to_conventions(pt)
        m.prepare(pt)
        res = []
        for o in obs:
            pt.FTn = o[1]
            res.append(getattr(m, o[0])(pt))
        sys.stdout.write('%.3f  %.3f  %.3f   % 8.5f % 8.5f % 8.5f % 8.5f % 8.5f % 8.5f\n' % tuple(list(kin)+res))







