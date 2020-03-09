
import logging, logzero, os
_lg = logzero.logger
logzero.loglevel(logging.INFO)
basename = os.path.splitext(os.path.basename(__file__))[0]
logzero.logfile("/home/kkumer/tmp/{}.log".format(basename), loglevel=logging.INFO)

import sys, shelve
GEPARD_DIR = '/home/kkumer/gepard'
sys.path.append(GEPARD_DIR+'/pype')
import Model, Approach, Fitter, Data, utils, plots
import numpy
from results import *
from abbrevs import *

HEpts = ( L4_ALUI+ L4_AC_0+ L4_AC_1+ L4_ALL_1+ L4_AUTI_1+ L4_AUTI_0+ L4_AUTI_m1+
        L4_AUTDVCS+ L4_ALTI_1+ L4_ALTI_0+ L4_ALTI_m1+ L4_ALTBHDVCS_0)  # HERMES as in KMM12

CLAS14pts = CLAS14TSApts + CLAS14BTSApts + CLAS14BSApts
CLAS15pts = C_BSDwpts + C_BSSw0pts + C_BSSw1pts
HA15pts = H_BSDwpts + H_BSSw0pts + H_BSSw1pts
HA17pts = H17_BSDwpts + H17_BSSw0pts + H17_BSSw1pts

#pts = BCApts + ALUIpts + HA15pts + C_BSDwpts[::4] + C_BSSw0pts[::4] + C_BSSw1pts[::4]  
# increase 89 --> 128 and introduce target polarization data
#pts = C_BSDwpts[::4]+ C_BSSw0pts[::4] + C_BSSw1pts[::4] + HA15pts + BCApts + ALUIpts + LPpoints[::2] + TPpoints

# 1. pts = GLO15new + data[136] - fail

pts = GLO15new

numn = 10

numpy.random.seed(68)

mNN1 = Model.ModelNN(
        output_layer=['ImH', 'ReH', 'ImE', 'ReE', 'ImHt', 'ReHt', 'ImEt', 'ReEt']
        )

tNN1 = Approach.BM10tw2(mNN1)
tNN1.name = "nfNN-eight-J15"
fNN1 = Fitter.FitterBrain(pts, tNN1, nnets=numn, nbatch=30, minprob=0.000001)
fNN1.verbose = 1
print('\n--- NN fit to {} data points ---\n'.format(len(pts)))
fNN1.fitgood()
#dbn = shelve.open(GEPARD_DIR+'/pype/nndr.db')
#tNN1.save(dbn)
#dbn.close()

print('{} NNs {} {:.1f}/{} p={:.3g}'.format(numn, tNN1.name, *tNN1.chisq(pts)))

