
import logging, logzero, os
_lg = logzero.logger
logzero.loglevel(logging.INFO)
basename = os.path.splitext(os.path.basename(__file__))[0]
logfilename="/home/kkumer/tmp/{}.log".format(basename)
logzero.logfile(logfilename, loglevel=logging.INFO)

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


pts = ALUIpts + HA15pts + C_BSDwpts[::4] + C_BSSw0pts[::4] + C_BSSw1pts[::4]  


numn = 10

numpy.random.seed(68)

m = Model.ModelNNKelly(zeropointpower=True,
        output_layer=['ImH', 'ReH', 'ImE', 'ReE', 'ImHt', 'ReHt', 'ImEt', 'ReEt'] )
_lg.info('New model created: {}'.format(m.__class__))

th = Approach.BM10tw2(m)
th.name = "nfNN"
f = Fitter.FitterBrain(pts, th, nnets=numn, nbatch=30, minprob=0.000001)
f.verbose = 1
_lg.info('Start fit to {} data points'.format(len(f.fitpoints)))
f.fitgood()

print('{} NNs {} {:.1f}/{} p={:.3g}'.format(numn, th.name, *th.chisq(f.fitpoints)))

_lg.info('Done. Emailing log file.')
utils.mailfile('kkumer@calculon.phy.hr', 'kkumer@calculon.phy.hr',
        '{} fit is done'.format(basename), logfilename)
