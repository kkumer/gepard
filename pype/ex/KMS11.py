#!/usr/bin/env python


import shelve

import plots, utils, Approach

print """
This will produce figures 3, 4 and 5 for KMS11 paper.
For paper, this script was run with pype ver. 149.
"""

## [1] Load experimental data and theoretical models

data = utils.loaddata('/home/kkumer/pype/data/ep2epgamma', approach=Approach.hotfixedBMK)  
data.update(utils.loaddata('/home/kkumer/pype/data/gammastarp2gammap', approach=Approach.hotfixedBMK))
db = shelve.open('/home/kkumer/pype/theories.db')


plots.HERMES09(lines=[db['KM09a'], db['KM09b']], bands=[db['KMS11-NN']])

plots.COMPASSt(lines=[db['KM09a'], db['KM09b']], bands=[db['KMS11-NN']])

plots.CFF2(lines=[db['KM09a'], db['KM09b']], bands=[db['KMS11-NN'], db['KMS11-DR']])

raw_input('Press ENTER to finish:')
