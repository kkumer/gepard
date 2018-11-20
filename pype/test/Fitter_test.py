
import copy
from nose.tools import *
import numpy as np


import logging 
mpl_logger = logging.getLogger('matplotlib') 
mpl_logger.setLevel(logging.WARNING)

import utils, Model, Approach, Fitter
from results import DMGLO1, DMepsGLO1  #use some testpars here?

data = utils.loaddata('data/ep2epgamma', approach=Approach.hotfixedBMK)  

m = Model.ModelDR()
m.parameters.update(DMGLO1)
m.parameters['limit_Mv'] = (0.9, 1.1)  # for compatibility with old g.
m.release_parameters('bS', 'Mv')
t = Approach.hotfixedBMK(m, optimization = False)

# Optimized model
mopt = Model.ModelDR(optimization = False)
mopt.parameters.update(DMepsGLO1)
mopt.parameters['Mv'] = 1.5
mopt.parameters['bv'] = 0.4
mopt.parameters['rv'] = 0.798
mopt.parameters['C'] = 5.36
mopt.parameters['MC'] = 2.
mopt.parameters['tMv'] = 2.
mopt.parameters['trv'] = 8.
mopt.parameters['tbv'] = 3.78
mopt.ndparameters = np.array([mopt.parameters[name] for name in mopt.parameter_names]+[0.,0.])
mopt.release_parameters('bS', 'Mv')
topt = Approach.BM10(mopt, optimization = False)


# testing data set
testpoints = [data[31][12]] + [data[8][1]] + [data[29][2]] + [data[30][3]]
# testing data set for fits
fitpoints = data[31][12:14] + data[8][1:3] + data[30][2:4]


def test_fit():
    """Testing set of fitting observables."""
    chisq = 0.
    for pt in testpoints:
        chisq = chisq + (
                (getattr(t, pt.yaxis)(pt) - pt.val)**2 / pt.err**2 )
    assert_almost_equal(chisq, 8.8479922732922276)

def test_fit2():
    """Testing fitting by FitterMinuit."""
    f = Fitter.FitterMinuit(fitpoints, t)
    f.fit()
    chisq = t.chisq(fitpoints)[0]
    #assert_almost_equal(chisq, 8.4891170857950087, 4)
    # After pyminuit --> iminuit
    assert_almost_equal(chisq,  8.139870272995228, 4)

test_fit2.long = 1
test_fit2.batch = 1

def test_fit2opt():
    """Testing optimized (not really) fitting by FitterMinuit."""
    fopt = Fitter.FitterMinuit(fitpoints, topt)
    fopt.fit()
    chisq = topt.chisq(fitpoints)[0]
    assert_almost_equal(chisq, 7.3648989640243645, 4)

test_fit2opt.long = 1
test_fit2opt.batch = 1
test_fit2opt.optimization = 1

