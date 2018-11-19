#!/usr/bin/env python

"""Plot predictions for COMPASS according to two NPB models."""

# Results sent to Nicole d'Hose produced with SVN revision 88.

import shelve

import numpy as np

import Model, Approach, Fitter
import utils 
import plots

from results import *

# [1] Load experimental data and theoretical models

data = utils.loaddata('data/ep2epgamma', approach=Approach.hotfixedBMK)
data.update(utils.loaddata('data/gammastarp2gammap', approach=Approach.hotfixedBMK))
db = shelve.open('theories.db')



## [2] Choose subset of datapoints for fitting

GLOpoints = data[31][12:] + data[8] + data[29]  # DM's GLO set
GLO1points = data[31][12:] + data[8] + data[29] + data[30]  # DM's GLO1 set


## [3] Create a theory


# DR only
mDRonly = Model.ModelDR()
tDR = Approach.hotfixedBMK(mDRonly)
tDR.name = 'HERMES+CLAS'

mDRonly1 = Model.ModelDR()
tDR1 = Approach.hotfixedBMK(mDRonly1)
tDR1.name = 'HERMES+CLAS+HALLA'


tDR.m.parameters.update(DMepsGLO)
tDR1.m.parameters.update(DMepsGLO1)


## [4] Check the fit

#print "fit to HERMES+CLAS ... "
#tDR.print_chisq(GLOpoints)
#print "fit to HERMES+CLAS+HALLA ..."
#tDR1.print_chisq(GLO1points)

print("Plotting predictions for COMPASS ...")

kins = [
(1.4, 0.007),
(1.4, 0.014),
(1.4, 0.024),
(1.4, 0.038),
(1.4, 0.077),
(2.4, 0.015),
(2.4, 0.025),
(2.4, 0.042),
(2.4, 0.080),
(5 , 0.026),
(5 , 0.047),
(5 , 0.085)]


for kin in kins:
    plots.COMPASS([tDR, tDR1], Q2=kin[0], xB=kin[1], path='.', fmt='eps')
