""" 
Outputing grids of cross-sections for a specific model.

"""

import sys, os, math
import numpy as np

import Data, utils
from constants import Mp, Mp2


# load experimental data
#data = utils.loaddata('data/ep2epgamma', approach=Approach.hotfixedBMK) 
#data.update(utils.loaddata('data/gammastarp2gammap', approach=Approach.hotfixedBMK)) 

def tmin(xB, Q2):
    """BMK Eq. (31)"""

    eps2 = (4. * xB**2 * Mp2) / Q2
    return min(-Q2 * ( 2. * (1.-xB)*(1. - math.sqrt(1.+eps2)) + eps2 ) / (
            4. * xB * (1.-xB) + eps2 ), -0.15)


def grid(th, in1energy, in2energy=Mp, filename='out.grid'):
    """Makes grid in ...
    
    """
    f = open(filename, 'w')
    if in2energy > Mp:
        # collider
        s = 2 * in1energy * (in2energy + math.sqrt(
            in2energy**2 - Mp2)) + Mp2
    else: 
        # fixed target
        s = 2 * Mp * in1energy + Mp2
    #xmin = -3  # this will be interpreted as 10^-3 !
    xmax = 0.4
    Q2min = 1.2
    Q2max = 8.0
    tmn = -0.0
    tmax = -0.6
    #xs = np.concatenate((np.logspace(-3, -1, 10), np.linspace(0.1, xmax, 10)))
    xs = np.linspace(1.5/s, xmax, 10)
    qs = np.linspace(Q2min, Q2max, 10)
    ts = np.linspace(tmax, tmn, 10)
    phis = np.linspace(0, 2*np.pi, 10)
    fs = "%5.3f  %6.3f  %6.3f  %5.3f    %8.3e  %8.3e  %8.3e  % 8.3e  % 8.3e\n"
    for xB in xs:
        for Q2 in qs:
            for t in ts:
                for phi in phis:
                    if (Q2 >= (s-Mp2)*xB) or (t >= tmin(xB,Q2)):
                        f.write(fs % (xB, Q2, t, phi, 0, 0, 0, 0, 0))
                    else:
                        pt = Data.DummyPoint()
                        pt.s = s
                        pt.in1charge = -1
                        pt.in1polarization = 1
                        pt.xB = xB
                        pt.Q2 = Q2
                        pt.t = t
                        pt.xi = xB/(2.-xB)
                        pt.phi = phi
                        pt.frame = 'BMK'
                        utils.fill_kinematics(pt)
                        th.prepare(pt)
                        #
                        pf = th.PreFacSigma(pt)
                        TBH = pf * th.TBH2unp(pt) 
                        TDVCSp = pf * th.TDVCS2unp(pt)
                        TINTp = pf * th.TINTunp(pt) 
                        pt.in1polarization = -1
                        TDVCSm = pf * th.TDVCS2unp(pt)
                        TINTm = pf * th.TINTunp(pt) 
                        D0 = (TDVCSp + TDVCSm)/2
                        D1 = (TDVCSp - TDVCSm)/2
                        I0 = (TINTp + TINTm)/2
                        I1 = (TINTp - TINTm)/2
                        #Xtest = TBH + (D0 - D1) + (I0 - I1)
                        #X = th.Xunp(pt)
                        f.write(fs % (xB, Q2, t, phi, TBH, D0, D1, I0, I1))
    f.close()

