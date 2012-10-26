#!/usr/bin/env python

#import pylab
#import matplotlib.pyplot as plt

import shelve, copy, sys, logging

import numpy as np
import scipy.stats

logging.basicConfig(level=logging.INFO)

import Model, Approach, Fitter, Data, utils, plots
from constants import Mp, Mp2

from results import *
from math import sqrt

from dispersion import *
from quadrature import rthtquadrature

#class HFilter(logging.Filter):

#    def __init__(self, treestring):
#        self.tree = treestring.split('.')
#        logging.Filter.__init__(self)

#    def filter(self, rec):
#        head = rec.name.split('.')[:len(self.tree)]
#        if head != self.tree: 
#            return 0
#        #sys.stderr.write(' --- HFilter hit by %s (%s)\n' % (rec.name, rec.levelname))
#        return 1

#_lg = logging.getLogger('p')

#hfil = HFilter('p')
#logging._handlerList[0].addFilter(hfil)

#lg = logging.Logger('A')
#lg.addHandler(logging.StreamHandler())
#fil = HFilter('A')
#lg.handlers[0].addFilter(hfil)
#lg.addFilter(hfil)
#lg.setLevel(logging.DEBUG)  #DEBUG, INFO, WARNING, ERROR, CRITICAL

## [1] Load experimental data and theoretical models

data = utils.loaddata('/home/kkumer/pype/data/ep2epgamma', approach=Approach.BMK)  
data.update(utils.loaddata('/home/kkumer/pype/data/gammastarp2gammap', approach=Approach.BMK))
#data.update(utils.loaddata('/home/kkumer/pype/data/gammastarp2gammap/EIC', approach=Approach.BMK))
#data.update(utils.loaddata('/home/kkumer/pype/data/ep2epgamma/EIC', approach=Approach.BMK))
db = shelve.open('/home/kkumer/pype/theories.db')
#dell = shelve.open('/home/kkumer/pype/dellB.db')

## [2] Choose subset of datapoints for fitting

####  --  Unpolarized target --
#
## H1 and ZEUS
#
DVCSpoints = data[36] + data[37] + data[38] + data[39] + \
  data[40] + data[41] + data[42] + data[43] + data[44] + \
  data[45]
H1ZEUSpoints = DVCSpoints + data[48]
#H1ZEUSindependent = data[45] + data[39] + data[36] + data[46]
H1ZEUSindependentNEW = data[45] + data[39] + data[63] + data[46]
H1ZEUS = H1ZEUSindependentNEW + utils.select(data[47], criteria=['Q2 >= 4.0'])
#
## HERMES 
#
ALUIpoints = utils.select(data[68], criteria=['FTn == -1'])  # HERMES
BCA0points = utils.select(data[67], criteria=['FTn == 0'])  # HERMES
BCA1points = utils.select(data[67], criteria=['FTn == 1'])  # HERMES
ALUIpts = ALUIpoints[:6]
BCApts = BCA0points[:6] + BCA1points[:6]

#
## CLAS
#
BSACLAS_KKpoints = data[25]
BSACLAS_DMpoints = data[8]
CLASpts = utils.select(data[8], criteria=['Q2 >= 2.0'])
#
# Hall A
#
BSDwpoints = utils.select(data[50], criteria=['FTn == -1'])
BSSwpoints = utils.select(data[51], criteria=['FTn>=0', 'FTn <= 1'])
HApts = BSDwpoints[::2] + BSSwpoints[::2]
#
# EIC mock
#
#EICX = data[2001]
#for n in range(2002,2024):
#    EICX = EICX + data[n]
#EICTSA = data[2102]
#for n in range(2103,2110) + range(2111,2118) + range(2119,2125):
#    EICTSA = EICTSA + data[n]
#EICmockkk = data[1002]


####  --  Longitudinally polarized target --
#
TSA1points = utils.select(data[52], criteria=['FTn == -1'])  # HERMES A_UL
TSApoints = TSA1points + data[54]  # HERMES+CLAS  A_UL
BTSApoints = utils.select(data[53], criteria=['FTn==0'])   # HERMES A_LL
LPpoints = TSApoints + BTSApoints  # total longitudinal target
AULpts = TSA1points[:4] + data[54]
ALLpts = BTSApoints[:4]

####  --  Transversally polarized target --
#
AUTIpoints = utils.select(data[66], criteria=['FTn==1'])  # HERMES A_UT_I
AUTDVCSpoints = data[65]  # HERMES A_UT_DVCS
TPpoints = AUTIpoints + AUTDVCSpoints  # total transversal target
AUTIpts = AUTIpoints[:4]

# Global combinations
#
GLOpoints = data[31][12:] + data[8] + data[29]  # DM's GLO set
#ALTGLOpoints = data[5] + data[25] + data[32][18:]  # KK's CLAS BSA
ALTGLO5points = data[5] + data[8] + data[32][18:]   # DM's CLAS BSA
#UNPpoints = ALTGLOpoints + BSSwpoints + BSDwpoints
UNP5points = ALTGLO5points + BSSwpoints + BSDwpoints
#
GLOall = H1ZEUS[::3] + ALUIpts + BCApts + CLASpts + HApts + AULpts + ALLpts + AUTIpts
GLOfull = (H1ZEUS + ALUIpts + BCApts + BSACLAS_DMpoints + BSSwpoints + BSDwpoints
            + LPpoints + TPpoints)
# Excluding LP
GLOnoL = H1ZEUS[::3] + ALUIpts + BCApts + CLASpts + HApts + AUTIpts
# Excluding problematic Hall A BSS:
GLOnoBSS = H1ZEUS[::3] + ALUIpts + BCApts + CLASpts + BSDwpoints[::2] + AULpts + ALLpts + AUTIpts
#               12         6         12        4           6              4+6      4         4
GLOnoBSS2 = H1ZEUS + ALUIpts + BCApts + CLASpts + BSDwpoints + AULpts + ALLpts + AUTIpts
unppts = [ALUIpts, BCApts[6:], CLASpts, BSSwpoints[::-2]]
polpts = [TSA1points[:4], data[54], BTSApoints[:4], AUTIpoints[:4], AUTDVCSpoints[:4]]

# Local 4-bin fits
# Updated data by Morgan and DM
L4_ALUI = utils.select(data[71], criteria=['FTn == -1'])
L4_AC_0 = utils.select(data[70], criteria=['FTn == 0'])
L4_AC_1 = utils.select(data[70], criteria=['FTn == 1'])
# polarized target data
L4_AUL = utils.select(data[52], criteria=['FTn == -1'])
L4_ALL_0 = utils.select(data[53], criteria=['FTn==0'])
L4_ALL_1 = utils.select(data[53], criteria=['FTn==1'])
L4_AUTI_1 = utils.select(data[66], criteria=['FTn==1'])
L4_AUTI_0 = utils.select(data[66], criteria=['FTn==0'])
L4_AUTI_m1 = utils.select(data[66], criteria=['FTn==-1'])
L4_AUTDVCS = data[65]

#bins = zip(L4_ALUI, L4_AC_0, L4_AC_1, L4_AUL, L4_ALL_0, L4_AUTI_1)
bins = zip(L4_ALUI, L4_AC_0, L4_AC_1, L4_AUL, L4_ALL_0, 
        L4_ALL_1, L4_AUTI_1, L4_AUTI_0, L4_AUTI_m1, L4_AUTDVCS)
#bins = zip(L4_ALUI, L4_AC_0, L4_AC_1, L4_AUL, L4_ALL_0, 
#        L4_ALL_1, L4_AUTI_1, L4_AUTI_0, L4_AUTI_m1)
#bins = zip(L4_ALUI, L4_AUL, L4_AUTI_1, L4_AUTI_m1)
#bins = zip(L4_ALUI, L4_AC_0, L4_AC_1, L4_AUL, L4_ALL_0, L4_ALL_1)

## [3] Create a theory


#thAFKM12 = db['AFKM12']
#th = thAFKM12
#Model.ComptonGepard.gepardPool.pop()
#thKM10 = db['KM10']
#Model.ComptonGepard.gepardPool.pop()
#theories = [thAFKM12, thKM10]

m = Model.ModelLocal()
th = Approach.BM10(m)
th.name = 'DMlocal'
th.m.parameters.update({'pImE': -50.02097709090083,
 'pImEt': -94.04223960082038,
 'pImH': 26.90785517414438,
 'pImHt': -4.4912953169971015,
 'pReE': -172.74446204986222,
 'pReEt': 16.30482089309816,
 'pReH': 2.0973238402092345,
 'pReHt': 2.03359772383489})



## [4] Do the fit
#th.model.fix_parameters('ALL')
#th.model.release_parameters(
#   'ALPS', 'M02S', 'SECS', 'THIS', 'ALPG', 'M02G', 'SECG', 'THIG')
#th.model.release_parameters(
#   'EAL0S', 'EALPS', 'EM02S', 'ESECS', 'ETHIS', 'KAPS',
#   'EAL0G', 'EM02G',  'ESECG')
#th.model.release_parameters(
#  'KAPS',  'M02S', 'M02G', 'SECS', 'SECG', 'EAL0S', 'EM02S', 'ESECS', 'EAL0G')
#th.model.release_parameters(
#    'M02S', 'M02G', 'SECS', 'SECG')
#th.model.release_parameters(
#   'rv', 'Mv', 'bv', 'C', 'MC', 'trv', 'tbv')
#th.model.release_parameters('M02S', 'SECS', 'SECG', 'THIS', 'THIG', 
#   'rv', 'bv', 'Mv', 'C', 'MC', 'trv', 'tbv', 'tMv', 'rpi', 'Mpi')
#f = Fitter.FitterMinuit(GLOnoBSS2+BSSwpoints, th)

#th.model.release_parameters('pImH', 'pReH', 'pImE', 'pReE', 'pImHt', 'pReHt')
#nbin = 1
#th.model.release_parameters('pImH', 'pImHt', 'pImE', 'pImEt')
#th.name = th.name + 'bin %s' % nbin
#f = Fitter.FitterMinuit(bins[nbin-1], th)
#f.minuit.tol = 80
#f.minuit.printMode = 2
#f.minuit.maxcalls = 100


#f.fit()
#fl = open('aux.par')
#fl.write(str(th.m.chisq(f.fitpoints)))
#fl.write(str(th.m.parameters))
#fl.write('\n')
#fl.write(str(th.m.covariance))
#fl.close()
#f.fit()


## [5] Some shortcuts ...

def ld(db):
    utils.listdb(db)

GLO12 = H1ZEUS + UNP5points + LPpoints + TPpoints

def pc(th):
    #exps = ['UNP5points', 'ALTGLO5', 'CLAS', 'CLASDM', 'BSDw', 'BSSw', 'TSA1', 'BTSA', 'TPpoints']
    #ptssets = [UNP5points, ALTGLO5points, data[25], data[8], BSDwpoints, BSSwpoints, TSA1points, BTSApoints, TPpoints]
    exps = ['H1ZEUS', 'ALUIpts', 'BCApts', 'CLASpts', 'BSDwpoints', 'BSSwpoints', 'AULpts', 'ALLpts', 'AUTIpts' ]
    ptssets = [H1ZEUS, ALUIpts, BCApts, CLASpts, BSDwpoints, BSSwpoints, AULpts, ALLpts, AUTIpts ]
    for name, pts in zip(exps,ptssets):
        print '%10s: chi/npts = %6.2f/%d' % (name, th.chisq(pts)[0], len(pts))
        #cutpts = utils.select(pts, criteria=['Q2>=1.6'])
        #print '%10s: chi/npts = %6.2f/%d (cut)' % (name, th.chisq(cutpts)[0], len(cutpts))

def pcs(th):
    exps = ['BSAs', 'BCAs']
    ptssets = [Hpoints[:6], Hpoints[18:24]]
    for name, pts in zip(exps,ptssets):
        print '%10s: chi/npts = %6.2f/%d' % (name, th.chisq(pts)[0], len(pts))
        #cutpts = utils.select(pts, criteria=['Q2>=1.6'])
        #print '%10s: chi/npts = %6.2f/%d (cut)' % (name, th.chisq(cutpts)[0], len(cutpts))
    
# fixed target datapoint FIXME: WRONG!
def ptfix(th, Q=1, pol=-1, Ee=160., xB=0.1, Q2=2.2, t=-0.1, phi=3.5, FTn=None):
    ptf = Data.DummyPoint()
    ptf.in1energy = Ee
    ptf.s = 2 * Mp * ptf.in1energy + Mp2
    ptf.in1charge = Q
    ptf.in1polarization = pol
    ptf.xB = xB
    ptf.Q2 = Q2
    ptf.t = t
    ptf.xi = ptf.xB/(2.-ptf.xB)
    ptf.phi = phi
    ptf.frame = 'Trento'
    ptf.units = {'phi': 'radian'}
    utils.fill_kinematics(ptf)
    th.to_conventions(ptf)
    th.prepare(ptf)
    return ptf


# collider datapoint
def ptcol(th, Q=-1, pol=0, Ee=20, Ep=250, xB=0.001998, Q2=4., t=-0.1, 
        phi=np.pi, varphi=-np.pi/2., FTn=None):
#def ptcol(th, Q=-1, pol=0, Ee=20, Ep=250, xB=0.002, Q2=7.3, t=-0.275, 
#        phi=np.pi, varphi=-np.pi/2., FTn=None):
    ptc = Data.DummyPoint()
    ptc.in1energy = Ee
    ptc.in2energy = Ep
    ptc.s = 2 * ptc.in1energy * (ptc.in2energy + sqrt(
                ptc.in2energy**2 - Mp2)) + Mp2
    ptc.in1charge = Q
    ptc.in1polarization = pol
    ptc.in2polarizationvector = 'T'
    ptc.in2polarization = 1 # relevant only for XLP and TSA
    ptc.xB = xB
    ptc.Q2 = Q2
    ptc.t = t
    ptc.xi = ptc.xB/(2.-ptc.xB)
    if phi:
        ptc.phi = phi
        ptc.units = {'phi': 'radian'}
    elif FTn:
        ptc.FTn = FTn
    ptc.varphi = varphi
    ptc.frame = 'Trento'
    utils.fill_kinematics(ptc)
    th.to_conventions(ptc)
    th.prepare(ptc)
    return ptc

#ptc = ptcol(th, Q=1, pol=0, Ee=5, Ep=100, phi=0.1, Q2=4.4, xB=8.2e-3, t=-0.25)
ptI = ptcol(th, Q=1, pol=0, Ee=5, Ep=100, phi=np.pi, Q2=4.4, xB=5.1e-3, t=-0.25)
ptII = ptcol(th, Q=1, pol=0, Ee=20, Ep=250, phi=np.pi, Q2=4.4, xB=5.1e-4, t=-0.25)

def ccals(th, pt):
    cals = ['DVCSunp', 'INTunp', 'INTunpV', 'INTunpA']
    for c in cals:
        print '%8s =  %10.5f + %10.5f * I ' % (c, getattr(th, 'CCAL'+c)(pt), 
                getattr(th, 'CCAL'+c)(pt, im=1))


def _derpt(th, p, pt, f=False, h=0.05):
    """Compute derivative of f w.r.t. model parameter p at point pt.
    
    Simple difference is used (f(p+h/2)-f(p-h/2))/h.
    f is string representing appropriate method of th, or
    observable will be taken as yaxis of pt

    """
    if f:
        fun = th.__getattribute__(f) 
    else:
        fun = th.__getattribute__(pt.yaxis)
    mem = th.m.parameters[p]
    th.m.parameters[p] = mem+h/2.
    up = fun(pt)
    th.m.parameters[p] = mem-h/2.
    down = fun(pt)
    th.m.parameters[p] = mem
    return (up-down)/h


def der(th, pars, pts, f=False,  h=0.05):
    """Compute average derivative at points for each par in pars."""

    for par in pars:
        ders = np.array([_derpt(th, par, pt, f, h) for pt in pts])
        print '%4s  |  %5.2f' % (par, ders.mean())

#pti = data[66][13]
#ptd = data[65][0]

ptH1ZEUS = H1ZEUS[8]
ptHallA = BSSwpoints[::-2][-1]
ptHERMES = TSA1points[1]
ptCLAS = data[54][1]

def CFFatpt(th, cffs=['ImH', 'ReH', 'ImE', 'ReE', 'ImHt', 'ReHt', 'ImEt', 'xReEt']):
    pts = [ptH1ZEUS, ptHallA, ptCLAS, ptHERMES]
    ptss = ['H1ZEUS', 'HallA', 'CLAS', 'HERMES']
    tmpl = '%5s | ' + len(pts)*' % 6.2f '
    print ('%5s | ' + len(pts)*' % 6s ') % tuple([''] + ptss)
    print '------+' + 33*'-'
    for kin in ['xB', 'Q2', 'tm']:
        vals = [getattr(pt, kin) for pt in pts]
        print tmpl % tuple([kin] + vals)
    print '------+' + 33*'-'
    for cff in cffs:
        if cff != 'xReEt':
            vals = []
            for pt in pts:
                if hasattr(th.m, 'g'):
                    th.m.g.newcall = 1
                vals.append(getattr(th.m, cff)(pt))
        else:
            vals = [pt.xi*th.m.ReEt(pt) for pt in pts]
        print tmpl % tuple([cff]+vals)
    print '------+' + 33*'-'



#def rth(m, pt, tht):
#    """Calculate <r(tht)> for model m."""
#    norm = rthtquadrature(lambda r: m.gpdHQbpolpol(pt, r, tht), 0.0, 2.5)
#    aux = rthtquadrature(lambda r: r*m.gpdHQbpolpol(pt, r, tht), 0.0, 2.5)
#    return aux/norm

#def bsq(m, pt):
#    """Calculate <b^2> for model m."""
#    norm = rthtquadrature(lambda b: b*m.gpdHQbIntg(pt, b), 0.0, 2.5)
#    aux = rthtquadrature(lambda b: b**3*m.gpdHQbIntg(pt, b), 0.0, 2.5)
#    return aux/norm
