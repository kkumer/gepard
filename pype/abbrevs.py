# Abbreviations for various collections of data
import os, utils, Approach

# load experimental data
DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')
data = utils.loaddata(os.path.join(DATA_DIR,'ep2epgamma'), 
       approach=Approach.BMK) 
data.update(utils.loaddata(os.path.join(DATA_DIR, 'gammastarp2gammap'),
       approach=Approach.BMK)) 
data.update(utils.loaddata(os.path.join(DATA_DIR, 'gammastarp2Mp'),
       approach=Approach.BMK)) 
data.update(utils.loaddata(os.path.join(DATA_DIR, 'DIS'),
       approach=Approach.BMK)) 
#data.update(utils.loaddata('/home/kkumer/pype/data/gammastarp2gammap/EIC', approach=Approach.BMK))
#data.update(utils.loaddata('/home/kkumer/pype/data/ep2epgamma/EIC', approach=Approach.BMK))


## [2] Choose subset of datapoints for fitting

DISpoints = data[201] + data[202] + data[203] + data[204] + \
            data[205] + data[206] + data[207] + data[208] + \
            data[209] + data[210] + data[211] + data[212]

H1all = data[36] + data[37] + data[38] + data[39] + \
  data[40] + data[41] + data[42] + data[43] + data[44] + \
  data[60] + data[61] + data[62] + data[63] + data[64] + \
  data[75] + data[76] + data[77] + data[78] + data[79] + data[80]
ZEUSall = data[45] + data[46] + data[47] + data[48] + data[49]
HERMESall = data[29] + data[31] + data[32] + data[52] + data[53] + \
            data[65] + data[66] + \
            data[67] + data[68] + data[69] + data[70] + data[71] + \
            data[72] + data[73] + data[74]
CLASall = data[7] + data[54]
HallAall = data[33] + data[34] + data[55] + data[56]
####  --  Unpolarized target --
#
## H1 and ZEUS
#
DVCSpoints = data[36] + data[37] + data[38] + data[39] + \
  data[40] + data[41] + data[42] + data[43] + data[44] + \
  data[45]
H1ZEUSpoints = DVCSpoints + data[48]
H1ZEUSindependent = data[45] + data[39] + data[36] + data[46]
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
CLASTSApts = data[54][:3]
# [-3:] means cut Q2>1.51 GeV^2
CLAS08pts = utils.select(data[81], criteria=['FTn == -1'])[-3:]
CLASptsOLD = utils.select(data[8], criteria=['Q2 >= 2.0'])
CLASpts = utils.select(data[8], criteria=['Q2 >= 1.5']) + CLAS08pts
CLASKKpts = utils.select(data[25], criteria=['Q2 >= 1.5']) + CLAS08pts
CLAS14BSApts = data[85]
CLAS14TSApts = utils.select(data[86], criteria=['FTn == -1']) 
CLAS14BTSApts = data[87]
C_BSDwpts = data[101]
C_BSSw0pts = data[102][:48]
C_BSSw1pts = data[102][48:]
C_BSD = utils.select(data[99], criteria=['val != 0'])
C_BSS = utils.select(data[100], criteria=['val != 0'])
#
# Hall A
#
BSDwpoints = utils.select(data[50], criteria=['FTn == -1'])
BSSwpoints = utils.select(data[51], criteria=['FTn>=0', 'FTn <= 1'])
HApts = BSDwpoints[::2] + BSSwpoints[::2]
H_BSDwpts = data[117]
H_BSSw0pts = data[116][:10]
H_BSSw1pts = data[116][10:]
H_BSDpts = data[121]
H_BSS0pts = data[120][:10]
H_BSS1pts = data[120][10:]
H_BSD = data[109]+data[110]+data[111]
H_BSS = data[107]+data[108]
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
#
# H1 DVMP points
#
H109XL = utils.select(data[76], criteria=['Q2 >= 4.0'])
H109tdep = utils.select(data[75], criteria=['Q2 >= 4.0'])


####  --  Longitudinally polarized target --
#
TSA1points = utils.select(data[52], criteria=['FTn == -1'])  # HERMES A_UL
TSApoints = TSA1points + data[54]  # HERMES+CLAS  A_UL
BTSApoints = utils.select(data[53], criteria=['FTn==0'])   # HERMES A_LL
LPpoints = TSApoints + BTSApoints  # total longitudinal target
H_AULpts = TSA1points[:4]
C_AULpts = data[54][:3]
AULpts = H_AULpts + C_AULpts
AULptsOLD = H_AULpts + data[54]  # CLAS points not independent
ALLpts = BTSApoints[:4]

####  --  Transversally polarized target --
#
AUTIpoints = utils.select(data[66], criteria=['FTn==1'])  # \sin\varphi\cos\phi
AUTICSpoints = utils.select(data[66], criteria=['FTn==-1'])  # \cos\varphi\sin\phi
AUTDVCSpoints = data[65]  # HERMES A_UT_DVCS
TPpoints = AUTIpoints + AUTDVCSpoints  # total transversal target
AUTIpts = AUTIpoints[:4]
AUTICSpts = AUTICSpoints[:4]
AUTDVCSpts = AUTDVCSpoints[:4]

# Global combinations
#
GLOpoints = data[31][12:] + data[8] + data[29]  # DM's GLO set
ALTGLOpoints = data[5] + data[25] + data[32][18:]  # KK's CLAS BSA
ALTGLO5points = data[5] + data[8] + data[32][18:]   # DM's CLAS BSA
#UNPpoints = ALTGLOpoints + BSSwpoints + BSDwpoints
UNP5points = ALTGLO5points + BSSwpoints + BSDwpoints
#
GLOall = H1ZEUS[::3] + ALUIpts + BCApts + CLASpts + HApts + AULpts + ALLpts + AUTIpts
GLOfull = (H1ZEUS + ALUIpts + BCApts + CLASKKpts + BSSwpoints + BSDwpoints
            + LPpoints + TPpoints)
GLOfixfullKK = (ALUIpts + BCApts + CLASKKpts + BSSwpoints + BSDwpoints
            + LPpoints + TPpoints)
GLOfixfull = (ALUIpts + BCApts + CLASpts + BSSwpoints + BSDwpoints
            + LPpoints + TPpoints)
# Excluding LP
GLOnoL = H1ZEUS[::3] + ALUIpts + BCApts + CLASpts + HApts + AUTIpts
# Excluding problematic Hall A BSS:
GLOnoBSS = H1ZEUS[::3] + ALUIpts + BCApts + CLASptsOLD + BSDwpoints[::2] + AULptsOLD + ALLpts + AUTIpts
#               12         6         12        4           6              4+6      4         4
GLOnoBSS2 = H1ZEUS + ALUIpts + BCApts + CLASptsOLD + BSDwpoints + AULptsOLD + ALLpts + AUTIpts
GLO15 = H1ZEUS + ALUIpts + BCApts + CLASpts + AULpts + ALLpts + AUTIpts
GLO15new = data[94]+data[95]+data[96]+data[101]+data[102]+data[116]+data[117]
GLO15newuw = data[94]+data[95]+data[96]+data[101]+data[102]+data[120]+data[121]
# Removing dependent CLAS07 data and adding all new 2015 data:
GLO15b = H1ZEUS + ALUIpts + BCApts + CLAS08pts + AULpts + ALLpts + AUTIpts + GLO15new
GLO15buw = H1ZEUS + ALUIpts + BCApts + CLAS08pts + AULpts + ALLpts + AUTIpts + GLO15newuw
unppts = [ALUIpts, BCApts[6:], CLASpts, BSSwpoints[::-2]]
polpts = [TSA1points[:4], data[54], BTSApoints[:4], AUTIpoints[:4], AUTDVCSpoints[:4]]
#
# Independent sets according to target polarization:
#
HERMES_U = BCApts +  ALUIpts 
HERMES_L = TSA1points[:4] + ALLpts
HERMES_T = AUTIpts + AUTICSpts + AUTDVCSpts
#
CLAS_U = CLASKKpts + CLAS14BSApts
CLAS_L = data[54][:3] + CLAS14TSApts + CLAS14BTSApts
#


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
# For ALT last point is overall
L4_ALTI_1 = utils.select(data[74], criteria=['FTn==1'])[:-1]
L4_ALTI_0 = utils.select(data[74], criteria=['FTn==0'])[:-1]
L4_ALTI_m1 = utils.select(data[74], criteria=['FTn==-1'])[:-1]
L4_ALTBHDVCS_0 = utils.select(data[73], criteria=['FTn==0'])[:-1]

bins = zip(L4_ALUI, L4_AC_0, L4_AC_1, L4_AUL, L4_ALL_0, 
        L4_ALL_1, L4_AUTI_1, L4_AUTI_0, L4_AUTI_m1, L4_AUTDVCS,
        L4_ALTI_m1, L4_ALTI_0, L4_ALTI_1, L4_ALTBHDVCS_0)

