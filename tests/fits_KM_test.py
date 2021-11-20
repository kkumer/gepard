"""Testing published KM fits."""

import gepard as g
from pytest import approx, mark

GLOpoints = g.data.dset[31][12:] + g.data.dset[8] + g.data.dset[29]

DVCSpoints = g.data.dset[36] + g.data.dset[37] + g.data.dset[38] + g.data.dset[39] + \
  g.data.dset[40] + g.data.dset[41] + g.data.dset[42] + g.data.dset[43] + \
  g.data.dset[44] + g.data.dset[45]

H1ZEUSpoints = DVCSpoints + g.data.dset[48]

ALTGLO5points = g.data.dset[5] + g.data.dset[8] + g.data.dset[32][18:]   # DM's CLAS BSA
BSDwpoints = g.utils.select(g.data.dset[50], criteria=['FTn == -1'])
BSSwpoints = g.utils.select(g.data.dset[51], criteria=['FTn>=0', 'FTn <= 1'])
UNP5points = ALTGLO5points + BSSwpoints + BSDwpoints

H1ZEUSindependentNEW = g.data.dset[45] + g.data.dset[39] + g.data.dset[63] + g.data.dset[46]
H1ZEUS = H1ZEUSindependentNEW + g.utils.select(g.data.dset[47], criteria=['Q2 >= 4.0'])

ALUIpoints = g.utils.select(g.data.dset[68], criteria=['FTn == -1'])  # HERMES
BCA0points = g.utils.select(g.data.dset[67], criteria=['FTn == 0'])  # HERMES
BCA1points = g.utils.select(g.data.dset[67], criteria=['FTn == 1'])  # HERMES
ALUIpts = ALUIpoints[:6]
BCApts = BCA0points[:6] + BCA1points[:6]

TSA1points = g.utils.select(g.data.dset[52], criteria=['FTn == -1'])  # HERMES A_UL
H_AULpts = TSA1points[:4]
C_AULpts = g.data.dset[54][:3]
AULpts = H_AULpts + C_AULpts

BTSApoints = g.utils.select(g.data.dset[53], criteria=['FTn==0'])   # HERMES A_LL
ALLpts = BTSApoints[:4]

AUTIpoints = g.utils.select(g.data.dset[66], criteria=['FTn==1'])  # \sin\varphi\cos\phi
AUTIpts = AUTIpoints[:4]

CLAS08pts = g.utils.select(g.data.dset[81], criteria=['FTn == -1'])[-3:]

GLO15new = g.data.dset[94]+g.data.dset[95]+g.data.dset[96]+g.data.dset[101]+g.data.dset[102]+g.data.dset[116]+g.data.dset[117]

GLO15b = H1ZEUS + ALUIpts + BCApts + CLAS08pts + AULpts + ALLpts + AUTIpts + GLO15new


par_KM09a = {'tMv': 0.8, 'rS': 1.0, 'alv': 0.43, 'Nv': 1.35,
             'rv': 0.9479, 'NS': 1.5, 'alS': 1.13, 'alpS': 0.15,
             'C': 0.24, 'tNv': 0.0, 'tbv': 3.0, 'bv': 0.4481, 'Mv': 0.8, 'alpv': 0.85,
             'MC': 0.5002, 'MS': 0.7071067811865476, 'bS': 3.0879, 'trv': 0.0}

par_KM09b = {'tMv': 0.8, 'rS': 1.0, 'alv': 0.43, 'Nv': 1.35, 'rv': 1.1064, 'NS': 1.5,
             'alS': 1.13, 'alpS': 0.15, 'C': 6.042, 'tNv': 0.6, 'tbv': 1.5, 'bv': 2.398,
             'Mv': 0.8, 'alpv': 0.85, 'MC': 1.5238, 'MS': 0.7071067811865476,
             'bS': 4.5845, 'trv': 4.8966}

par_KM10 = {'tMv': 0.88386035557556641, 'rS': 1.0, 'rpi': 3.5355742996824659,
            'alv': 0.43, 'Mpi': 0.72673907227274226, 'Nv': 1.35,
            'rv': 0.61981779615528687, 'Nsea': 0.0, 'alS': 1.13, 'alpS': 0.15,
            'C': 8.7771705279055432, 'tNv': 0.6, 'tbv': 2.0496515976221952,
            'bv': 0.40375381018570683, 'Mv': 3.999996566773552, 'alpv': 0.85,
            'MC': 0.9746453693853796, 'bS': 2.0, 'trv': 7.7592155354071064,
            'PG': 2.0, 'PS': 2.0,
            'AL0S': 1.1575060246398083, 'SKEWG': 0.0, 'AL0G': 1.2473167010704711,
            'SKEWS': 0.0, 'M02G': 0.7,
            'DELM2S': 0.0, 'DELM2G': 0.0,
            # This is duplicated, consolidate please
            'MS': 0.707, 'M02S': 0.51318802297356658,
            # These are renamed already:
            'ns':  0.15203911208796006, 'al0s': 1.1575060246398083, 'alps': 0.15,
            'ms2': 0.51318802297356658, 'secs': 0.27800486286895021, 'ng': 0.5,
            'this': -0.12999801345190121,
            'al0g': 1.247316701070471, 'alpg': 0.15, 'mg2': 0.7,
            'secg': -2.787206139384161, 'thig': 0.89633675296304127,
            'kaps': 0.0, 'kapg': 0.0}

par_KM10b = {'tMv': 0.8, 'rS': 1.0, 'rpi': 4.0201, 'alv': 0.43, 'Nsea': 0.0,
             'Nv': 1.35, 'rv': 0.8081, 'Mpi': 1.5369, 'alS': 1.13, 'alpS': 0.15,
             'C': 5.4259, 'tNv': 0.6, 'bS': 2.0, 'bv': 0.7706, 'Mv': 0.8,
             'tbv': 1.0, 'alpv': 0.85, 'MC': 1.3305, 'MS': 0.707, 'trv': 3.2931,
             'EAL0G': 1.1, 'ESECS': 0.0, 'EDELM2S': 0.0, 'EPS': 2.0, 'ETHIS': 0.0,
             'ESECG': 0.0, 'EPG': 2.0, 'EDELM2G': 0.0, 'PS': 2.0, 'EALPG': 0.15,
             'EKAPG': 0.0, 'ESKEWG': 0.0, 'M02S': 0.49754317018981614,
             'EALPS': 0.15, 'EKAPS': 0.0, 'DELB': 0.0, 'ESKEWS': 0.0, 'SKEWS': 0.0,
             'ETHIG': 0.0, 'EM02G': 0.7, 'EAL0S': 1.0, 'DELM2S': 0.0, 'EM02S': 1.0,
             'SKEWG': 0.0, 'PG': 2.0, 'DELM2G': 0.0, 'ns': 0.15203911208796006,
             'al0s': 1.1575060246398083, 'alps': 0.15, 'al0g': 1.247316701070471,
             'secs': -0.4600511871918772, 'this': 0.09351798951979662,
             'alpg': 0.15, 'mg2': 0.7, 'secg': -2.5151319493485427,
             'ms2': 0.49754317018981614,
             'thig': 0.8915757559175185, 'kaps': 0.0, 'kapg': 0.0}

par_KM15 = {'tMv': 3.992860161655587, 'rS': 1.0, 'alv': 0.43, 'tal': 0.43,
            'Mpi': 3.999999852084612, 'Nv': 1.35, 'MS': 0.707, 'rv': 0.918393047884448,
            'Nsea': 0.0, 'alS': 1.13, 'rpi': 2.6463144464701536,  'alpS': 0.15,
            'C': 2.7678681812890016, 'tNv': 0.6, 'bS': 2.0, 'tbv': 0.4000000003259146,
            'bv': 0.40000206775282354, 'Mv': 0.7892078914770488, 'alpv': 0.85,
            'talp': 0.85, 'MC': 1.2040519863464496, 'trv': 0.881085721967267,
            'EPS': 2.0,  'EPG': 2.0, 'EM02G': 0.7,  'EDELM2S': 0.0,
            'SKEWG': 0.0, 'SKEWS': 0.0, 'PS': 2.0, 'EALPG': 0.15,
            'M02D': 1.0, 'M02G': 0.7, 'EALPS': 0.15, 'PG': 2.0, 'EDELM2G': 0.0,
            'AL0D': 0.5, 'DELM2S': 0.0,  'DELM2G': 0.0, 'EAL0G': 1.1,
            'EAL0S': 1.0, 'ND': 1.0, 'ESKEWG': 0.0, 'ETHIG': 0.0,
            'ESECS': 0.0, 'ETHIS': 0.0, 'ESECG': 0.0, 'ESKEWS': 0.0, 'DELB': 0.0,
            'EM02S': 1.0, 'ALPD': 1.0,
            # This is duplicated, consolidate please
            'MS': 0.707, 'M02S': 0.4818827240886959,
            # These are renamed already:
            'ns':  0.15203911208796006, 'al0s': 1.1575060246398083, 'alps': 0.15,
            'ms2': 0.4818827240886959, 'secs': 1.0707825621025808, 'ng': 0.5,
            'this': -0.36618269477432946,
            'al0g': 1.247316701070471, 'alpg': 0.15, 'mg2': 0.7,
            'secg': -2.990809378821039, 'thig': 0.9052207712570559,
            'kaps': 0.0, 'kapg': 0.0}


def test_KM09a():
    """Test model: KM09a."""
    m = g.model.ModelDR()
    m.parameters.update(par_KM09a)
    th = g.theory.hotfixedBMK(m)
    chisq = th.chisq(GLOpoints)
    assert chisq == approx(32.044069303618073)


def test_KM09b():
    """Test model: KM09b."""
    m = g.model.ModelDR()
    m.parameters.update(par_KM09b)
    th = g.theory.hotfixedBMK(m)
    pts = GLOpoints + g.data.dset[30]
    chisq = th.chisq(pts)
    assert chisq == approx(33.36747338543438)


@mark.slow
def test_KM10():
    """Test model: KM10."""
    fit_gpd = g.model.PWNormGPD()
    mMB = g.model.MellinBarnesModel(gpds=fit_gpd)
    mDR = g.model.ComptonModelDRPP()
    m = g.model.HybridDipole(mMB, mDR)
    m.parameters.update(par_KM10)
    th = g.theory.BM10(m)
    pts = H1ZEUSpoints + UNP5points
    chisq = th.chisq(pts)
    # assert_almost_equal(chisq, 135.85499940324056)
    # # from BSA to ALUI
    # assert_almost_equal(chisq, 135.72627340347722)
    # # unclear why the diff
    assert chisq == approx(135.7650869105709)


@mark.skip('KM10a not yet transferred')
def test_KM10a():
    """Test model: KM10a."""
    th = db['KM10a']
    pts = H1ZEUSpoints+ALTGLOpoints
    chisq = th.chisq(pts)[0]
    #assert_almost_equal(chisq, 132.14636420551949)
    # from BSA to ALUI
    assert chisq == approx(129.18281370844684)


@mark.slow
def test_KM10b():
    """Test model: KM10b."""
    fit_gpd = g.model.PWNormGPD()
    mMB = g.model.MellinBarnesModel(gpds=fit_gpd)
    mDR = g.model.ComptonModelDRPP()
    m = g.model.HybridKelly(mMB, mDR)
    m.parameters.update(par_KM10b)
    th = g.theory.BM10(m)
    pts = DVCSpoints+GLOpoints+g.data.dset[30]
    chisq = th.chisq(pts)
    assert chisq == approx(115.54198973827977)


@mark.skip('KMM12 not yet transferred')
def test_KMM12():
    """Test model: KMM12."""
    th = db['KMM12']
    AULptsOLD = TSA1points[:4] + data[54]  # data[54] are not independent
    pts = H1ZEUS + ALUIpts + BCApts + CLASptsOLD + BSDwpoints\
            + AULptsOLD + ALLpts + AUTIpts + BSSwpoints
    chisq = th.chisq(pts)[0]
    # For difference to EPJA review, see comment for KM15 below
    assert chisq == approx(123.53321520926985)


@mark.slow
def test_KM15():
    """Test model: KM15."""
    fit_gpd = g.model.PWNormGPD()
    mMB = g.model.MellinBarnesModel(gpds=fit_gpd)
    mDR = g.model.ComptonModelDRPP()
    m = g.model.HybridKelly(mMB, mDR)
    m.parameters.update(par_KM15)
    th = g.theory.BM10tw2(m)
    pts = GLO15b
    #H1ZEUS + ALUIpts + BCApts + CLAS08pts + AULpts + ALLpts + AUTIpts + GLO15new
    chisq = th.chisq(pts)
    # Number is slightly different from EPJA review due to change in
    # treatment of uncertainties of Hall A 2015 data 
    # in commit 4033d7accb5bea7c371e4145343ef650cb38b6b9
    # on Feb 14 2017 
    assert chisq == approx(245.82778194537178)
