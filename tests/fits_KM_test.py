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
            'ms': 0.51318802297356658, 'secs': 0.27800486286895021, 'ng': 0.5,
            'this': -0.12999801345190121,
            'al0g': 1.247316701070471, 'alpg': 0.15, 'mg': 0.7,
            'secg': -2.787206139384161, 'thig': 0.89633675296304127,
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
    fit_gpd = g.model.Fit()
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


@mark.skip('KM10 not yet transferred')
def test_KM10a():
    """Test model: KM10a."""
    th = db['KM10a']
    pts = H1ZEUSpoints+ALTGLOpoints
    chisq = th.chisq(pts)[0]
    #assert_almost_equal(chisq, 132.14636420551949)
    # from BSA to ALUI
    assert chisq == approx(129.18281370844684)


@mark.skip('KM10 not yet transferred')
def test_KM10b():
    """Test model: KM10b."""
    th = db['KM10b']
    pts = DVCSpoints+GLOpoints+data[30]
    chisq = th.chisq(pts)[0]
    assert chisq == approx(115.54198973827977)


@mark.skip('KM10 not yet transferred')
def test_KMM12():
    """Test model: KMM12."""
    th = db['KMM12']
    AULptsOLD = TSA1points[:4] + data[54]  # data[54] are not independent
    pts = H1ZEUS + ALUIpts + BCApts + CLASptsOLD + BSDwpoints\
            + AULptsOLD + ALLpts + AUTIpts + BSSwpoints
    chisq = th.chisq(pts)[0]
    # For difference to EPJA review, see comment for KM15 below
    assert chisq == approx(123.53321520926985)


@mark.skip('KM10tw2 not yet transferred')
def test_KM15():
    """Test model: KM15."""
    th = db['KM15']
    pts = GLO15b
    chisq = th.chisq(pts)[0]
    # Number is slightly different from EPJA review due to change in
    # treatment of uncertainties of Hall A 2015 data 
    # in commit 4033d7accb5bea7c371e4145343ef650cb38b6b9
    # on Feb 14 2017 
    assert chisq == approx(245.82778194537178)
