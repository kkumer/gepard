"""Testing model database theories.db"""

import gepard as g
from pytest import approx, mark

GLOpoints = g.data.dset[31][12:] + g.data.dset[8] + g.data.dset[29]

par_KM09a = {'tMv': 0.8, 'rS': 1.0, 'alv': 0.43, 'Nv': 1.35,
             'rv': 0.9479, 'NS': 1.5, 'alS': 1.13, 'alpS': 0.15,
             'C': 0.24, 'tNv': 0.0, 'tbv': 3.0, 'bv': 0.4481, 'Mv': 0.8, 'alpv': 0.85,
             'MC': 0.5002, 'MS': 0.7071067811865476, 'bS': 3.0879, 'trv': 0.0}

par_KM09b = {'tMv': 0.8, 'rS': 1.0, 'alv': 0.43, 'Nv': 1.35, 'rv': 1.1064, 'NS': 1.5,
             'alS': 1.13, 'alpS': 0.15, 'C': 6.042, 'tNv': 0.6, 'tbv': 1.5, 'bv': 2.398,
             'Mv': 0.8, 'alpv': 0.85, 'MC': 1.5238, 'MS': 0.7071067811865476,
             'bS': 4.5845, 'trv': 4.8966}


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


@mark.skip('KM10 not yet transferred')
def test_KM10():
    """Test model: KM10."""
    th = db['KM10']
    pts = H1ZEUSpoints + UNP5points
    chisq = th.chisq(pts)[0]
    #assert_almost_equal(chisq, 135.85499940324056)
    # from BSA to ALUI
    #assert_almost_equal(chisq, 135.72627340347722)
    # unclear why the diff
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
