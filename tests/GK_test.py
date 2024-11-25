"""Testing the Goloskokov-Kroll GPD/CFF model.

Numbers for comparison are obtained by running
PARTONS version 84bb0db 2023-05-22

on xml file containing:

<task service="GPDService" method="computeSingleKinematic" storeInDB="0">

      <!-- Define GPD kinematics -->
      <kinematics type="GPDKinematic">
         <param name="x" value="0.1" />
         <param name="xi" value="0.25" />
         <param name="t" value="-0.1" />
         <param name="MuF2" value="4." />
         <param name="MuR2" value="4." />
      </kinematics>

      <!-- Define physics assumptions -->
       <computation_configuration>

         <!-- Select GPD model -->
         <module type="GPDModule" name="GPDGK16">

        <!-- Select GPD evolution model -->
            <module type="GPDEvolutionModule" name="GPDEvolutionApfel">

           <!-- Select pQCD order -->
               <param name="qcd_order_type" value="LO" />

           <!-- Select alpha_s model -->
               <module type="RunningAlphaStrongModule" name="RunningAlphaStrongApfel">
                  <param name="qcd_order_type" value="LO" />
          <param name="alphasRef" value="0.118" />
          <param name="muRef" value="91.1876" />
          <param name="thresholds" value="0 0 0 1.e6 1.e6 1.e6" />
               </module>

           <!-- Select model defining number of active flavors -->
               <module type="ActiveFlavorsThresholdsModule" name="ActiveFlavorsThresholdsConstant">

          <!-- In this ActiveFlavorsThresholdsModule model we put the value by hand -->
          <param name="nFlavors" value="3" />

               </module>

            </module>

         </module>

      </computation_configuration>

   </task>


with appropriate changes of the value of x.

Note that only input scale values are tested.

up, um below correspond to PARTONS' u(+), u(-)

"""

import gepard as g
from pytest import approx, fixture, mark

dglap = (0.3, 0.25, -0.1, 4.)
erbl = (0.1, 0.25, -0.1, 4.)
merbl = (-0.1, 0.25, -0.1, 4.)
mdglap = (-0.3, 0.25, -0.1, 4.)

tilde_erbl = (0.1, 0.25, -0.3, 4.)
mtilde_erbl = (-0.1, 0.25, -0.3, 4.)

pt = g.DataPoint({'xi': 0.25, 't': -0.1, 'Q2': 4.})

@fixture
def thGK():
    th = g.cff.GoloskokovKrollCFF()
    return th

def test_GK_GPD_H(thGK):
    """Test GK model's GPD H"""
    # -- u quark ---
	# DGLAP
    Hxp = thGK.Huval(*dglap) + thGK.Hudsea(*dglap)
    Hxm = thGK.Huval(*mdglap) + thGK.Hudsea(*mdglap)
    up = Hxp - Hxm
    um = Hxp + Hxm
    assert up == approx(2.20540466070162)
    assert um == approx(2.01002230133684)
	# ERBL (x=-0.1 !)
    Hxp = thGK.Huval(*merbl) + thGK.Hudsea(*merbl)
    Hxm = thGK.Huval(*erbl) + thGK.Hudsea(*erbl)
    up = Hxp - Hxm
    um = Hxp + Hxm
    assert up == approx(-3.39228378371742)
    assert um == approx(4.96820106695312)
    # -- s quark ---
	# ERBL (x=-0.1 !)
    Hxp = thGK.Hs(*merbl)
    Hxm = thGK.Hs(*erbl)
    up = Hxp - Hxm
    um = Hxp + Hxm
    assert up == approx(-0.795100066988823, 6)
    assert um == approx(0)
    # -- gluon ---
	# ERBL (x=-0.1 !)
    Hxp = thGK.Hg(*merbl)
    Hxm = thGK.Hg(*erbl)
    assert Hxp == approx(1.25938307923902, 6)
    assert Hxm == approx(1.25938307923902, 6)
    
def test_GK_GPD_E(thGK):
    """Test GK model's GPD E"""
    # -- u quark ---
	# DGLAP
    Exp = thGK.Euval(*dglap) + thGK.Esea(*dglap)
    Exm = thGK.Euval(*mdglap) + thGK.Esea(*mdglap)
    up = Exp - Exm
    um = Exp + Exm
    assert up == approx(1.00472696853448)
    assert um == approx(1.1553414061214)
    # -- d quark ---
	# ERBL (x=0.1)
    Exp = thGK.Edval(*erbl) + thGK.Esea(*erbl)
    Exm = thGK.Edval(*merbl) + thGK.Esea(*merbl)
    up = Exp - Exm
    um = Exp + Exm
    assert up == approx(-3.13166152120034, 6)
    assert um == approx(-6.124044762973, 6)
    # -- s quark ---
	# ERBL (x=0.1)
    Exp = thGK.Es(*erbl)
    Exm = thGK.Es(*merbl)
    up = Exp - Exm
    um = Exp + Exm
    assert up == approx(-1.31445903393, 6)
    assert um == approx(0)
    # -- gluon ---
	# ERBL (x=0.1)
    Exp = thGK.Eg(*erbl)
    Exm = thGK.Eg(*merbl)
    assert Exp == approx(0.441547825524533, 6)
    assert Exm == approx(0.441547825524533, 6)

def test_GK_GPD_Ht(thGK):
    """Test GK model's GPD H tilde"""
    # -- u quark ---
	# DGLAP
    Hxp = thGK.Htu(*dglap)
    Hxm = thGK.Htu(*mdglap)
    up = Hxp + Hxm
    um = Hxp - Hxm
    assert up == approx(1.21590783919976)
    assert um == approx(1.21590783919976)
    # -- d quark ---
	# DGLAP
    Hxp = thGK.Htd(*dglap)
    Hxm = thGK.Htd(*mdglap)
    up = Hxp + Hxm
    um = Hxp - Hxm
    assert up == approx(-0.393040371782663)
    assert um == approx(-0.393040371782663)
    # -- gluon ---
	# DLAP (x=-0.3 !)
    Htg = thGK.Htg(*mdglap)
    # This is MINUS PARTONS number. There is a likely
    #   bug in PARTONS since this GPD should be odd in x
    assert Htg == approx(-0.240948963124176, 6)

def test_GK_GPD_Et(thGK):
    """Test GK model's GPD E tilde"""
    # -- u quark ---
	# erbl
    Exp = thGK.Etu(*tilde_erbl)
    Exm = thGK.Etu(*mtilde_erbl)
    up = Exp + Exm
    um = Exp - Exm
    assert up == approx(44.9152604654097)
    assert um == approx(6.06306948294471)
    # -- d quark ---
    Exp = thGK.Etd(*tilde_erbl)
    Exm = thGK.Etd(*mtilde_erbl)
    up = Exp + Exm
    um = Exp - Exm
    assert up == approx(-21.2010402326042)
    assert um == approx(1.73230556655563)


def test_GK_CFF_H(thGK):
    """Test GK model: CFF H."""
    imh = thGK.ImH(pt)
    reh = thGK.ReH(pt)
    assert imh == approx(4.95687064454363, 6)
    assert reh == approx(1.47975836538955, 6)

def test_GK_CFF_E(thGK):
    """Test GK model: CFF E."""
    ime = thGK.ImE(pt)
    ree = thGK.ReE(pt)
    assert ime == approx(1.01847574299162, 6)
    assert ree == approx(-0.937841730264802, 6)

def test_GK_CFF_Ht(thGK):
    """Test GK model: CFF Ht."""
    imht = thGK.ImHt(pt)
    reht = thGK.ReHt(pt)
    assert imht == approx(1.84416391852452, 6)
    assert reht == approx(1.67047040922322, 6)

def test_GK_CFF_Et_nopole(thGK):
    """Test GK model: CFF Ht. (no pion pole)"""
    imet = thGK.ImEt(pt)
    reet = thGK.ReEt(pt)
    assert imet == approx(14.080820529834, 6)
    assert reet == approx(34.0325125264094, 6)

def test_GK_CFF_Et_pole(thGK):
    """Test GK model: CFF Ht. (+ pion pole)"""
    ptp = g.DataPoint({'xi': 0.25, 't': -0.3, 'Q2': 4.})
    imet = thGK.ImEt(ptp)
    reet = thGK.ReEt(ptp)
    assert imet == approx(9.76198175159559, 6)
    assert reet == approx(41.1394045960092, 6)
