C     ****h* gepard/init.F
C  FILE DESCRIPTION
C    initialization routine
C
C    $Id$
C     *******

C     ****s* init.F/INIT
C  NAME
C        INIT  --  basic initialization
C  DESCRIPTION
C     Puts values of abscissas and weights of Mellin-Barnes
C     integration contour on common blocks, as well as
C     corresponding values of Wilson coefficients and
C     anomalous dimensions
C  SYNOPSIS

      SUBROUTINE INIT

C  CHILDREN
C      BETAF, INTEGRAF, WgammaV*F, WcV*F, CDVCSF, MSBARF, BIGCNSF
C  PARENTS
C      FIT, FITINIT, cffHInternal, cffEInternal
C  SOURCE
C

      IMPLICIT NONE
      INTEGER NINTGMAX, K
      INTEGER L, ORD, SEC, PIDK
      INTEGER EFFACC, NBNDMAX
      DOUBLE PRECISION NFD
      DOUBLE COMPLEX J, Z, EPH, FSHU
      DOUBLE COMPLEX HS1, HS2, HS3, HS4
      DOUBLE COMPLEX CF2, CFL
      DOUBLE COMPLEX WCTMP(0:2,4)
      DOUBLE COMPLEX CLNGAMMA
      PARAMETER ( NINTGMAX = 12 )
*     PARAMETER ( NINTGMAX = 17 )
      DOUBLE PRECISION DOWN(NINTGMAX+1)
#include "header.F"
      PARAMETER (NBNDMAX = 1)
      DOUBLE PRECISION BND(NBNDMAX+1), MTSEXP(4), MTSAUX(MTINDMAX)

*  Initialization of MTS array, carrying all needed -t values

*     [1] MT=0 for forward (DIS) case
      MTS(0) = 0.0D0

*     [2] MT's occuring as abscissae of Gaussian integration of PARSIGMA(MT)
*     Integration is done on the interval:
      DATA BND /0.0D0, 1.0D0/ 
*     Integrand PARSIGMA(MT) is nicely-behaved function and we
*     integrate it by equidistant Gauss-Legendre integration.
*     Effective numerical accuracy used for this integration is:
      EFFACC = ACC - SPEED + 1
*     which means that resulting number of integration points is:
      NMTS = 2**EFFACC
*     Calculating abscissas and weights
      CALL INTEGRAF(EFFACC, 1, BND, NBNDMAX, MTSAUX, MTWG)
*     Puting abscissas into MTS array
      DO 4 K = 1, NMTS
        MTS(K) = MTSAUX(K)
  4   CONTINUE

*     [3] MT's occuring in experimental data
      NMTSEXP = 4
      DATA MTSEXP / 0.1d0, 0.3d0, 0.5d0, 0.8d0 /
      DO 5 K = 1, NMTSEXP
        MTS(NMTS+K) = MTSEXP(K)
  5   CONTINUE

*  (Re-)initialization of NQS array, carrying number of
*  Q2 values for which WCE has been calculated already

      DO 6 K = -4, PIDMAX
        NQS(K) = 0
  6   CONTINUE


*  Initialization of N array, carrying all needed j values on MB contour

*   (2**ACC)-point integration is to be performed on each
*   interval defined by taking points (1, 1 + SPEED,
*   1 + 2*SPEED, ...) from the list DOWN
      
      DATA DOWN 
     & / 0.D0, 0.01D0, 0.08D0, 0.15D0, 0.3D0, 
     &         0.5D0, 1.0D0, 1.5D0, 2.0D0, 
     &         4.0D0, 6.0D0, 8.0D0, 10.0D0 /


*  Same as Pegasus:
*     DATA DOWN  / 0.D0,  0.1D0, 0.3D0, 0.6D0, 1.0D0, 1.6D0, 2.4D0,
*    1             3.5D0, 5.0D0, 7.0D0, 10.D0, 14.D0, 19.D0, 25.D0,
*    2             32.D0, 40.D0, 50.D0, 63.D0 /


*   Gaussian abscissae and weights:

      CALL INTEGRAF(ACC, SPEED, DOWN, NINTGMAX, Y, WG)


      NPTS = 2**ACC * NINTGMAX / SPEED

*   Calculating actual Mellin-Barnes contour points from Gaussian abscissae

      EPH = EXP ( DBLCMPLX(0.D0, PHI) )
      DO 10 K = 1, NPTS
        N(1,K) = C + 1.0d0 + Y(K) * EPH 
        N(2,K) = C + 3.0d0 + Y(K) * EPH 
 10     N(3,K) = C + 5.0d0 + Y(K) * EPH 

*   ----  Initialization of common blocks ----


*   1. Initialization of QCD beta function coefficients

      CALL BETAF

*   Making everything for three MB contours, first shifted (j+4),
*   then shifted (j+2),
*   and then the original one, so that SEC-independent stuff
*   ends up with original values
      DO 100 SEC = 3, 1, -1
*   Now looping over MB contour points and initializing
      DO 100 K = 1, NPTS

      Z = N(SEC, K)
      J = N(SEC, K) - 1

*   2. ADACF initialization

*    NB: adacf routines want double precision NF
      NFD = DBLE(NF)

*   2.a Harmonic sums

      S1 =  HS1(Z)
      IF (P .GE. 1) THEN
        S2 =  HS2(Z)
        IF (P .GE. 2) THEN
            S3 = HS3(Z)
            S4 = HS4(Z)
        END IF
      END IF
              
*   2.b Anomalous dimensions matrices: LO, NLO, and NNLO

      CALL WgammaVQQ0F(NFD, Z, GAM(SEC, K, 0, 1, 1))
      CALL WgammaVQG0F(NFD, Z, GAM(SEC, K, 0, 1, 2))
      CALL WgammaVGQ0F(NFD, Z, GAM(SEC, K, 0, 2, 1))
      CALL WgammaVGG0F(NFD, Z, GAM(SEC, K, 0, 2, 2))
      CALL WgammaVNSP0F(NFD, Z, GAMNS(SEC, K, 0, 0))
      CALL WgammaVNSP0F(NFD, Z, GAMNS(SEC, K, 0, 1))

      IF (P .GE. 1) THEN
        CALL WgammaVQQ1F(NFD, Z, GAM(SEC, K, 1, 1, 1))
        CALL WgammaVQG1F(NFD, Z, GAM(SEC, K, 1, 1, 2))
        CALL WgammaVGQ1F(NFD, Z, GAM(SEC, K, 1, 2, 1))
        CALL WgammaVGG1F(NFD, Z, GAM(SEC, K, 1, 2, 2))
        CALL WgammaVNSP1F(NFD, Z, GAMNS(SEC, K, 1, 0))
        CALL WgammaVNSM1F(NFD, Z, GAMNS(SEC, K, 1, 1))

        IF (P .GE. 2) THEN
          CALL WgammaVQQ2F(NFD, Z, GAM(SEC, K, 2, 1, 1))
          CALL WgammaVQG2F(NFD, Z, GAM(SEC, K, 2, 1, 2))
          CALL WgammaVGQ2F(NFD, Z, GAM(SEC, K, 2, 2, 1))
          CALL WgammaVGG2F(NFD, Z, GAM(SEC, K, 2, 2, 2))
          CALL WgammaVNSP2F(NFD, Z, GAMNS(SEC, K, 2, 0))
c         CALL WgammaVNSM2F(NFD, Z, GAMNS(SEC, K, 2, 1))
        END IF
      END IF



*   2.c DIS Wilson coefficients: LO, NLO, and NNLO

      CDIS1(K, 0, 1) = (1.0d0, 0.0d0)
      CDIS1(K, 0, 2) = (0.0d0, 0.0d0)
      CDIS2(K, 0, 1) = (1.0d0, 0.0d0)
      CDIS2(K, 0, 2) = (0.0d0, 0.0d0)
      CDISNS1(K, 0) = (1.0d0, 0.0d0)

      IF (P .GE. 1) THEN
        CALL WcVF2Q1F(NFD, Z, CF2)
        CALL WcVFLQ1F(NFD, Z, CFL)
        CDIS1(K, 1, 1) = CF2 -  CFL
        CDIS2(K, 1, 1) = CF2
        CALL WcVF2G1F(NFD, Z, CF2)
        CALL WcVFLG1F(NFD, Z, CFL)
        CDIS1(K, 1, 2) = CF2 -  CFL
        CDIS2(K, 1, 2) = CF2
        CALL WcVF2NSP1F(NFD, Z, CF2)
        CALL WcVFLNSP1F(NFD, Z, CFL)
        CDISNS1(K, 1) = CF2 -  CFL

        IF (P .GE. 2) THEN
          CALL WcVF2Q2F(NFD, Z, CF2)
          CALL WcVFLQ2F(NFD, Z, CFL)
          CDIS1(K, 2, 1) = CF2 -  CFL
          CALL WcVF2G2F(NFD, Z, CF2)
          CALL WcVFLG2F(NFD, Z, CFL)
          CDIS1(K, 2, 2) = CF2 -  CFL
          CALL WcVF2NSP2F(NFD, Z, CF2)
          CALL WcVFLNSP2F(NFD, Z, CFL)
          CDISNS1(K, 2) = CF2 -  CFL
        END IF
      END IF

*   2.d  "Shuvaev" factor
*     2^(J+1) Gamma(5/2+J) / Gamma(3/2) / Gamma(3+J)
      FSHU = 2.0d0**(J+1.d0)
     &        * EXP(CLNGAMMA(2.5d0 + J) - CLNGAMMA(3.0d0 + J))
     &        / 0.886226925452758014d0


*   3. Wilson coefficients for all needs

*  --- first initialize everything to zero ...

        DO 20 L = 1, 2
        DO 20 ORD = 1, 2
        WCTMP(ORD, L) = (0.0d0, 0.0d0)
        DO 20 PIDK = -4, PIDMAX
        WC(PIDK, SEC, K, ORD, L) = (0.0d0, 0.0d0)
 20     CONTINUE

*  --- ... then change what needs to be changed

*     PIDK = -1: Quark PDF or GPD at eta=0 
        WC(-1, SEC, K, 0, 1) = (1.0d0, 0.0d0)
        WC(-1, SEC, K, 0, 3) = (1.0d0, 0.0d0)
        WC(-1, SEC, K, 0, 4) = (1.0d0, 0.0d0)

*     PIDK = -2: Gluon PDF or GPD at eta=0 
        WC(-2, SEC, K, 0, 2) = (1.0d0, 0.0d0)

*     PIDK = -3: Quark GPD at eta=x 
        WC(-3, SEC, K, 0, 1) = (1.0d0, 0.0d0) * FSHU
        WC(-3, SEC, K, 0, 3) = (1.0d0, 0.0d0) * FSHU
        WC(-3, SEC, K, 0, 4) = (1.0d0, 0.0d0) * FSHU

*     PIDK = -4: Gluon GPD at eta=x 
        WC(-4, SEC, K, 0, 2) = FSHU * 2.0d0 / (3.0d0 + J)


*     PIDK = 0: DIS F2 
*       here only SEC=1 case makes sense
        IF (SEC .EQ. 1) THEN
        CALL WCCSBAR(0, SEC, K, WCTMP)
          DO 30 L = 1,2
          DO 30 ORD = 0, P
 30         WC(0, SEC, K, ORD, L) = WCTMP(ORD, L) 
        END IF

*     PIDK = 1: DVCS 
        IF (SCHEME .EQ. 'CSBAR') THEN
          CALL WCCSBAR(1, SEC, K, WCTMP)
        ELSE
          CALL WCMSBAR(SEC, K, WCTMP)
        END IF
        DO 40 L = 1,4
        DO 40 ORD = 0, P
 40       WC(1, SEC, K, ORD, L) = WCTMP(ORD, L) * FSHU

*     PIDK = 2: DVMP
        CALL WCDVMP(SEC, K, WCTMP)
*       WC = T_jk from Eq. (3.61c) in big DVMP draft

        DO 60 L = 1,4
        DO 60 ORD = 0, P
        IF (L .EQ. 2) THEN
*         --- gluons get additional factor
          WC(2, SEC, K, ORD, L) = ( WCTMP(ORD, L) * FSHU * 3.d0
     &                * 2.d0/CF/(J+3.d0) )
        ELSE
*         --- quarks, singlet and non-singlet
          WC(2, SEC, K, ORD, L) = WCTMP(ORD, L) * FSHU * 3.d0
        END IF
 60     CONTINUE


100   CONTINUE

      RETURN
      END 
C     ***
