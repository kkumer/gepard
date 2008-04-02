C     ****h* gepard/init.f
C  FILE DESCRIPTION
C    initialization routine
C
C    $Id: dvcs.f 11 2006-07-12 07:50:28Z kuk05260 $
C     *******

C     ****s* init.f/INIT
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
C  BUGS
C       For large xi>0.7 integration should be extended to larger Y
C  SOURCE
C

      IMPLICIT NONE
      INTEGER NINTGMAX, K, K1, K2, K3
      INTEGER L, ORD
      INTEGER EFFACC, NBNDMAX
      DOUBLE PRECISION NFD
      DOUBLE COMPLEX J, Z, EPH
      DOUBLE COMPLEX HS1, HS2, HS3, HS4
      DOUBLE COMPLEX CF2, CFL
      DOUBLE COMPLEX C1F2(2)
      DOUBLE COMPLEX BIGC0(2), BIGC1(2), BIGC2(2)
      DOUBLE COMPLEX BIGCTMP(0:2,2)
      DOUBLE COMPLEX CLNGAMMA, PREFACT
      PARAMETER ( NINTGMAX = 12 )
      DOUBLE PRECISION DOWN(NINTGMAX+1)
      INCLUDE 'header.f'
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



*  Initialization of N array, carrying all needed j values on MB contour

*   (2**ACC)-point integration is to be performed on each
*   interval defined by taking points (1, 1 + SPEED,
*   1 + 2*SPEED, ...) from the list DOWN
      
      DATA DOWN 
     & / 0.D0, 0.01D0, 0.025D0, 0.067D0, 0.18D0, 
     &         0.5D0, 1.3D0, 3.7D0, 10.D0, 
     &         25.D0, 0.67D2, 1.8D2, 5.0D2 /


*   For fast calculation it's better to compress integration region
*   The value 1.3 below is optimized for small xi of cca. 10^-4
! (but this was OK and needed when NINTG was 8, and integration went to
!  10. It should be changed, if it's needed at all. FIXME)
!      IF ( SPEED .GE. 4) THEN
!        DOWN (NINTG + 1) = 1.3d0
!      END IF

*   Gaussian abscissae and weights:

      CALL INTEGRAF(ACC, SPEED, DOWN, NINTGMAX, Y, WG)


*   If you suspect numerical problems and want to extend Mellin-Barnes
*   integration fom 0..10  to 0..500, put 4 -> 0 below
*   FIXME: this should be automatic

      NPTS = 2**ACC * (NINTGMAX - 4) / SPEED

*   Calculating actual Mellin-Barnes contour points from Gaussian abscissae

      EPH = EXP ( COMPLEX(0.D0, PHI) )
      DO 10 K = 1, NPTS
 10     N(K) = C + 1.0d0 + Y(K) * EPH 

*   ----  Initialization of common blocks ----

*   1. Initialization of QCD beta function coefficients

      CALL BETAF

*   Now looping over MB contour points and initializing
      DO 100 K = 1, NPTS

      Z = N(K)
      J = N(K) - 1

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


      IF ( FFTYPE(:7) .EQ. 'SINGLET' ) THEN

*   -----------  SINGLET  ------------------
              
*   2.b Anomalous dimensions matrices: LO, NLO, and NNLO

      CALL WgammaVQQ0F(NFD, Z, GAM(K, 0, 1, 1))
      CALL WgammaVQG0F(NFD, Z, GAM(K, 0, 1, 2))
      CALL WgammaVGQ0F(NFD, Z, GAM(K, 0, 2, 1))
      CALL WgammaVGG0F(NFD, Z, GAM(K, 0, 2, 2))

      IF (P .GE. 1) THEN
        CALL WgammaVQQ1F(NFD, Z, GAM(K, 1, 1, 1))
        CALL WgammaVQG1F(NFD, Z, GAM(K, 1, 1, 2))
        CALL WgammaVGQ1F(NFD, Z, GAM(K, 1, 2, 1))
        CALL WgammaVGG1F(NFD, Z, GAM(K, 1, 2, 2))

        IF (P .GE. 2) THEN
          CALL WgammaVQQ2F(NFD, Z, GAM(K, 2, 1, 1))
          CALL WgammaVQG2F(NFD, Z, GAM(K, 2, 1, 2))
          CALL WgammaVGQ2F(NFD, Z, GAM(K, 2, 2, 1))
          CALL WgammaVGG2F(NFD, Z, GAM(K, 2, 2, 2))
        END IF
      END IF


*   2.c DIS Wilson coefficients ("small" c)

      CDIS1(K, 0, 1) = (1.0d0, 0.0d0)
      CDIS1(K, 0, 2) = (0.0d0, 0.0d0)
      CDIS2(K, 0, 1) = (1.0d0, 0.0d0)
      CDIS2(K, 0, 2) = (0.0d0, 0.0d0)

      IF (P .GE. 1) THEN
        CALL WcVF2Q1F(NFD, Z, CF2)
        CALL WcVFLQ1F(NFD, Z, CFL)
        CDIS1(K, 1, 1) = CF2 -  CFL
        CDIS2(K, 1, 1) = CF2
        CALL WcVF2G1F(NFD, Z, CF2)
        CALL WcVFLG1F(NFD, Z, CFL)
        CDIS1(K, 1, 2) = CF2 -  CFL
        CDIS2(K, 1, 2) = CF2

        IF (P .GE. 2) THEN
          CALL WcVF2Q2F(NFD, Z, CF2)
          CALL WcVFLQ2F(NFD, Z, CFL)
          CDIS1(K, 2, 1) = CF2 -  CFL
          CALL WcVF2G2F(NFD, Z, CF2)
          CALL WcVFLG2F(NFD, Z, CFL)
          CDIS1(K, 2, 2) = CF2 -  CFL
        END IF
      END IF

*   3. "Big C' Wilson coefficients

      IF ((SCHEME .EQ. 'CSBAR') .OR. (PROCESS(:3) .EQ. 'DIS')) THEN
          CALL CDVCSF(K, BIGCTMP)
      ELSE IF (SCHEME(:3) .EQ. 'MSB') THEN
          CALL MSBARF(K, BIGCTMP)
      END IF

*     Writing this to BIGC or BIGCF2 common blocks.
*     "Big C" Wilson coefficients of DVCS (in BIGC) have to be multiplied by
*     Gamma(5/2+J) / Gamma(3+J)

      IF ( PROCESS(:3) .EQ. 'DVC' ) THEN
        DO 30 L = 1,2
        DO 30 ORD = 0, P
          BIGC(K, ORD, L) = BIGCTMP(ORD, L) * 
     &             EXP(CLNGAMMA(2.5d0 + J) - CLNGAMMA(3.0d0 + J))
 30     CONTINUE
      ELSE
*        -- DIS --
        DO 40 L = 1,2
        DO 40 ORD = 0, P
          BIGCF2(K, ORD, L) = BIGCTMP(ORD, L)
 40     CONTINUE
      END IF


      ELSE
*  -------------  NON - SINGLET  ----------------

*   2.bNS Anomalous dimensions matrices: LO, NLO, and NNLO

      CALL WgammaVNSP0F(NFD, Z, GAMNS(K, 0))
      IF (P .GE. 1) THEN
        CALL WgammaVNSP1F(NFD, Z, GAMNS(K, 1))
        IF (P .GE. 2) THEN
          CALL WgammaVNSP2F(NFD, Z, GAMNS(K, 2))
        END IF
      END IF

*   2.cNS DIS Wilson coefficients ("small" c)

      CDISNS1(K, 0) = (1.0d0, 0.0d0)
      IF (P .GE. 1) THEN
        CALL WcVF2NSP1F(NFD, Z, CF2)
        CALL WcVFLNSP1F(NFD, Z, CFL)
        CDISNS1(K, 1) = CF2 -  CFL
        IF (P .GE. 2) THEN
          CALL WcVF2NSP2F(NFD, Z, CF2)
          CALL WcVFLNSP2F(NFD, Z, CFL)
          CDISNS1(K, 2) = CF2 -  CFL
        END IF
      END IF

*   3.NS "Big C' Wilson coefficients

*     non-singlet DIS is not implemented yet, so here we call
*     BIGCNSF which directly writes to BIGCNS common block

      CALL BIGCNSF(K)

*     "Big C" Wilson coefficients of DVCS have to be multiplied by
*     Gamma(5/2+J) / Gamma(3+J)

      IF ( PROCESS(:3) .EQ. 'DVC' ) THEN
        DO 50 ORD = 0, P
          BIGCNS(K, ORD) = BIGCNS(K, ORD) * 
     &             EXP(CLNGAMMA(2.5d0 + J) - CLNGAMMA(3.0d0 + J))
 50     CONTINUE
      END IF

      END IF

100   CONTINUE

      RETURN
      END 
C     ***
