C     ****h* gepard/init.f
C  FILE DESCRIPTION
C    initialization routine
C
C    $Id: dvcs.f 11 2006-07-12 07:50:28Z kuk05260 $
C     *******

C     ****s* init.f/INIT
C  NAME
C        INIT  --  initialization (singlet)
C  DESCRIPTION
C     Puts values of abscissas and weights of Mellin-Barnes
C     integration contour on common blocks, as well as
C     corresponding values of Wilson coefficients and
C     anomalous dimensions
C  SYNOPSIS
C     SUBROUTINE INIT
C  OUTPUT
C     BIGC  -- values of DVCS Wilson coefficients for
C              singlet CFF form factor \mathcal{H}
C     BIGC  -- values of DVCS Wilson coefficients for
C              singlet CFF form factor \mathcal{H}
C              but on a shifted contour C_{j+2}
C      GAM  -- singlet anomalous dimensions
C    GAMNS  -- non-singlet anomalous dimensions
C   BIGCF2  -- values of DIS Wilson coefficients for
C              singlet form factor F2
C  IDENTIFIERS
C       SPEED -- speed of evaluations
C           P -- approximation order, which is N^{P}LO
C           C -- intersection point of Mellin- Barnes integration 
C                path with real axis
C         PHI -- angle between Mellin-Barnes contour and Re(J) axis
C        GAM? -- ?=0,1,2  singlet anomalous dimensions
C      GAMNS? -- ?=0,1,2  non-singlet anomalous dimensions
C          C? -- ?=0,1,2  DIS Wilson coefficients
C  CHILDREN
C      BETAF, WgammaV*F, WcV*F, CDVCSF, MSBARF
C  PARENTS
C      AUXTEST, TEST, RADCORR, SCALEDEP, FIT
C  BUGS
C       For large xi>0.7 integration should be extended to larger Y
C  SOURCE
C

      SUBROUTINE INIT
* 
      IMPLICIT NONE
      INTEGER NINTGMAX, K, K1, K2, K3
      INTEGER L, ORD, SEC
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


*   If you suspect numerical problems and want extend Mellin-Barnes
*   integration fom 0..10  to 0..500, put 4 -> 0 below

      NPTS = 2**ACC * (NINTGMAX - 4) / SPEED

*   Calculating actual Mellin-Barnes contour points from Gaussian abscissae

      EPH = EXP ( COMPLEX(0.D0, PHI) )
      DO 10 K = 1, NPTS
        N(1,K) = C + 1.0d0 + Y(K) * EPH 
 10     N(2,K) = C + 3.0d0 + Y(K) * EPH 

*   ----  Initialization of common blocks ----

*   0. Initialization of charge factor

      IF ( FFTYPE(:7) .EQ. 'SINGLET' ) THEN
        IF (NF .EQ. 3) THEN
           CHARGEFAC =  2.0D0 / 9.0D0
        ELSE IF (NF .EQ. 4) THEN
           CHARGEFAC =  5.0D0 / 18.0D0
        ELSE
           CALL ERROR ('GeParD', ' INIT',
     & 'NF is not integer equal to 3 or 4!                          ',
     & 5, 3)
        END IF
      ELSE
*         -- NONSINGLET --
        IF (NF .EQ. 3) THEN
           CHARGEFAC =  1.0D0 / 9.0D0
        ELSE IF (NF .EQ. 4) THEN
           CHARGEFAC =  1.0D0 / 6.0D0
        ELSE
           CALL ERROR ('GeParD', ' INIT',
     & 'NF is not integer equal to 3 or 4!                          ',
     & 5, 3)
        END IF
      END IF

*         -- when calculating just GPDs --

      IF ( ( SCHEME(:4) .EQ. 'ZERO' ) .OR. 
     &     ( SCHEME(:4) .EQ. 'TRAJ' ) ) THEN
        CHARGEFAC = 1.0d0
      END IF

*   1. Initialization of QCD beta function coefficients

      CALL BETAF

*   Making everything for two MB contours, first shifted (j+2),
*   and then the original one, so that SEC-independent stuff
*   ends up with original values
      DO 100 SEC = 2, 1, -1
*   Now looping over MB contour points and initializing
      DO 100 K = 1, NPTS

      Z = N(SEC,K)
      J = N(SEC,K) - 1

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

*****************  SINGLET  *****************
              
*   2.b Anomalous dimensions matrices: LO, NLO, and NNLO

      CALL WgammaVQQ0F(NFD, Z, GAM(SEC, K, 0, 1, 1))
      CALL WgammaVQG0F(NFD, Z, GAM(SEC, K, 0, 1, 2))
      CALL WgammaVGQ0F(NFD, Z, GAM(SEC, K, 0, 2, 1))
      CALL WgammaVGG0F(NFD, Z, GAM(SEC, K, 0, 2, 2))

      IF (P .GE. 1) THEN
        CALL WgammaVQQ1F(NFD, Z, GAM(SEC, K, 1, 1, 1))
        CALL WgammaVQG1F(NFD, Z, GAM(SEC, K, 1, 1, 2))
        CALL WgammaVGQ1F(NFD, Z, GAM(SEC, K, 1, 2, 1))
        CALL WgammaVGG1F(NFD, Z, GAM(SEC, K, 1, 2, 2))

        IF (P .GE. 2) THEN
          CALL WgammaVQQ2F(NFD, Z, GAM(SEC, K, 2, 1, 1))
          CALL WgammaVQG2F(NFD, Z, GAM(SEC, K, 2, 1, 2))
          CALL WgammaVGQ2F(NFD, Z, GAM(SEC, K, 2, 2, 1))
          CALL WgammaVGG2F(NFD, Z, GAM(SEC, K, 2, 2, 2))
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

      IF (SCHEME .EQ. 'CSBAR') THEN
          CALL CDVCSF(SEC, K, BIGCTMP)
      ELSE IF (SCHEME(:3) .EQ. 'MSB') THEN
          CALL MSBARF(SEC, K, BIGCTMP)
      ELSE IF ( (SCHEME .EQ. 'ZEROQ') .OR. (SCHEME .EQ. 'TRAJQ') ) THEN
        BIGCTMP(0, 1) = (1.0d0, 0.0d0)
        BIGCTMP(0, 2) = (0.0d0, 0.0d0)
        DO 15 L = 1, 2
        DO 15 ORD = 1, 2
 15     BIGCTMP(ORD, L) = (0.0d0, 0.0d0)
      ELSE IF ( (SCHEME .EQ. 'ZEROG') .OR. (SCHEME .EQ. 'TRAJG') .OR.
     &          (SCHEME .EQ. 'EVOLG') ) THEN
        BIGCTMP(0, 1) = (0.0d0, 0.0d0)
*  Gluons need special treatment here and for DVCS additional factor of xi later
        IF ( SCHEME .EQ. 'TRAJG' ) THEN
*         GPD(x, eta=x,t)                                              
          BIGCTMP(0, 2) = 2.0d0 / (3.0d0 + J)
        ELSE
*         PDF(x) and GPD(x, eta=0,t)                                              
          BIGCTMP(0, 2) = (1.0d0,0.0d0)
        END IF
        DO 20 L = 1, 2
        DO 20 ORD = 1, 2
 20     BIGCTMP(ORD, L) = (0.0d0, 0.0d0)
      END IF

*     Writing this to BIGC or BIGCF2 common blocks.
*     "Big C" Wilson coefficients of DVCS (in BIGC) have to be multiplied by
*     2^(J+1) Gamma(5/2+J) / Gamma(3/2) / Gamma(3+J)
*     BTW,  Gamma(3/2) = 0.88622...

      DO 30 L = 1,2
      DO 30 ORD = 0, P
        IF ( PROCESS(:3) .EQ. 'DIS' ) THEN
*          -- F2 or PDFs --
          BIGCF2(K, ORD, L) = BIGCTMP(ORD, L)
        ELSE IF ( SCHEME(:4) .EQ. 'ZERO' ) THEN
*          -- GPDs at eta=0 trajectory --
          BIGC(SEC, K, ORD, L) = BIGCTMP(ORD, L)
        ELSE
*          -- CFFs or GPDs at eta=x trajectory --
          BIGC(SEC, K, ORD, L) = BIGCTMP(ORD, L) * 2.0d0**(J+1.d0)
     &           * EXP(CLNGAMMA(2.5d0 + J) - CLNGAMMA(3.0d0 + J))
     &           / 0.886226925452758014d0
        END IF
 30   CONTINUE


      ELSE
*****************  NON - SINGLET  *****************

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

*     "Big C" Wilson coefficients of DVCS (in BIGC) have to be multiplied by
*     2^(J+1) Gamma(5/2+J) / Gamma(3/2) / Gamma(3+J)
*     BTW,  Gamma(3/2) = 0.88622...

      IF ( PROCESS(:3) .EQ. 'DVC' ) THEN
        DO 50 ORD = 0, P
          BIGCNS(K, ORD) = BIGCNS(K, ORD) * 2.0d0**(J+1.d0) 
     &          *  EXP(CLNGAMMA(2.5d0 + J) - CLNGAMMA(3.0d0 + J))
     &          / 0.886226925452758014d0
 50     CONTINUE
      END IF

      END IF

100   CONTINUE

      RETURN
      END 
C     ***
