C     ****h* gepard/init.f
C  FILE DESCRIPTION
C    initialization routine
C
C    $Id: dvcs.f 11 2006-07-12 07:50:28Z kuk05260 $
C     *******

C     ****s* init.f/INIT
C  NAME
C        INIT  --  initialization
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
C     NGAM  -- singlet anomalous dimensions
C   NGAMNS  -- non-singlet anomalous dimensions
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
      INTEGER ACCMAX, NGAUSSMAX, NGAUSS, NINTG, K, K1, K2, K3
      DOUBLE PRECISION NFD, SUMM, DIFF, YI
      DOUBLE COMPLEX J, Z, EPH
      DOUBLE COMPLEX HS1, HS2, HS3, HS4
      DOUBLE COMPLEX CF2, CFL
      DOUBLE COMPLEX C1F2(2)
      DOUBLE COMPLEX BIGC0(2), BIGC1(2), BIGC2(2)
      DOUBLE COMPLEX CLNGAMMA, PREFACT
      PARAMETER ( NGAUSSMAX = 64, ACCMAX = 6, NINTG = 12 )
      DOUBLE PRECISION ABSC(ACCMAX,NGAUSSMAX), WGHT(ACCMAX,NGAUSSMAX)
      DOUBLE PRECISION DOWN(NINTG+1), UP(NINTG)
      INCLUDE 'header.f'

*   Output common-blocks


*   NGAUSS-point integration is to be performed on each
*   interval defined by taking points (1, 1 + SPEED,
*   1 + 2*SPEED, ...) from the list DOWN
      
      NGAUSS = 2**ACC

      DATA DOWN 
     & / 0.D0, 0.01D0, 0.025D0, 0.067D0, 0.18D0, 
     &         0.5D0, 1.3D0, 3.7D0, 10.D0, 
     &         20.D0, 50.D0,  1.0D2,  5.0D2 /

      CALL GAUSS (ABSC, WGHT)


*   For fast calculation it's better to compress integration region
*   The value 1.3 below is optimized for small xi of cca. 10^-4
! (but this was OK and needed when NINTG was 8, and integration went to
!  10. It should be changed, if it's needed at all. FIXME)
!      IF ( SPEED .GE. 4) THEN
!        DOWN (NINTG + 1) = 1.3d0
!      END IF

*    NB: adacf routines want double precision NF
      NFD = DBLE(NF)

*   Mellin-Barnes contour points are C + YI * EPH
      EPH = EXP ( COMPLEX(0.D0, PHI) )

*   Setting up upper limits of integration intervals

      DO 10 K = 1, NINTG, SPEED
 10     UP(K) = DOWN(K+SPEED)

*   Abscissas N(K) and weights WG(K) for the whole contour. Note that
*   the factor (upper limit - lower limit)/2 is put into the weights

      K = 0
      DO 20 K2 = 1, NINTG, SPEED
        SUMM = UP(K2) + DOWN(K2) 
        DIFF = UP(K2) - DOWN(K2) 
      DO 20 K3 = 1, NGAUSS
        K = K + 1
        YI = (DIFF * ABSC(ACC, K3) + SUMM) * 0.5D0
        N(K) = (C + 1.0D0) + YI * EPH
        WG(K) = 0.5D0 * DIFF * WGHT(ACC, K3) 
 20   CONTINUE

*   For everything apart from ND evolution one should use 0->4 below
*   and increase speed by 30 %

      NPTS = NGAUSS * (NINTG-0) / SPEED


*   1. Initialization of QCD beta function coefficients

      CALL BETAF

*   Now looping over contour points and initializing

      DO 100 K = 1, NPTS

      Z = N(K)
      J = N(K) - 1

*   2. Auxilliary stuff not needed outside of this subroutine

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

      CALL WgammaVQQ0F(NFD, Z, GAM0(1,1))
      CALL WgammaVQG0F(NFD, Z, GAM0(1,2))
      CALL WgammaVGQ0F(NFD, Z, GAM0(2,1))
      CALL WgammaVGG0F(NFD, Z, GAM0(2,2))

      IF (P .GE. 1) THEN
        CALL WgammaVQQ1F(NFD, Z, GAM1(1,1))
        CALL WgammaVQG1F(NFD, Z, GAM1(1,2))
        CALL WgammaVGQ1F(NFD, Z, GAM1(2,1))
        CALL WgammaVGG1F(NFD, Z, GAM1(2,2))

        IF (P .GE. 2) THEN
          CALL WgammaVQQ2F(NFD, Z, GAM2(1,1))
          CALL WgammaVQG2F(NFD, Z, GAM2(1,2))
          CALL WgammaVGQ2F(NFD, Z, GAM2(2,1))
          CALL WgammaVGG2F(NFD, Z, GAM2(2,2))
        END IF
      END IF

      IF ( ANSATZ(:2) .EQ. 'NS' ) THEN
        CALL WgammaVNSP0F(NFD, Z, GAMNS0)
        CALL WgammaVNSP1F(NFD, Z, GAMNS1)
!        CALL WgammaVNS2PF(NFD, Z, GAMNS2)
      END IF


*   2.c DIS Wilson coefficients

      C0(1) = (1.0d0, 0.0d0)
      C0(2) = (0.0d0, 0.0d0)

      IF (P .GE. 1) THEN
        CALL WcVF2Q1F(NFD, Z, CF2)
        CALL WcVFLQ1F(NFD, Z, CFL)
        C1(1) = CF2 -  CFL
        C1F2(1) = CF2
        CALL WcVF2G1F(NFD, Z, CF2)
        CALL WcVFLG1F(NFD, Z, CFL)
        C1(2) = CF2 -  CFL
        C1F2(2) = CF2

        IF (P .GE. 2) THEN
          CALL WcVF2Q2F(NFD, Z, CF2)
          CALL WcVFLQ2F(NFD, Z, CFL)
          C2(1) = CF2 -  CFL
          CALL WcVF2G2F(NFD, Z, CF2)
          CALL WcVFLG2F(NFD, Z, CFL)
          C2(2) = CF2 -  CFL
        END IF
      END IF

*   3. Stuff used also outside

*   3.a Anomalous dimensions
      DO 30 K1 = 1,2
      DO 30 K2 = 1,2
        NGAM(K, 0, K1, K2) = GAM0(K1, K2)
        NGAM(K, 1, K1, K2) = GAM1(K1, K2)
 30     NGAM(K, 2, K1, K2) = GAM2(K1, K2)

      IF ( ANSATZ(:2) .EQ. 'NS' ) THEN
        NGAMNS(K, 0) = GAMNS0
        NGAMNS(K, 1) = GAMNS1
!        NGAMNS(K, 2) = GAMNS2
      END IF

*   3.b "Big C" Wilson coefficients of DVCS [multiplied by
*           PREFACT = Gamma(5/2+J) / Gamma(3+J)]

      IF (SCHEME .EQ. 'CSBAR') THEN
          CALL CDVCSF (J, BIGC0, BIGC1, BIGC2, 'DVCS')
      ELSE IF (SCHEME(:3) .EQ. 'MSB') THEN
          CALL MSBARF (J, BIGC0, BIGC1)
      END IF

      PREFACT = EXP(CLNGAMMA(2.5d0 + J) - CLNGAMMA(3.0d0 + J))

      DO 40 K1 = 1,2
        BIGC(K, 0, K1) = BIGC0(K1) * PREFACT
        BIGC(K, 1, K1) = BIGC1(K1) * PREFACT
        BIGC(K, 2, K1) = BIGC2(K1) * PREFACT
 40   CONTINUE

*   3.c "Big C" Wilson coefficients of DIS 

*     First put C_1->C_2 on WC block
      DO 50 K1 = 1,2
 50     C1(K1) = C1F2(K1)
      
      CALL CDVCSF (J, BIGC0, BIGC1, BIGC2, 'DIS')

      DO 60 K1 = 1,2
        BIGCF2(K, 0, K1) = BIGC0(K1)
 60     BIGCF2(K, 1, K1) = BIGC1(K1)

100   CONTINUE

      RETURN
      END 
C     ***
