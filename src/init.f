*
      SUBROUTINE INIT
* 
      IMPLICIT NONE
      INTEGER NGAUSS, NINTG, NPTS, K, K1, K2, K3, NF, PMAX, SPEED
      INTEGER NPTSMAX
      DOUBLE PRECISION NFD, PI, C, PHI, SUMM, DIFF, YI
      DOUBLE PRECISION RF2, RR2
      DOUBLE COMPLEX J, Z, EPH
      DOUBLE COMPLEX S1, S2, S3, S4
      DOUBLE COMPLEX HS1, HS2, HS3, HS4
      DOUBLE COMPLEX F2, FL
      DOUBLE COMPLEX GAM0(2,2), GAM1(2,2), GAM2(2,2)
      DOUBLE COMPLEX C0(2), C1(2), C2(2)
      DOUBLE COMPLEX BIGC0(2), BIGC1(2), BIGC2(2)
      DOUBLE COMPLEX CLNGAMMA
      PARAMETER ( NGAUSS = 8, NINTG = 8, NPTSMAX = 64 )
      PARAMETER ( PI = 3.1415 92653 58979 D0 )
      CHARACTER SCHEME*5, ANSATZ*6
*
      DOUBLE PRECISION ABSCISSAS(NGAUSS), WEIGHTS(NGAUSS)
      DOUBLE PRECISION DOWN(NINTG+1), UP(NINTG)
      DOUBLE PRECISION Y(NPTSMAX), WG(NPTSMAX)
      DOUBLE COMPLEX N(NPTSMAX)
      DOUBLE COMPLEX BIGC(NPTSMAX,0:2,2), NGAM(NPTSMAX,0:2,2,2)

*   Input common-blocks
      COMMON / LABELS   /  SCHEME, ANSATZ
      COMMON / INITPAR  /  SPEED, PMAX

*   Output common-blocks

*     - Mellin-Barnes integration contour points
      COMMON / POINTS   /  Y, WG
      COMMON / CPHI     /  C, PHI
      COMMON / NPOINTS  /  N
      
*     - Values on a particular contour point
      COMMON / HARMONIC /  S1, S2, S3, S4
      COMMON / WGAMMA   /  GAM0, GAM1, GAM2
      COMMON / WC       /  C0, C1, C2

*     - Values on the whole contour
      COMMON / BIGC     /  BIGC
      COMMON / NGAM     /  NGAM


*   Abscissas and weights for 8 point Gauss quadrature 
*   according to Abramowitz and Stegun Eq. (25.4.30)

       DATA ABSCISSAS
     &  /-0.96028 98564 97536 D0, -0.79666 64774 13627 D0, 
     &   -0.52553 24099 16329 D0, -0.18343 46424 95650 D0, 
     &    0.18343 46424 95650 D0,  0.52553 24099 16329 D0,
     &    0.79666 64774 13627 D0,  0.96028 98564 97536 D0/

       DATA WEIGHTS
     &  / 0.10122 85362 90376 D0,  0.22238 10344 53374 D0, 
     &    0.31370 66458 77887 D0,  0.36268 37833 78362 D0, 
     &    0.36268 37833 78362 D0,  0.31370 66458 77887 D0,
     &    0.22238 10344 53374 D0,  0.10122 85362 90376 D0/

*   NGAUSS = 8 point integration is to be performed on each
*   interval defined by taking points (1, 1 + SPEED,
*   1 + 2*SPEED, ...) from the list DOWN

      DATA DOWN 
     & / 0.D0, 0.01D0, 0.025D0, 0.067D0, 0.18D0, 
     &         0.5D0, 1.3D0, 3.7D0, 10.D0 /

*   For fast calculation it's better to compress integration region
*   The value 1.3 below is optimized for small xi of cca. 10^-4
      IF ( SPEED .GE. 4) THEN
        DOWN (NINTG + 1) = 1.3d0
      END IF

*   Read fixed initialization parameters
      OPEN (UNIT = 61, FILE = "INIT.DAT", STATUS = "OLD")
      READ (61, *) SPEED
      READ (61, *) PMAX
      CLOSE (61)

*    (FIXME: Should be determined elsewhere and propagated here!)
      RF2 = 1.0d0
      RR2 = 1.0d0
      NF = 3

*    NB: adacf routines want double precision NF
      NFD = DBLE(NF)

*   Parameters of the Mellin-Barnes contour

      C   = 0.5D0
      PHI = 3.D0/4.D0 * PI
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
        YI = (DIFF * ABSCISSAS(K3) + SUMM) * 0.5D0
        N(K) = (C + 1) + YI * EPH
        WG(K) = 0.5D0 * DIFF * WEIGHTS(K3) 
 20   CONTINUE

*    After the above loop, total number of points on the contour is
*    (NB: Integer division!)

      NPTS = NPTSMAX / SPEED


*   1. Initialization of QCD beta function coefficients

      CALL BETAF

*   Now looping over contour points and initializing

      DO 100 K = 1, NPTS

      Z = N(K)
      J = N(K) - 1

*   2. Auxilliary stuff not needed outside of this subroutine

*   2.a Harmonic sums

      S1 =  HS1(Z)
      IF (PMAX .GE. 1) THEN
        S2 =  HS2(Z)
        IF (PMAX .GE. 2) THEN
            S3 = HS3(Z)
            S4 = HS4(Z)
        END IF
      END IF

*   2.b Anomalous dimensions matrices: LO, NLO, and NNLO

      CALL WgammaVQQ0F(NFD, Z, GAM0(1,1))
      CALL WgammaVQG0F(NFD, Z, GAM0(1,2))
      CALL WgammaVGQ0F(NFD, Z, GAM0(2,1))
      CALL WgammaVGG0F(NFD, Z, GAM0(2,2))

      IF (PMAX .GE. 1) THEN
        CALL WgammaVQQ1F(NFD, Z, GAM1(1,1))
        CALL WgammaVQG1F(NFD, Z, GAM1(1,2))
        CALL WgammaVGQ1F(NFD, Z, GAM1(2,1))
        CALL WgammaVGG1F(NFD, Z, GAM1(2,2))

        IF (PMAX .GE. 2) THEN
          CALL WgammaVQQ2F(NFD, Z, GAM2(1,1))
          CALL WgammaVQG2F(NFD, Z, GAM2(1,2))
          CALL WgammaVGQ2F(NFD, Z, GAM2(2,1))
          CALL WgammaVGG2F(NFD, Z, GAM2(2,2))
        END IF
      END IF

*   2.c DIS Wilson coefficients

      C0(1) = (1.0d0, 0.0d0)
      C0(2) = (0.0d0, 0.0d0)

      IF (PMAX .GE. 1) THEN
        CALL WcVF2Q1F(NFD, Z, F2)
        CALL WcVFLQ1F(NFD, Z, FL)
        C1(1) = F2 -  FL
        CALL WcVF2G1F(NFD, Z, F2)
        CALL WcVFLG1F(NFD, Z, FL)
        C1(2) = F2 -  FL

        IF (PMAX .GE. 2) THEN
          CALL WcVF2Q2F(NFD, Z, F2)
          CALL WcVFLQ2F(NFD, Z, FL)
          C2(1) = F2 -  FL
          CALL WcVF2G2F(NFD, Z, F2)
          CALL WcVFLG2F(NFD, Z, FL)
          C2(2) = F2 -  FL
        END IF
      END IF

*   3. Stuff used also outside

*   3.a Anomalous dimensions
      DO 30 K1 = 1,2
      DO 30 K2 = 1,2
        NGAM(K, 0, K1, K2) = GAM0(K1, K2)
        NGAM(K, 1, K1, K2) = GAM1(K1, K2)
 30     NGAM(K, 2, K1, K2) = GAM2(K1, K2)

*   3.b "Big C" Wilson coefficients of DVCS

      IF (SCHEME .EQ. 'CSBAR') THEN
          CALL CDVCSF (NF, J, RF2, RR2, BIGC0, BIGC1, BIGC2)
      ELSE IF (SCHEME .EQ. 'MSBAR') THEN
          CALL MSBARF (NF, J, RF2, RR2, BIGC0, BIGC1)
      END IF

      DO 40 K1 = 1,2
        BIGC(K, 0, K1) = BIGC0(K1) *
     &       EXP(CLNGAMMA(2.5d0 + J)) / EXP(CLNGAMMA(3.0d0 + J))
        BIGC(K, 1, K1) = BIGC1(K1) *
     &       EXP(CLNGAMMA(2.5d0 + J)) / EXP(CLNGAMMA(3.0d0 + J))
        BIGC(K, 2, K1) = BIGC2(K1) *
     &       EXP(CLNGAMMA(2.5d0 + J)) / EXP(CLNGAMMA(3.0d0 + J))
 40   CONTINUE

100   CONTINUE

      RETURN
      END 
