*
      SUBROUTINE INIT
* 
      IMPLICIT NONE
      INTEGER NGAUSS, NINTG, NPTS, K, K2, K3, NF, P
      DOUBLE PRECISION PI, C, PHI, SUMM, DIFF, YI
      DOUBLE COMPLEX Z, EPH
      DOUBLE COMPLEX HS1, HS2, HS3, HS4
      PARAMETER ( NGAUSS = 8, NINTG = 4, NPTS = NGAUSS*NINTG )
      PARAMETER ( PI = 3.1415 92653 58979 D0 )
*
      DOUBLE PRECISION ABSCISSAS(NGAUSS), WEIGHTS(NGAUSS)
      DOUBLE PRECISION DOWN(NINTG+1), UP(NINTG)
      DOUBLE PRECISION Y(NPTS), WG(NPTS)
      DOUBLE COMPLEX GAM(NPTS,0:2,2,2), CF2(NPTS,0:2,2), CFL(NPTS,0:2,2)
      DOUBLE COMPLEX HARMS(4,NPTS), N(NPTS)

*   Output common-blocks

      COMMON / POINTS     /  Y, WG
      COMMON / NPOINTS    /  N
      COMMON / VALUES     /  HARMS, GAM, CF2, CFL
      COMMON / CPHI       /  C, PHI


*   Abscissas and weights for 8 point Gauss quadrature 
*   according to Abramowitz and Stegun Eq. (25.4.30)
*
       DATA ABSCISSAS
     &  /-0.96028 98564 97536 D0, -0.79666 64774 13627 D0, 
     &   -0.52553 24099 16329 D0, -0.18343 46424 95650 D0, 
     &    0.18343 46424 95650 D0,  0.52553 24099 16329 D0,
     &    0.79666 64774 13627 D0,  0.96028 98564 97536 D0/
*
       DATA WEIGHTS
     &  / 0.10122 85362 90376 D0,  0.22238 10344 53374 D0, 
     &    0.31370 66458 77887 D0,  0.36268 37833 78362 D0, 
     &    0.36268 37833 78362 D0,  0.31370 66458 77887 D0,
     &    0.22238 10344 53374 D0,  0.10122 85362 90376 D0/
*
*   8 point integration is to be performed between following points
*
      DATA DOWN / 0.D0, 0.18D0, 0.5D0, 1.34D0, 3.6D0 /

*   Initialization of QCD beta function coefficients

      CALL BETAF

*   Parameters
      NF = 3
*   Just a local value of P for this subroutine!
      P = 2

* 
*   Parameters of the Mellin inversion contour
*
      C   = 0.5 D0
      PHI = 3.0D0/4.D0 * PI
      EPH = EXP ( COMPLEX(0.0d0, PHI) )
*
*   Setting upper limits of integrations
*
      DO 10 K = 1, NINTG
 10     UP(K) = DOWN(K+1)

*   Abscissas N(K) and weights WG(K) for the whole contour. Note that
*   the factor (upper limit - lower limit)/2 is put into the weights

      K = 0
      DO 20 K2 = 1, NINTG
        SUMM = UP(K2) + DOWN(K2) 
        DIFF = UP(K2) - DOWN(K2) 
      DO 20 K3 = 1, NGAUSS
        K = K + 1
        YI = (DIFF * ABSCISSAS(K3) + SUMM) * 0.5D0
        N(K) = (C + 1) + YI * EPH
        WG(K) = 0.5D0 * DIFF * WEIGHTS(K3) 
 20   CONTINUE


*   Putting values of harmonic summs, anomalous dimensions and Wilson coefficients
*   on integration points
      DO 30 K = 1, NPTS

      Z = N(K)

*   Harmonic sum initialization

      HARMS(1,K) = HS1(Z)
      IF (P .GE. 1) THEN
        HARMS(2,K) = HS2(Z)
        IF (P .GE. 2) THEN
            HARMS(3,K) = HS3(Z)
            HARMS(4,K) = HS4(Z)
        END IF
      END IF

*  Initializing anomalous dimensions matrices: LO, NLO, and NNLO
*    NB: adacf routines want double precision NF

      CALL WgammaVQQ0F(DBLE(NF), Z, GAM(K,0,1,1))
      CALL WgammaVQG0F(DBLE(NF), Z, GAM(K,0,1,2))
      CALL WgammaVGQ0F(DBLE(NF), Z, GAM(K,0,2,1))
      CALL WgammaVGG0F(DBLE(NF), Z, GAM(K,0,2,2))

      IF (P .GE. 1) THEN
      CALL WgammaVQQ1F(DBLE(NF), Z, GAM(K,1,1,1))
      CALL WgammaVQG1F(DBLE(NF), Z, GAM(K,1,1,2))
      CALL WgammaVGQ1F(DBLE(NF), Z, GAM(K,1,2,1))
      CALL WgammaVGG1F(DBLE(NF), Z, GAM(K,1,2,2))

      IF (P .GE. 2) THEN
      CALL WgammaVQQ2F(DBLE(NF), Z, GAM(K,2,1,1))
      CALL WgammaVQG2F(DBLE(NF), Z, GAM(K,2,1,2))
      CALL WgammaVGQ2F(DBLE(NF), Z, GAM(K,2,2,1))
      CALL WgammaVGG2F(DBLE(NF), Z, GAM(K,2,2,2))
      END IF
      END IF

*  Initializing Wilson coefficients
*    NB: adacf routines want double precision NF

      CF2(K,0,1) = (1.0d0, 0.0d0)
      CF2(K,0,2) = (0.0d0, 0.0d0)
      CFL(K,0,1) = (0.0d0, 0.0d0)
      CFL(K,0,2) = (0.0d0, 0.0d0)

      IF (P .GE. 1) THEN
      CALL WcVF2Q1F(DBLE(NF), Z, CF2(K,1,1))
      CALL WcVFLQ1F(DBLE(NF), Z, CFL(K,1,1))
      CALL WcVF2G1F(DBLE(NF), Z, CF2(K,1,2))
      CALL WcVFLG1F(DBLE(NF), Z, CFL(K,1,2))

      IF (P .GE. 2) THEN
      CALL WcVF2Q2F(DBLE(NF), Z, CF2(K,2,1))
      CALL WcVFLQ2F(DBLE(NF), Z, CFL(K,2,1))
      CALL WcVF2G2F(DBLE(NF), Z, CF2(K,2,2))
      CALL WcVFLG2F(DBLE(NF), Z, CFL(K,2,2))
      END IF
      END IF

 30   CONTINUE

      RETURN
      END 
