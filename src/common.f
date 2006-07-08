C     ****h* gepard/common.f
C  FILE DESCRIPTION
C    initialization of common blocks with anomalous dimensions
C    and Wilson coefficients of DIS
C
C    $Id$
C     *******


C     ****s* common.f/COMMONF
C  NAME
C     COMMONF  --   initialization of WGAMMA and WC common blocks
C  DESCRIPTION
C    initialization of common blocks with anomalous dimensions
C    (WGAMMA) and Wilson coefficients (WC) of DIS
C  SYNOPSIS
C     SUBROUTINE COMMONF (Z, NF)
C
C     INTEGER NF
C     DOUBLE COMPLEX Z 
C     DOUBLE COMPLEX S1, S2, S3, S4
C     DOUBLE COMPLEX GAM0(2,2), GAM1(2,2), GAM2(2,2)
C     DOUBLE COMPLEX C0(2), C1(2), C2(2)
C
C     COMMON / HARMONIC /  S1, S2, S3, S4
C     COMMON / WGAMMA   /  GAM0, GAM1, GAM2
C     COMMON / WC       /  C0, C1, C2
C  INPUTS
C           P -- approximation order, which is N^{P}LO
C          NF -- number of active flavours
C           Z -- complex Mellin moment, Z = N = J + 1
C                where J is complex conformal moment
C  PARENTS
C      PARWAVF
C  CHILDREN
C      WgammaV???F, WcVF???F, HS?  --  from adacf library
C  NOTES
C    - initialization of HARMONIC common block is required
C      by adacf library
C
C    - PSI, DPSI_1, DPS_2, DPSI_3 should probably also be
C      initialized here and put in a common block.
C
C    - One idea for speed improvement is to evaluate these
C      on a FIXED GRID of J values and later use quadratures
C      with fixed abscissas.
C  SOURCE
C

      SUBROUTINE COMMONF (K)

      IMPLICIT NONE
      INTEGER K, P, NPTS
      DOUBLE COMPLEX S1, S2, S3, S4
      DOUBLE COMPLEX GAM0(2,2), GAM1(2,2), GAM2(2,2)
      DOUBLE COMPLEX C0(2), C1(2), C2(2)
      DOUBLE COMPLEX F2, FL
      PARAMETER ( NPTS = 32 )
      DOUBLE COMPLEX GAM(NPTS,0:2,2,2), CF2(NPTS,0:2,2), CFL(NPTS,0:2,2)
      DOUBLE COMPLEX HARMS(4,NPTS), N(NPTS)

*   Input common-blocks 

      COMMON / APPROX     /  P
      COMMON / NPOINTS    /  N
      COMMON / VALUES     /  HARMS, GAM, CF2, CFL

*   Output common-blocks

      COMMON / HARMONIC /  S1, S2, S3, S4
      COMMON / WGAMMA   /  GAM0, GAM1, GAM2
      COMMON / WC       /  C0, C1, C2

*   Harmonic sum initialization

      S1 = HARMS(1, K)
      IF (P .GE. 1) THEN
        S2 = HARMS(2, K)
        IF (P .GE. 2) THEN
            S3 = HARMS(3, K)
            S4 = HARMS(4, K)
        END IF
      END IF

*  Initializing anomalous dimensions matrices: LO, NLO, and NNLO
*    NB: adacf routines want double precision NF

      GAM0(1,1) = GAM(K, 0, 1, 1)
      GAM0(1,2) = GAM(K, 0, 1, 2)
      GAM0(2,1) = GAM(K, 0, 2, 1)
      GAM0(2,2) = GAM(K, 0, 2, 2)

      IF (P .GE. 1) THEN
      GAM1(1,1) = GAM(K, 1, 1, 1)
      GAM1(1,2) = GAM(K, 1, 1, 2)
      GAM1(2,1) = GAM(K, 1, 2, 1)
      GAM1(2,2) = GAM(K, 1, 2, 2)

      IF (P .GE. 2) THEN
      GAM2(1,1) = GAM(K, 2, 1, 1)
      GAM2(1,2) = GAM(K, 2, 1, 2)
      GAM2(2,1) = GAM(K, 2, 2, 1)
      GAM2(2,2) = GAM(K, 2, 2, 2)
      END IF
      END IF

*  Initializing Wilson coefficients
*    NB: adacf routines want double precision NF

      C0(1) = (1.0d0, 0.0d0)
      C0(2) = (0.0d0, 0.0d0)

      IF (P .GE. 1) THEN
      F2 = CF2(K, 1, 1)
      FL = CFL(K, 1, 1)
      C1(1) = F2 -  FL
      F2 = CF2(K, 1, 2)
      FL = CFL(K, 1, 2)
      C1(2) = F2 -  FL

      IF (P .GE. 2) THEN
      F2 = CF2(K, 2, 1)
      FL = CFL(K, 2, 1)
      C2(1) = F2 -  FL
      F2 = CF2(K, 2, 2)
      FL = CFL(K, 2, 2)
      C2(2) = F2 -  FL
      END IF
      END IF

      RETURN
      END
C     *******
