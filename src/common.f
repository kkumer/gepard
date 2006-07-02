C     ****h* gepard/common.f
C  FILE DESCRIPTION
C    initialization of common blocks with anomalous dimensions
C    and Wilson coefficients of DIS
C
C    $Id: common.f,v 1.1 2006-04-12 19:11:15+02 kkumer Exp kkumer $
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

      SUBROUTINE COMMONF (Z, NF)

      IMPLICIT NONE
      INTEGER NF, P
      DOUBLE COMPLEX Z 
      DOUBLE COMPLEX S1, S2, S3, S4
      DOUBLE COMPLEX GAM0(2,2), GAM1(2,2), GAM2(2,2)
      DOUBLE COMPLEX C0(2), C1(2), C2(2)
      DOUBLE COMPLEX HS1, HS2, HS3, HS4
      DOUBLE COMPLEX F2, FL

*   Input common-blocks 

      COMMON / APPROX     /  P

*   Output common-blocks

      COMMON / HARMONIC /  S1, S2, S3, S4
      COMMON / WGAMMA   /  GAM0, GAM1, GAM2
      COMMON / WC       /  C0, C1, C2

*   Harmonic sum initialization

      S1 = HS1(Z)
      IF (P .GE. 1) THEN
        S2 = HS2(Z)
        IF (P .GE. 2) THEN
            S3 = HS3(Z)
            S4 = HS4(Z)
        END IF
      END IF

*  Initializing anomalous dimensions matrices: LO, NLO, and NNLO
*    NB: adacf routines want double precision NF

      CALL WgammaVQQ0F(DBLE(NF), Z, GAM0(1,1))
      CALL WgammaVQG0F(DBLE(NF), Z, GAM0(1,2))
      CALL WgammaVGQ0F(DBLE(NF), Z, GAM0(2,1))
      CALL WgammaVGG0F(DBLE(NF), Z, GAM0(2,2))

      IF (P .GE. 1) THEN
      CALL WgammaVQQ1F(DBLE(NF), Z, GAM1(1,1))
      CALL WgammaVQG1F(DBLE(NF), Z, GAM1(1,2))
      CALL WgammaVGQ1F(DBLE(NF), Z, GAM1(2,1))
      CALL WgammaVGG1F(DBLE(NF), Z, GAM1(2,2))

      IF (P .GE. 2) THEN
      CALL WgammaVQQ2F(DBLE(NF), Z, GAM2(1,1))
      CALL WgammaVQG2F(DBLE(NF), Z, GAM2(1,2))
      CALL WgammaVGQ2F(DBLE(NF), Z, GAM2(2,1))
      CALL WgammaVGG2F(DBLE(NF), Z, GAM2(2,2))
      END IF
      END IF

*  Initializing Wilson coefficients
*    NB: adacf routines want double precision NF

      C0(1) = (1.0d0, 0.0d0)
      C0(2) = (0.0d0, 0.0d0)

      IF (P .GE. 1) THEN
      CALL WcVF2Q1F(DBLE(NF), Z, F2)
      CALL WcVFLQ1F(DBLE(NF), Z, FL)
      C1(1) = F2 -  FL
      CALL WcVF2G1F(DBLE(NF), Z, F2)
      CALL WcVFLG1F(DBLE(NF), Z, FL)
      C1(2) = F2 -  FL

      IF (P .GE. 2) THEN
      CALL WcVF2Q2F(DBLE(NF), Z, F2)
      CALL WcVFLQ2F(DBLE(NF), Z, FL)
      C2(1) = F2 -  FL
      CALL WcVF2G2F(DBLE(NF), Z, F2)
      CALL WcVFLG2F(DBLE(NF), Z, FL)
      C2(2) = F2 -  FL
      END IF
      END IF

      RETURN
      END
C     *******
