C     ****h* gepard/cdvcs.f
C  FILE DESCRIPTION
C    calculation of  Wilson coefficients for DVCS in CSbar scheme
C    according to KMKPS06 paper
C
C    $Id$
C     *******


C     ****s* cdvcs.f/CDVCSF
C  NAME
C     CDVCSF  --   Wilson coefficients (in CSbar scheme)
C  DESCRIPTION
C    calculates Wilson coefficients for DVCS in CSbar scheme
C    according to KMKPS06 paper. Also used for DIS where
C    CSbar=MSbar, and for special SCHEME='UNITC'.
C  SYNOPSIS
C     SUBROUTINE CDVCSF (J, BIGC0, BIGC1, BIGC2)
C
C     DOUBLE COMPLEX J, BIGC0(2), BIGC1(2), BIGC2(2)
C     CHARACTER PROCESS*4
C  INPUTS
C           J -- conformal moment
C  OUTPUT
C       BIGC0 -- vector of two C^{(0)} (quark and gluon) Wilson 
C                coefficients (trivially equal to (1,0))
C       BIGC1 -- vector  C^{(1)} 
C       BIGC2 -- vector  C^{(2)} 
C  IDENTIFIERS
C       BETABLK, WGAMMA, WC -- common blocks with beta function coefficients
C       of QCD, and J-th moments of anomalous dimensions and Wilson 
C       coefficients of DIS, initialized by subroutine INIT
C  PARENTS
C      INIT
C  CHILDREN
C      VECMAT, HS1, HS2
C  SOURCE
C

      SUBROUTINE CDVCSF (J, BIGC0, BIGC1, BIGC2, PROCESS)

      IMPLICIT NONE
      DOUBLE COMPLEX J, BIGC0(2), BIGC1(2), BIGC2(2)
      CHARACTER PROCESS*4
      INTEGER SPEED, ACC, P, NF
      DOUBLE PRECISION AS0, RF2, RR2
      CHARACTER SCHEME*5, ANSATZ*6
      INTEGER NFMIN, NFMAX, K
      DOUBLE PRECISION LRF2, LRR2
      DOUBLE PRECISION BETA0, BETA1, BETA2, BETA3
      DOUBLE COMPLEX GAM0(2,2), GAM1(2,2), GAM2(2,2)
      DOUBLE COMPLEX HS1, HS2
      DOUBLE COMPLEX S1, S2
      DOUBLE COMPLEX C0(2), C1(2), C2(2), F2, FL
      DOUBLE COMPLEX VM00(2), VM01(2), VM10(2), VM000(2)
      PARAMETER (NFMIN = 3, NFMAX = 6)

*   Input common-blocks

      COMMON / PARINT /  SPEED, ACC, P, NF
      COMMON / PARFLT /  AS0, RF2, RR2
      COMMON / PARCHR /  SCHEME, ANSATZ

      COMMON / BETABLK / BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX),
     &                   BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)
      COMMON / WGAMMA  /  GAM0, GAM1, GAM2
      COMMON / WC      /  C0, C1, C2


      LRR2 = LOG(RR2)
      LRF2 = LOG(RF2)

      IF ( PROCESS .EQ. 'DVCS') THEN
        S1 = HS1(J + 1.5d0) - HS1(J + 2.0d0) + 2.0d0 * LOG(2.0d0)
     &        - LRF2
        S2 = S1*S1 - HS2(J + 1.5d0) + HS2(J + 2.0d0)
      ELSE
        S1 =  - LRF2
        S2 = S1*S1
      END IF

      CALL VECMAT(C0, GAM0, VM00)
      IF (P .GE. 2) THEN
      CALL VECMAT(C0, GAM1, VM01)
      CALL VECMAT(C1, GAM0, VM10)
      CALL VECMAT(VM00, GAM0, VM000)
      END IF

      IF (SCHEME .EQ. 'EVOLQ') THEN
        BIGC0(1) = (1.0d0,0.0d0)
        BIGC0(2) = (0.0d0,0.0d0)
        DO 8 K = 1, 2
        BIGC1(K) = (0.0d0, 0.0d0)
  8     BIGC2(K) = (0.0d0, 0.0d0)
        RETURN
      ELSE IF (SCHEME .EQ. 'EVOLG') THEN
        BIGC0(1) = (0.0d0,0.0d0)
        BIGC0(2) = (1.0d0,0.0d0)
        DO 9 K = 1, 2
        BIGC1(K) = (0.0d0, 0.0d0)
  9     BIGC2(K) = (0.0d0, 0.0d0)
        RETURN
      END IF

      BIGC0(1) = (1.0d0,0.0d0)
      BIGC0(2) = (0.0d0,0.0d0)

      DO 10 K = 1, 2
      BIGC1(K) = C1(K) + 0.5d0 * S1 * VM00(K)
      IF (P .GE. 2) THEN
      BIGC2(K) = C2(K) + 0.5d0 * S1 * (VM01(K) + VM10(K)) +
     &      0.125d0 * S2 * VM000(K) + 0.5d0 * BETA0(NF) * (
     &      BIGC1(K) * LRR2 + 0.25d0 * VM00(K) * LRF2**2 )
      END IF
 10   CONTINUE


      RETURN
      END
C     *******


C     ****s* cdvcs.f/VECMAT
C  NAME
C     VECMAT  --   multiplies row vector with matrix
C  SYNOPSIS
C     SUBROUTINE VECMAT (VEC, MAT, VM)
C
C     DOUBLE COMPLEX MAT(2,2), VEC(2), VM(2)
C  INPUTS
C         VEC -- row vector
C         MAT -- 2x2 matrix
C  OUTPUT
C          VM -- resulting vector VM = VEC . MAT
C  PARENTS
C     CDVCSF
C  SOURCE
C

      SUBROUTINE VECMAT (VEC, MAT, VM)

      IMPLICIT NONE
      DOUBLE COMPLEX MAT(2,2), VEC(2), VM(2)

      VM(1) = VEC(1) * MAT(1,1) + VEC(2) * MAT(2,1)
      VM(2) = VEC(1) * MAT(1,2) + VEC(2) * MAT(2,2)

      RETURN
      END
C     *******
