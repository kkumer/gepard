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

      SUBROUTINE CDVCSF (  K, BIGCTMP )

      IMPLICIT NONE
      INTEGER  K
      DOUBLE COMPLEX BIGCTMP(0:2,2)
      INTEGER L
      DOUBLE PRECISION LRF2, LRR2
      DOUBLE COMPLEX J, HS1, HS2
      DOUBLE COMPLEX SHIFT1, SHIFT2
      DOUBLE COMPLEX C0(2), C1(2), C2(2)
      DOUBLE COMPLEX VM00(2), VM01(2), VM10(2), VM000(2)
      INCLUDE 'header.f'

      J = N(K) - 1

      LRF2 = LOG(RF2)
      LRR2 = LOG(RR2)

      IF (SCHEME .EQ. 'EVOLQ') THEN
        BIGCTMP(0, 1) = (1.0d0,0.0d0)
        BIGCTMP(0, 2) = (0.0d0,0.0d0)
        DO 10 L = 1, 2
        BIGCTMP(1, L) = (0.0d0, 0.0d0)
 10     BIGCTMP(2, L) = (0.0d0, 0.0d0)
        RETURN
      ELSE IF (SCHEME .EQ. 'EVOLG') THEN
        BIGCTMP(0, 1) = (0.0d0,0.0d0)
        BIGCTMP(0, 2) = (1.0d0,0.0d0)
        DO 20 L = 1, 2
        BIGCTMP(1, L) = (0.0d0, 0.0d0)
 20     BIGCTMP(2, L) = (0.0d0, 0.0d0)
        RETURN
      END IF

      IF ( PROCESS(:3) .EQ. 'DVC' ) THEN
        SHIFT1 = HS1(J + 1.5d0) - HS1(J + 2.0d0) + 2.0d0 * LOG(2.0d0)
     &        - LRF2
        SHIFT2 = SHIFT1*SHIFT1 - HS2(J + 1.5d0) + HS2(J + 2.0d0)
      ELSE
        SHIFT1 =  - LRF2
        SHIFT2 = SHIFT1*SHIFT1
      END IF

*     auxilliary W-coeffs "small C" at given point K
      IF ( PROCESS(:3) .EQ. 'DVC' ) THEN
        DO 30 L = 1, 2
          C0(L) = CDIS1(K, 0, L) 
          C1(L) = CDIS1(K, 1, L) 
          C2(L) = CDIS1(K, 2, L) 
 30     CONTINUE
      ELSE
        DO 40 L = 1, 2
          C0(L) = CDIS2(K, 0, L) 
          C1(L) = CDIS2(K, 1, L) 
          C2(L) = CDIS2(K, 2, L) 
 40     CONTINUE
      END IF

      CALL VECMAT(C0, 0, K, VM00)
      IF (P .GE. 2) THEN
      CALL VECMAT(C0, 1, K, VM01)
      CALL VECMAT(C1, 0, K, VM10)
      CALL VECMAT(VM00, 0, K, VM000)
      END IF


      BIGCTMP(0, 1) = (1.0d0,0.0d0)
      BIGCTMP(0, 2) = (0.0d0,0.0d0)

      DO 50 L = 1, 2
      BIGCTMP(1, L) = C1(L) + 0.5d0 * SHIFT1 * VM00(L)
      IF (P .GE. 2) THEN
      BIGCTMP(2, L) = C2(L) + 0.5d0 * SHIFT1 * (VM01(L) + VM10(L)) +
     &      0.125d0 * SHIFT2 * VM000(L) + 0.5d0 * BETA0(NF) * (
     &      BIGCTMP(1, L) * LRR2 + 0.25d0 * VM00(L) * LRF2**2 )
      END IF
 50   CONTINUE


      RETURN
      END
C     *******


C     ****s* cdvcs.f/VECMAT
C  NAME
C     VECMAT  --   multiplies row vector with anomalous matrix
C  SYNOPSIS
C     SUBROUTINE VECMAT (VEC, ORD, VM)
C
C     DOUBLE COMPLEX VEC(2), VM(2)
C  INPUTS
C         VEC -- row vector
C  OUTPUT
C          VM -- resulting vector VM = VEC . GAMMA^(ORD)_J(K)
C  PARENTS
C     CDVCSF
C  SOURCE
C

      SUBROUTINE VECMAT (VEC, ORD, K, VM)

      IMPLICIT NONE
      INTEGER  K, ORD
      DOUBLE COMPLEX VEC(2), VM(2)
      INCLUDE 'header.f'

      VM(1) = VEC(1) * GAM(K,ORD,1,1) + VEC(2) * GAM(K,ORD,2,1)
      VM(2) = VEC(1) * GAM(K,ORD,1,2) + VEC(2) * GAM(K,ORD,2,2)

      RETURN
      END
C     *******
