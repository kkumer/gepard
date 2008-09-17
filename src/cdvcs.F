C     ****h* gepard/cdvcs.F
C  FILE DESCRIPTION
C    calculation of singlet  Wilson coefficients for DVCS in CSbar scheme
C    and DIS (where MSbar=CSbar)
C
C    $Id$
C     *******


C     ****s* cdvcs.F/CDVCSF
C  NAME
C     CDVCSF  --   singlet Wilson coefficients  "big C"
C  DESCRIPTION
C    calculates Wilson coefficients for DVCS 
C    according to  formulas from [Kumericki:2007sa]
C       (101b-c) - CSbar;   (127,129) - MSbar
C    Also used for DIS where
C    CSbar=MSbar, and for special 'EVOLQ' and 'EVOLG'
C    schemes used for evolution of GPDs/PDFs without
C    multiplying with Wilson coefficients..
C  SYNOPSIS

      SUBROUTINE CDVCSF (  K, BIGCTMP )

      IMPLICIT NONE
      INTEGER  K
      DOUBLE COMPLEX BIGCTMP(0:2,2)

C  INPUTS
C           K -- conformal moment index
C  OUTPUT
C       BIGCTMP(P, flavour) --  "big  C^{(P)}"
C  PARENTS
C      INIT
C  CHILDREN
C      VECMAT, HS1, HS2
C  SOURCE
C

      INTEGER L
      DOUBLE PRECISION LRF2, LRR2
      DOUBLE COMPLEX J, HS1, HS2
      DOUBLE COMPLEX SHIFT1, SHIFT2
      DOUBLE COMPLEX C0(2), C1(2), C2(2)
      DOUBLE COMPLEX VM00(2), VM01(2), VM10(2), VM000(2)
#include "header.F"

      J = N(K) - 1

      LRF2 = LOG(RF2)
      LRR2 = LOG(RR2)


      IF ( PROCESS(:4) .EQ. 'DVCS' ) THEN
        SHIFT1 = HS1(J + 1.5d0) - HS1(J + 2.0d0) + 2.0d0 * LOG(2.0d0)
     &        - LRF2
        SHIFT2 = SHIFT1*SHIFT1 - HS2(J + 1.5d0) + HS2(J + 2.0d0)
      ELSE
        SHIFT1 =  - LRF2
        SHIFT2 = SHIFT1*SHIFT1
      END IF

*     auxilliary W-coeffs "small C" at given point K
      IF ( PROCESS(:4) .EQ. 'DVCS' ) THEN
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


C     ****s* cdvcs.F/VECMAT
C  NAME
C     VECMAT  --   multiplies row vector with anomalous matrix
C  SYNOPSIS

      SUBROUTINE VECMAT (VEC, ORD, K, VM)

      IMPLICIT NONE
      INTEGER  K, ORD
      DOUBLE COMPLEX VEC(2), VM(2)

C  INPUTS
C         VEC(flavour) -- row vector
C                  ORD -- order of perturbation expansion i.e. P
C                    K -- conformal moment index
C  OUTPUT
C          VM(flavour) -- resulting vector VM = VEC . GAMMA^(ORD)_J(K)
C  PARENTS
C     CDVCSF
C  SOURCE
C

#include "header.F"

      VM(1) = VEC(1) * GAM(K,ORD,1,1) + VEC(2) * GAM(K,ORD,2,1)
      VM(2) = VEC(1) * GAM(K,ORD,1,2) + VEC(2) * GAM(K,ORD,2,2)

      RETURN
      END
C     *******
