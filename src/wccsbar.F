C     ****h* gepard/wccsbar.F
C  FILE DESCRIPTION
C    calculation of Wilson coefficients for DVCS in CSbar scheme
C    and for DIS (where MSbar=CSbar)
C
C    $Id$
C     *******


C     ****s* wccsbar.F/WCCSBAR
C  NAME
C     WCCSBAR  --   singlet Wilson coefficients  "big C"
C  DESCRIPTION
C    calculates Wilson coefficients for DVCS 
C    according to  formulas from [Kumericki:2007sa]
C       (101b-c) - CSbar;   (127,129) - MSbar
C    Also used for DIS where CSbar=MSbar, 
C  SYNOPSIS

      SUBROUTINE WCCSBAR (PIDK, SEC, K, WCTMP )

      IMPLICIT NONE
      INTEGER  PIDK, SEC, K
      DOUBLE COMPLEX WCTMP(0:2, 4)

C  INPUTS
C        PIDK -- process ID (0 for DIS, 1 for DVCS)
C         SEC -- SO(3) PW index
C           K -- conformal moment index
C  OUTPUT
C       WCTMP(P, flavour) --  "big  C^{(P)}"
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

      J = N(SEC, K) - 1

      LRF2 = LOG(RF2)
      LRR2 = LOG(RR2)


      IF ( PIDK .EQ. 1 ) THEN
*       DVCS (CSBAR)
        SHIFT1 = HS1(J + 1.5d0) - HS1(J + 2.0d0) + 2.0d0 * LOG(2.0d0)
     &        - LRF2
        SHIFT2 = SHIFT1*SHIFT1 - HS2(J + 1.5d0) + HS2(J + 2.0d0)
      ELSE
*       DIS
        SHIFT1 =  - LRF2
        SHIFT2 = SHIFT1*SHIFT1
      END IF

*     auxilliary W-coeffs "small C" at given point K

      IF ( PIDK .EQ. 1 ) THEN
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

      CALL VECMAT(C0, 0, SEC, K, VM00)
      CALL VECMAT(C0, 1, SEC, K, VM01)
      CALL VECMAT(C1, 0, SEC, K, VM10)
      CALL VECMAT(VM00, 0, SEC, K, VM000)


*     --- LO ---
      WCTMP(0, 1) = (1.0d0, 0.0d0)
      WCTMP(0, 2) = (0.0d0, 0.0d0)
      WCTMP(0, 3) = (1.0d0, 0.0d0)
      WCTMP(0, 4) = (1.0d0, 0.0d0)

*     --- NLO and NNLO ---
      DO 50 L = 1, 2
        WCTMP(1, L) = C1(L) + 0.5d0 * SHIFT1 * VM00(L)
        WCTMP(2, L) = C2(L) + 0.5d0 * SHIFT1 * (VM01(L) + VM10(L)) +
     &        0.125d0 * SHIFT2 * VM000(L) + 0.5d0 * BETA0(NF) * (
     &        WCTMP(1, L) * LRR2 + 0.25d0 * VM00(L) * LRF2**2 )
 50   CONTINUE
*     TODO: implement NS- i.e. NSM
      WCTMP(1, 3) = CDISNS1(K, 1) + 0.5d0 * SHIFT1 * CDISNS1(K, 0) *
     &                GAMNS(SEC,K,0,0)
      WCTMP(2, 3) = CDISNS1(K, 2) + 0.5d0 * SHIFT1 * (CDISNS1(K, 0)
     &       * GAMNS(SEC,K,1,0) + CDISNS1(K, 1) * GAMNS(SEC,K,0,0)) +
     &      0.125d0 * SHIFT2 * CDISNS1(K, 0) * GAMNS(SEC,K,0,0)**2 + 
     &      0.5d0 * BETA0(NF) * ( WCTMP(1, 3) * LRR2 + 0.25d0 * 
     &        CDISNS1(K, 0) * GAMNS(SEC,K,0,0) * LRF2**2 )


      RETURN
      END
C     *******


C     ****s* wccsbar.F/VECMAT
C  NAME
C     VECMAT  --   multiplies row vector with anomalous matrix
C  SYNOPSIS

      SUBROUTINE VECMAT (VEC, ORD, SEC, K, VM)

      IMPLICIT NONE
      INTEGER  K, ORD, SEC
      DOUBLE COMPLEX VEC(2), VM(2)

C  INPUTS
C         VEC(flavour) -- row vector
C                  ORD -- order of perturbation expansion i.e. P
C                  SEC -- SO(3) PW
C                    K -- conformal moment index
C  OUTPUT
C          VM(flavour) -- resulting vector VM = VEC . GAMMA^(ORD)_J(K)
C  PARENTS
C     WCCSBAR
C  SOURCE
C

#include "header.F"

      VM(1) = VEC(1) * GAM(SEC,K,ORD,1,1) + VEC(2) * GAM(SEC,K,ORD,2,1)
      VM(2) = VEC(1) * GAM(SEC,K,ORD,1,2) + VEC(2) * GAM(SEC,K,ORD,2,2)

      RETURN
      END
C     *******
