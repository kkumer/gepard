C     ****h* gepard/nd.f
C  FILE DESCRIPTION
C     Calculation of contribution of non-diagonal evolution 
C
C    $Id: rnnlo.f 33 2006-07-26 18:42:37Z kuk05260 $
C     *******

C     ****s* rnnlo.f/NDNS
C  NAME
C        NDNS  --   R_1
C  DESCRIPTION
C     Calculation non-singlet R_1  (cf. my DIS-p61  for singlet case)
C     - combination of gamma^(n) and beta_m needed in evolution operator 
C  SYNOPSIS
C     SUBROUTINE NDNS (K, R, CNDNS)
C
C     INTEGER K
C     DOUBLE PRECISION R
C     DOUBLE COMPLEX CNDNS
C  INPUTS
C           K -- Mellin-Barnes integration point index
C           R -- ratio of astrong(mu)/astrong(mu0)
C  OUTPUT
C       CNDNS -- non-diagonally evolved Wilson coefficient
C  PARENTS
C      PARWAV
C  CHILDREN
C      CB1F, DCTAN
C  SOURCE
C

      SUBROUTINE NDNSF (K, R, CNDNS)

      IMPLICIT NONE
      INTEGER K
      DOUBLE PRECISION R
      DOUBLE COMPLEX CNDNS
      INTEGER L
      DOUBLE COMPLEX PIHALF, EPH, EPHND
      DOUBLE COMPLEX J, ZK, CB1, CB1C
      DOUBLE COMPLEX DCTAN
      INCLUDE 'header.f'

      PIHALF = (1.5707963267948966d0, 0.0d0)

      EPHND = EXP ( COMPLEX(0.D0, PHIND) )
      EPH = EXP ( COMPLEX(0.D0, PHI) )


*   Integrating over L->ZK, with fixed K->J

      J = N(K) - 1.0d0
      CNDNS = (0.0d0, 0.0d0)
      DO 100 L = 1, NPTS
*   Mellin-Barnes contour for ND evolution intersects real axis at CND
*   so contour saved in common blocks needs to be shifted by CND-C.
*   Note that N(L)-1-C = Y*EPH
        ZK = (N(L) - 1.0d0 - C)*CONJG(EPH)*EPHND + CND
        CALL CB1F (J + ZK + 2.0d0, J, R, CB1)
!        CNDNS = CNDNS + WG(L)*IMAGPART(EPH*CB1/DCTAN(PIHALF*ZK))
        CALL CB1F (J + CONJG(ZK) + 2.0d0, J, R, CB1C)
        CNDNS = CNDNS + WG(L)*( EPHND*CB1/DCTAN(PIHALF*ZK)
     &           -CONJG(EPHND)*CB1C/DCTAN(PIHALF*CONJG(ZK)) )
 100  CONTINUE

!      CNDNS = -0.5 * CNDNS
      CNDNS = (0.0d0, 0.25d0) * CNDNS

      RETURN
      END
C     ****


C     ****s* rnnlo.f/CB1
C  NAME
C        CB1  --  
C  DESCRIPTION
C     NLO non-diagonally evolved non-singlet LO Wilson coefficient.
C     It's multiplied by GAMMA(3/2) GAMMA(K+3) / (2^(K+1) GAMMA(K+5/2))
C     so it's ready to be combined with diagonally evolved C_K = 1 + ...
C     where in the end everything will be multiplied by 
C     (2^(K+1) GAMMA(K+5/2))  / ( GAMMA(3/2) GAMMA(K+3) )
C  SYNOPSIS
C     SUBROUTINE CB1F (ZN, ZK, R,  CB1)
C
C     DOUBLE COMPLEX ZN, ZK, CB1
C     DOUBLE PRECISION R
C  INPUTS
C          ZN -- non-diagonal evolution Mellin-Barnes integration point
C          ZK -- COPE Mellin-Barnes integration point
C           R -- ratio of astrong(mu)/astrong(mu0)
C  OUTPUT
C         CB1 -- Cnormalization(K)^-1 C_N  B^{(1)}_{N, K}
C  PARENTS
C      NDNS
C  CHILDREN
C      HS1, WgammaVNSP0F
C  SOURCE
C

      SUBROUTINE CB1F (ZN, ZK, R, CB1)

      IMPLICIT NONE
      DOUBLE COMPLEX ZN, ZK, CB1
      DOUBLE PRECISION R
      DOUBLE PRECISION NFD
      DOUBLE COMPLEX S1ORIG
      DOUBLE COMPLEX R1, A, GOD
      DOUBLE COMPLEX HS1, GAMN, GAMK
      INCLUDE 'header.f'

      NFD = DBLE(NF)

*   S1 needs to be reinitialized and then returned to original value
      S1ORIG = S1
      S1 = HS1(ZN+1.0d0)
      CALL WgammaVNSP0F(NFD, ZN+1.0d0 , GAMN)
      S1 = HS1(ZK+1.0d0)
      CALL WgammaVNSP0F(NFD, ZK+1.0d0 , GAMK)
      S1 =  S1ORIG

      R1 = ( 1.0d0 - (1.0d0/R)**( (BETA0(NF)+GAMN-GAMK)/BETA0(NF)) )
     &       / ( BETA0(NF)+GAMN-GAMK )

*   A_nk
      A = HS1((ZN+ZK+2.0d0)/2.0d0) - HS1((ZN-ZK-2.0d0)/2.0d0)
     &      + 2.0d0*HS1(ZN-ZK-1.0d0) - HS1(ZN+1.0d0)

*   GOD = "G Over D" = g_jk / d_jk
      GOD = (8.0d0/3.0d0) * ( 2.0d0*A +(A - HS1(ZN+1.0d0)
     &      )*(ZN-ZK)*(ZN+ZK+3.0d0)/(ZK+1.0d0)/(ZK+2.0d0) )

*   \tilde{B}_nk
      CB1 = R1 * (GAMN-GAMK) * (BETA0(NF) - GAMK + GOD)

*   together with prefactors ...
      CB1 = CB1 * (ZK+1.0d0)*(ZK+2.0d0)*(2.0d0*ZN+3.0D0)
     &           /(ZN+1.0d0)/(ZN+2.0d0)/(ZN-ZK)/(ZN+ZK+3.0d0)
      RETURN
      END
C     ****
