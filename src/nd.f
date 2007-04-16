C     ****h* gepard/nd.f
C  FILE DESCRIPTION
C     Calculation of contribution of non-diagonal evolution 
C
C    $Id$
C     *******

C     ****s* nd.f/NDNS
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

      SUBROUTINE NDINTF (SEC, K, R, NDINT, NI, NJ)

      IMPLICIT NONE
      INTEGER SEC, K, NI, NJ, ACCND, SPEEDND
      DOUBLE PRECISION R
      DOUBLE COMPLEX NDINT
      INTEGER L, NINTGNDMAX, NPTSND, NPTSNDMAX
      DOUBLE COMPLEX PIHALF, EPH, EPHND
      DOUBLE COMPLEX J, ZK, CB1, CB1C
      DOUBLE COMPLEX DCTAN
      PARAMETER ( NINTGNDMAX = 12, NPTSNDMAX = 768)
      DOUBLE PRECISION DOWNND(NINTGNDMAX+1)
      DOUBLE PRECISION YND(NPTSNDMAX), WGND(NPTSNDMAX)
      DOUBLE COMPLEX ZND(NPTSNDMAX)
      INCLUDE 'header.f'

      PIHALF = (1.5707963267948966d0, 0.0d0)

      EPHND = EXP ( COMPLEX(0.D0, PHIND) )
      EPH = EXP ( COMPLEX(0.D0, PHI) )

      DATA DOWNND 
     & / 0.D0, 0.01D0, 0.025D0, 0.067D0, 0.18D0, 
     &         0.5D0, 1.3D0, 3.7D0, 10.D0, 
     &         25.D0, 0.67D2, 1.8D2, 5.0D2 /

      ACCND = ACC
      SPEEDND = SPEED

      CALL INTEGRAF(ACCND, SPEEDND, DOWNND, NINTGNDMAX, YND, WGND)

      NPTSND = 2**ACCND * NINTGNDMAX / SPEEDND

*   Calculating actual Mellin-Barnes contour points from Gaussian abscissae

      EPHND = EXP ( COMPLEX(0.D0, PHIND) )
      DO 10 L = 1, NPTSND
 10     ZND(L) = CND + YND(L) * EPHND 

*   Integrating over L->ZK, with fixed K->J

      J = N(SEC,K) - 1.0d0
      NDINT = (0.0d0, 0.0d0)
      DO 100 L = 1, NPTSND

*       Mellin-Barnes contour for ND evolution intersects real axis at CND
*       so contour saved in common blocks needs to be shifted by CND-C.
*       Note that N(L)-1-C = Y*EPH
        ZK = ZND(L)

        CALL CB1F (J + ZK + 2.0d0, J, R, CB1, NI, NJ)
        CALL CB1F (J + CONJG(ZK) + 2.0d0, J, R, CB1C, NI, NJ)


        NDINT = NDINT + WGND(L)*( EPHND * CB1 / DCTAN(PIHALF * ZK)
     &      - CONJG(EPHND) * CB1C / DCTAN(PIHALF * CONJG(ZK)) )

 100  CONTINUE

      NDINT = (0.0d0, 0.25d0) * NDINT

      RETURN
      END
C     ****


C     ****s* nd.f/CB1F
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

      SUBROUTINE CB1F (ZN, ZK, R, CB1, NI, NJ)

      IMPLICIT NONE
      INTEGER NI, NJ
      DOUBLE COMPLEX ZN, ZK, CB1
      DOUBLE PRECISION R
      INTEGER A, B
      DOUBLE PRECISION NFD
      DOUBLE COMPLEX S1ORIG
      DOUBLE COMPLEX R1, AAA
      DOUBLE COMPLEX GOD(2,2), DM(2,2)
      DOUBLE COMPLEX PROJGOD(2,2,2,2), PROJDM(2,2,2,2)
      DOUBLE COMPLEX HS1, POCHHAMMER
      DOUBLE COMPLEX GAUX, GAMN(2,2), GAMK(2,2)
      DOUBLE COMPLEX LAMN(2), LAMK(2)
      DOUBLE COMPLEX PN(2,2,2), PK(2,2,2)
      DOUBLE COMPLEX ERFUNC1(2,2), ERFUNC2(2,2)
      INCLUDE 'header.f'

      NFD = DBLE(NF)

*   S1 needs to be reinitialized at correct point ... 
      S1ORIG = S1
      S1 = HS1(ZN+1.0d0)
      CALL WgammaVQQ0F(NFD, ZN+1.0d0 , GAUX)
      GAMN(1, 1) = GAUX
      CALL WgammaVQG0F(NFD, ZN+1.0d0 , GAUX)
      GAMN(1, 2) = GAUX
      CALL WgammaVGQ0F(NFD, ZN+1.0d0 , GAUX)
      GAMN(2, 1) = GAUX
      CALL WgammaVGG0F(NFD, ZN+1.0d0 , GAUX)
      GAMN(2, 2) = GAUX
      S1 = HS1(ZK+1.0d0)
      CALL WgammaVQQ0F(NFD, ZK+1.0d0 , GAUX)
      GAMK(1, 1) = GAUX
      CALL WgammaVQG0F(NFD, ZK+1.0d0 , GAUX)
      GAMK(1, 2) = GAUX
      CALL WgammaVGQ0F(NFD, ZK+1.0d0 , GAUX)
      GAMK(2, 1) = GAUX
      CALL WgammaVGG0F(NFD, ZK+1.0d0 , GAUX)
      GAMK(2, 2) = GAUX
*   ... and then returned to original value
      S1 =  S1ORIG

*     First calculate eigenvalues, projectors and R-function

      CALL LAMBDANDF(GAMN, GAMK, LAMN, LAMK)
      CALL PROJECTORSNDF(GAMN, GAMK, LAMN, LAMK, PN, PK)
      CALL ERFUNCF(R, LAMN, LAMK, ERFUNC1, ERFUNC2)

*   A_nk
      AAA = HS1((ZN+ZK+2.0d0)/2.0d0) - HS1((ZN-ZK-2.0d0)/2.0d0)
     &      + 2.0d0*HS1(ZN-ZK-1.0d0) - HS1(ZN+1.0d0)

*   GOD = "G Over D" = g_jk / d_jk

      GOD(1,1) = 2.0d0 * CF * ( 2.0d0*AAA +(AAA - HS1(ZN+1.0d0)
     &      )*(ZN-ZK)*(ZN+ZK+3.0d0)/(ZK+1.0d0)/(ZK+2.0d0) )

      GOD(1,2) = 0.0d0

      GOD(2,1) = 2.0d0*CF*(ZN-ZK)*(ZN+ZK+3.0d0)/ZN/(ZK+1.0d0)/(ZK+2.0d0)

      GOD(2,2) = 2.0d0 * CA * ( 2.0d0*AAA +(AAA - HS1(ZN+1.0d0)
     &  ) * ( POCHHAMMER(ZN,4) / POCHHAMMER(ZK,4) - 1.0d0 )
     &  + 2.0d0 * (ZN-ZK) * (ZN+ZK+3.0d0) / POCHHAMMER(ZK,4) ) * ZK / ZN

*   Diagonal matrix multiplying d_jk:

      DM(1,1) = 1.0d0
      DM(1,2) = 0.0d0
      DM(2,1) = 0.0d0
      DM(2,2) = ZK / ZN

*   Projecting P.DM.P  and P.GOD.P

      CALL PROJECTION(PN, DM, PK, PROJDM)
      CALL PROJECTION(PN, GOD, PK, PROJGOD)


      IF ( ANSATZ(:2) .EQ. 'NS' ) THEN
*       non-singlet

        R1 = ( 1.0d0 - (1.0d0/R)**((BETA0(NF) + GAMN(1,1) 
     &    - GAMK(1,1))/BETA0(NF)) ) / ( BETA0(NF)+GAMN(1,1)-GAMK(1,1) )
*      \tilde{B}_nk
        CB1 = R1 * (GAMN(1,1)-GAMK(1,1)) * (BETA0(NF) 
     &    - GAMK(1,1) + GOD(1,1)) * R**(-GAMK(1,1)/BETA0(NF))

      ELSE
*       singlet

*       Doing the summation over A and B to get \tilde{B}_nk
      CB1 = (0.0d0, 0.0d0)
      DO 20 B = 1, 2
        DO 10 A = 1, 2
          CB1 = CB1 + ERFUNC1(A,B) * (LAMN(A) - LAMK(B)) * (
     &         (BETA0(NF) - LAMK(B)) * PROJDM(A,B,NI,NJ) 
     &         + PROJGOD(A,B,NI,NJ) ) / BETA0(NF)
 10     CONTINUE
        CB1 = CB1 * R**( - LAMK(B) / BETA0(NF) )
 20   CONTINUE

      ENDIF


*   together with prefactors ...
      CB1 = CB1 * (ZK+1.0d0)*(ZK+2.0d0)*(2.0d0*ZN+3.0D0)
     &           /(ZN+1.0d0)/(ZN+2.0d0)/(ZN-ZK)/(ZN+ZK+3.0d0)
      RETURN
      END
C     ****
