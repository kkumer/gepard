C     ****h* gepard/rnnlo.f
C  FILE DESCRIPTION  - THIS FILE IS NOT USED
C     Calculation of matrices R_1 and R_2, from my DIS-p61 
C     - combination of gamma^(n) and beta_m needed in evolution operator 
C     - projected on +/- directions
C
C    $Id$
C     *******

C     ****s* rnnlo.f/RNNLONSF
C  NAME
C     RNNLONSF  --   R_1
C  DESCRIPTION
C     Calculation non-singlet R_1  (cf. my DIS-p61  for singlet case)
C     - combination of gamma^(n) and beta_m needed in evolution operator 
C  SYNOPSIS
C     SUBROUTINE RNNLONSF (K, R1)
C
C     INTEGER K
C     DOUBLE COMPLEX R1
C  INPUTS
C           K -- Mellin-Barnes integration point index
C  OUTPUT
C          R1
C  IDENTIFIERS
C       BETABLK, NGAM -- common blocks with beta function coefficients
C       of QCD, and moments of anomalous dimensions of DIS
C  PARENTS
C      EVOLNSF
C  SOURCE
C

      SUBROUTINE RNNLONSF (K, R1, R2)

      IMPLICIT NONE
      INTEGER K
      DOUBLE COMPLEX R1, R2
      INTEGER K1, L
      DOUBLE PRECISION INV
      INCLUDE 'header.f'


*     Inverse beta_0 is often needed below
      INV = 1.0d0 / BETA0(NF)

      R1 = INV * (GAMNS(K,1) 
     &            - 0.5d0 * INV * BETA1(NF) * GAMNS(K,0))

      IF (P. GE. 2) THEN

        R2 = INV * ( GAMNS(K, 2) - 0.5d0 * INV * BETA1(NF) * 
     &     GAMNS(K, 1) - 0.25d0 * INV * BETA2(NF) * GAMNS(K, 0) +
     &     0.25d0 * INV**2 * BETA1(NF)**2 * GAMNS(K, 0) )

      END IF

      RETURN
      END
C     ****


C     ****s* rnnlo.f/RNNLOF
C  NAME
C     RNNLOF  --   Matrices R_1 and R_2
C  DESCRIPTION
C     Calculation of matrices R_1 and R_2, from my DIS-p61 
C     - combination of gamma^(n) and beta_m needed in evolution operator 
C     - projected on +/- directions
C  SYNOPSIS
C     SUBROUTINE RNNLOF (K, PR, R1PROJ, R2PROJ)
C
C     INTEGER K
C     DOUBLE COMPLEX PR(2,2,2), R1PROJ(2,2,2,2), R2PROJ(2,2,2,2)
C  INPUTS
C           K -- Mellin-Barnes integration point index
C         LAM -- eigenvalues of LO evolution operator
C  OUTPUT
C          PR -- Two projector matrices PR(1,a,b)=P^+
C                and PR(2,a,b)=P^-
C      R1PROJ -- R1PROJ(a,b,i,j) = P_a . R_1(i,j) . P_b
C                a,b \in {+,-};  i,j \in {Q, G}
C      R2PROJ -- R2PROJ(a,b,i,j) = P_a . R_2(i,j) . P_b
C                a,b \in {+,-};  i,j \in {Q, G}
C  IDENTIFIERS
C       BETABLK, NGAM -- common blocks with beta function coefficients
C       of QCD, and moments of anomalous dimensions of DIS
C  PARENTS
C      EVOLF
C  CHILDREN
C      PROJECTORSF, PROJECTION
C  SOURCE
C

      SUBROUTINE RNNLOF (SEC, K, LAM, PR, R1PROJ, R2PROJ)

      IMPLICIT NONE
      INTEGER SEC, K
      DOUBLE COMPLEX LAM(2)
      DOUBLE COMPLEX PR(2,2,2), R1PROJ(2,2,2,2), R2PROJ(2,2,2,2)
      INTEGER K1, L
      DOUBLE PRECISION INV
      DOUBLE COMPLEX R1(2,2), R2(2,2)
      INCLUDE 'header.f'

      CALL PROJECTORSF(SEC, K, LAM, PR)

*     Inverse beta_0 is often needed below
      INV = 1.0d0 / BETA0(NF)

      DO 10 K1 = 1,2
      DO 10 L = 1,2
      R1(K1,L) = INV * (GAM(SEC,K,1,K1,L) 
     &            - 0.5d0 * INV * BETA1(NF) * GAM(SEC,K,0,K1,L))
      IF (P. GE. 2) THEN
        R2(K1,L) = INV * (GAM(SEC,K,2,K1,L) 
     &            - 0.5d0 * INV * BETA1(NF) * GAM(SEC,K,1,K1,L) 
     &            - 0.25d0 * INV * (BETA2(NF) 
     &            - INV * BETA1(NF)**2) * GAM(SEC,K,0,K1,L))
      END IF
 10   CONTINUE

      CALL PROJECTION(PR, R1, PR, R1PROJ)
      IF (P. GE. 2) THEN
      CALL PROJECTION(PR, R2, PR, R2PROJ)
      END IF

      RETURN
      END
C     ****

C     ****s* rnnlo.f/PROJECTION
C  NAME
C     PROJECTION  --  matrix multiplication  PL . M . PR
C  DESCRIPTION
C     Auxilliary subroutine
C  SYNOPSIS
C     SUBROUTINE PROJECTION (PL, M, PR, PROJ)
C
C     DOUBLE COMPLEX PL(2,2,2), M(2,2), PR(2,2,2)
C     DOUBLE COMPLEX PROJ(2,2,2,2)
C  INPUTS
C          PR -- projection matrices P(1)=P^+, P(2)=P^-
C           M -- matrix to be projected
C  OUTPUT
C        PROJ -- 4 possible projections Pa . M . Pb
C                a,b \in {+,-}
C  PARENTS
C     RNNLOF
C  SOURCE
C

      SUBROUTINE PROJECTION (PL, M, PR, PROJ)

      IMPLICIT NONE
      DOUBLE COMPLEX PL(2,2,2), M(2,2), PR(2,2,2)
      DOUBLE COMPLEX PROJ(2,2,2,2)
      INTEGER A, B, C, D, L, R
      DOUBLE COMPLEX AUX(2,2)

      DO 35 L = 1,2
      DO 35 R = 1,2

      DO 15 A = 1,2
      DO 15 C = 1,2
      AUX(A,C) = (0.0d0, 0.0d0)
      DO 15 B = 1,2
 15   AUX(A,C) = AUX(A,C) + PL(L,A,B) * M(B,C) 

      DO 25 A = 1,2
      DO 25 D = 1,2
      PROJ(L,R,A,D) = (0.0d0, 0.0d0)
      DO 25 C = 1,2
 25   PROJ(L,R,A,D) = PROJ(L,R,A,D) + AUX(A,C) * PR(R,C,D)
 35   CONTINUE

      RETURN
      END
C     ****
