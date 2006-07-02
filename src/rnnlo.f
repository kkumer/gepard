C     ****h* gepard/rnnlo.f
C  FILE DESCRIPTION
C     Calculation of matrices R_1 and R_2, from my DIS-p61 
C     - combination of gamma^(n) and beta_m needed in evolution operator 
C     - projected on +/- directions
C
C    $Id$
C     *******


C     ****s* rnnlo.f/RNNLOF
C  NAME
C     RNNLOF  --   Matrices R_1 and R_2
C  DESCRIPTION
C     Calculation of matrices R_1 and R_2, from my DIS-p61 
C     - combination of gamma^(n) and beta_m needed in evolution operator 
C     - projected on +/- directions
C  SYNOPSIS
C     SUBROUTINE RNNLOF (NF, PR, R1PROJ, R2PROJ)
C
C     INTEGER NF
C     DOUBLE COMPLEX PR(2,2,2), R1PROJ(2,2,2,2), R2PROJ(2,2,2,2)
C  INPUTS
C          NF -- number of active flavours
C  OUTPUT
C           P -- approximation order, which is N^{P}LO
C          PR -- Two projector matrices PR(1,a,b)=P^+
C                and PR(2,a,b)=P^-
C      R1PROJ -- R1PROJ(a,b,i,j) = P_a . R_1(i,j) . P_b
C                a,b \in {+,-};  i,j \in {Q, G}
C      R2PROJ -- R2PROJ(a,b,i,j) = P_a . R_2(i,j) . P_b
C                a,b \in {+,-};  i,j \in {Q, G}
C  IDENTIFIERS
C       BETABLK, WGAMMA -- common blocks with beta function coefficients
C       of QCD, and moments of anomalous dimensions of DIS
C  PARENTS
C      EVOLF
C  CHILDREN
C      PROJECTORSF, PROJECTION
C  NOTES
C      One needs to call BETAF and COMMONF first, to initialize
C      common blocks BETABLK and WGAMMA
C  SOURCE
C

      SUBROUTINE RNNLOF (NF, PR, R1PROJ, R2PROJ)

      IMPLICIT NONE
      INTEGER NF
      DOUBLE COMPLEX PR(2,2,2), R1PROJ(2,2,2,2), R2PROJ(2,2,2,2)
      INTEGER NFMIN, NFMAX, K, J, P
      DOUBLE PRECISION BETA0, BETA1, BETA2, BETA3, INV
      DOUBLE COMPLEX GAM0(2,2), GAM1(2,2), GAM2(2,2)
      DOUBLE COMPLEX R1(2,2), R2(2,2)
      PARAMETER (NFMIN = 3, NFMAX = 6)

*   Input common-blocks

      COMMON / BETABLK / BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX),
     &                   BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)
      COMMON / WGAMMA  / GAM0, GAM1, GAM2
      COMMON / APPROX     /  P

      CALL PROJECTORSF(NF, PR)

*     Inverse beta_0 is often needed below
      INV = 1.0d0 / BETA0(NF)

      DO 10 K = 1,2
      DO 10 J = 1,2
      R1(K,J) = INV * (GAM1(K,J) - 0.5d0 * INV * BETA1(NF) * GAM0(K,J))
      IF (P. GE. 2) THEN
      R2(K,J) = INV * (GAM2(K,J) - 0.5d0 * INV * BETA1(NF) * GAM1(K,J) 
     &   - 0.25d0 * INV * (BETA2(NF) - INV * BETA1(NF)**2) * GAM0(K,J))
      END IF
 10   CONTINUE

      CALL PROJECTION(PR, R1, R1PROJ)
      IF (P. GE. 2) THEN
      CALL PROJECTION(PR, R2, R2PROJ)
      END IF

      RETURN
      END
C     ****

C     ****s* rnnlo.f/PROJECTION
C  NAME
C     PROJECTION  --  matrix multiplication  PR . M . PR
C  DESCRIPTION
C     Auxilliary subroutine
C  SYNOPSIS
C     SUBROUTINE PROJECTION (PR, M, PROJ)
C
C     DOUBLE COMPLEX PR(2,2,2), M(2,2)
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

      SUBROUTINE PROJECTION (PR, M, PROJ)

      IMPLICIT NONE
      DOUBLE COMPLEX PR(2,2,2), M(2,2)
      DOUBLE COMPLEX PROJ(2,2,2,2)
      INTEGER A, B, C, D, L, R
      DOUBLE COMPLEX AUX(2,2)

      DO 35 L = 1,2
      DO 35 R = 1,2

      DO 15 A = 1,2
      DO 15 C = 1,2
      AUX(A,C) = (0.0d0, 0.0d0)
      DO 15 B = 1,2
 15   AUX(A,C) = AUX(A,C) + PR(L,A,B) * M(B,C) 

      DO 25 A = 1,2
      DO 25 D = 1,2
      PROJ(L,R,A,D) = (0.0d0, 0.0d0)
      DO 25 C = 1,2
 25   PROJ(L,R,A,D) = PROJ(L,R,A,D) + AUX(A,C) * PR(R,C,D)
 35   CONTINUE

      RETURN
      END
C     ****
