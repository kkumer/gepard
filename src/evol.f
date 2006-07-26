C     ****h* gepard/evol.f
C  FILE DESCRIPTION
C     calculates evolution operator
C
C    $Id$
C     *******


C     ****s* evol.f/EVOLF
C  NAME
C     EVOLF  --   evolution operator
C  DESCRIPTION
C    calculates evolution operator
C    according to KMKPS06 paper
C  SYNOPSIS
C     SUBROUTINE EVOLF (K, R, EVOLA)
C
C     INTEGER K
C     DOUBLE PRECISION R
C     DOUBLE COMPLEX EVOLA(3,2,2), 
C  INPUTS
C           K -- Mellin-Barnes integration point index
C           R -- ratio of astrong(mu)/astrong(mu0)
C  OUTPUT
C       EVOLA -- three 2x2 matrices representing three coefficients 
C                 {\cal A}^{0,1,2}(a_s/a0_s)^{-\gamma^0/beta_0} in
C                 alpha_s/(2\pi) expansion of evolution operator 
C                 \Epsilon(mu, mu0)
C  PARENTS
C      PARWAVF
C  CHILDREN
C      ERFUNCF, RNNLOF, KRONECKER
C  SOURCE
C

      SUBROUTINE EVOLF (K, R, EVOLA)

      IMPLICIT NONE
      INTEGER K
      DOUBLE PRECISION R
      DOUBLE COMPLEX EVOLA(3,2,2)
      INTEGER SPEED, ACC, P, NF
      INTEGER A, B, C, I, J, K1, ORD
      DOUBLE PRECISION RINV
      DOUBLE COMPLEX LAMB(2), PR(2,2,2)
      DOUBLE COMPLEX ERFUNC1(2,2), ERFUNC2(2,2)
      DOUBLE COMPLEX R1PROJ(2,2,2,2), R2PROJ(2,2,2,2)
      DOUBLE COMPLEX AUX, DINV, IJAUX(3), R1R1(2,2)
      INTEGER KRONECKER

*   Input common-blocks 

      COMMON / PARINT /  SPEED, ACC, P, NF


      CALL ERFUNCF (K, R, LAMB, ERFUNC1, ERFUNC2)
      CALL RNNLOF (K, PR, R1PROJ, R2PROJ)

      RINV = 1.0d0 / R

      DO 10 I = 1, 2
      DO 10 J = 1, 2
      DO 10 ORD = 1, P+1
 10      EVOLA(ORD,I,J) = (0.0d0, 0.0d0)

      DO 40 I = 1, 2
      DO 40 J = 1, 2

      DO 40 B = 1, 2
      DO 15 ORD = 1, P+1
 15      IJAUX(ORD) = (0.0d0, 0.0d0)
      DO 30 A = 1, 2
         IJAUX(1) = IJAUX(1) + KRONECKER(A,B) * PR(A,I,J)
         IJAUX(2) = IJAUX(2) - ERFUNC1(A,B) * R1PROJ(A,B,I,J)
      IF (P. GE. 2) THEN
         AUX = (0.0d0, 0.0d0)
         DO 20 C = 1, 2
            DINV = 1.0d0 + LAMB(C) - LAMB(B)
            R1R1(I,J) = (0.0d0, 0.0d0)
            DO 17 K1 = 1, 2
 17         R1R1(I,J) = R1R1(I,J) + R1PROJ(A,C,I,K1) * R1PROJ(C,B,K1,J)
 20         AUX = AUX + (ERFUNC2(A,B) - ERFUNC1(A,C) * RINV**DINV)/DINV
     &            * R1R1(I,J)
         IJAUX(3) = IJAUX(3) + AUX - ERFUNC2(A,B) * R2PROJ(A,B,I,J)
      END IF
 30   CONTINUE
      DO 35 ORD = 1, P+1
 35   EVOLA(ORD,I,J) = EVOLA(ORD,I,J) + IJAUX(ORD) * R**(-LAMB(B))
 40   CONTINUE

      RETURN
      END
C     *****


C     ****f* evol.f/KRONECKER
C  NAME
C    KRONECKER -- kronecker delta symbol
C  SYNOPSIS
C    INTEGER FUNCTION KRONECKER (A, B)
C
C    INTEGER A, B
C  PARENTS
C       EVOLF
C  SOURCE
C

      INTEGER FUNCTION KRONECKER (A, B)

      IMPLICIT NONE
      INTEGER A, B

      IF (A .EQ. B) THEN
            KRONECKER = 1
      ELSE
            KRONECKER = 0
      END IF

      RETURN
      END
C     *****
