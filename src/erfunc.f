C     ****h* gepard/erfunc.f
C  FILE DESCRIPTION
C    calculation of scale dependent part of the evolution operator
C
C    $Id$
C     *******


C     ****s* erfunc.f/ERFUNCF
C  NAME
C     ERFUNCF  --   scale dependent part of the evolution operator
C  DESCRIPTION
C    calculates scale dependent part of the evolution operator
C    i.e.  {\cal R}^{(n)}_{ab} - functions
C    |latex \begin{eqnarray*}
C    |latex {\cal R}^{(n)}_{ab}\equiv \beta_0{^{ab}\!R}_j(\mu,\mu_0|n)=
C    |latex \frac{\beta_0}{
C    |latex n \beta_0+{^{a}\! \lambda}_j-{^{b}\! \lambda}_j}\left[
C    |latex 1- \left(\frac{\alpha_s(\mu_0)}{\alpha_s(\mu)}\right)^{\frac{n
C    |latex \beta_0+{^{a}\! \lambda}_j-{^{b}\! \lambda}_j}{\beta_0}}
C    |latex \right]
C    |latex \end{eqnarray*}
C    For convenience, also eigenvalues of LO evolution operator
C    divided by beta_0 are returned.
C  SYNOPSIS
C     SUBROUTINE ERFUNCF (K, R, LAMB, ERFUNC1, ERFUNC2)
C
C     INTEGER K
C     DOUBLE PRECISION R
C     DOUBLE COMPLEX LAMB(2), ERFUNC1(2,2), ERFUNC2(2,2)
C  INPUTS
C           K -- Mellin-Barnes integration point index
C           R -- ratio of astrong(mu)/astrong(mu0)
C  OUTPUT
C        LAMB -- \lambda/\beta_0 where \lambda is eigenvalue
C                of LO evolution operator
C     ERFUNC1 -- matrix {\mathcal R}(mu, mu0 | 1}_{ab}
C     ERFUNC2 -- matrix {\mathcal R}(mu, mu0 | 2}_{ab}
C  PARENTS
C      EVOLF
C  CHILDREN
C      LAMBDAF
C  SOURCE
C

      SUBROUTINE ERFUNCF (K, R, LAMB, ERFUNC1, ERFUNC2)

      IMPLICIT NONE
      INTEGER K
      DOUBLE PRECISION R
      DOUBLE COMPLEX LAMB(2), ERFUNC1(2,2), ERFUNC2(2,2)
      INTEGER SPEED, ACC, P, NF
      INTEGER NFMIN, NFMAX, I, J
      DOUBLE PRECISION RINV, BETA0, BETA1, BETA2, BETA3
      DOUBLE COMPLEX LAM(2)
      PARAMETER (NFMIN = 3, NFMAX = 6)

*   Input common-blocks

      COMMON / PARINT /  SPEED, ACC, P, NF

      COMMON / BETABLK / BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX),
     &                   BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)

      RINV = 1.0d0 / R
      CALL LAMBDAF (K, LAM)

      DO 10 I = 1, 2
 10   LAMB(I) = LAM(I) / BETA0(NF)

      DO 20 I = 1, 2
      DO 20 J = 1, 2
      ERFUNC1(I,J)=( 1.0d0 - RINV**(1.0d0 + LAMB(I) - LAMB(J)) )/
     &      (1.0d0 + LAMB(I) - LAMB(J))
      IF (P. GE. 2) THEN
      ERFUNC2(I,J)=( 1.0d0 - RINV**(2.0d0 + LAMB(I) - LAMB(J)) )/
     &      (2.0d0 + LAMB(I) - LAMB(J))
      END IF
 20   CONTINUE

      RETURN
      END
C     *****
