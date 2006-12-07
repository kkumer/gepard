C     ****h* gepard/erfunc.f
C  FILE DESCRIPTION
C    calculation of scale dependent part of the evolution operator
C
C    $Id$
C     *******

C     ****s* erfunc.f/ERFUNCNSF
C  NAME
C     ERFUNCNSF  --   scale dependent part of the non-singlet evolution operator
C  DESCRIPTION
C    calculates scale dependent part of the non-singlet evolution operator
C    i.e.  {\cal R}^{(n)} - functions
C    For convenience, also \gamma_{NS}^{(0)}
C    divided by beta_0 is returned.
C  SYNOPSIS
C     SUBROUTINE ERFUNCF (K, R, GAMB, ERFUNCNS1)
C
C     INTEGER K
C     DOUBLE PRECISION R
C     DOUBLE COMPLEX GAMB, ERFUNCNS1
C  INPUTS
C           K -- Mellin-Barnes integration point index
C           R -- ratio of astrong(mu)/astrong(mu0)
C  OUTPUT
C        GAMB -- \gamma/\beta_0
C     ERFUNCNS1 --  {\mathcal R}(mu, mu0 | 1}
C  PARENTS
C      EVOLNSF
C  SOURCE
C

      SUBROUTINE ERFUNCNSF (K, R, GAMB, ERFUNCNS1)

      IMPLICIT NONE
      INTEGER K
      DOUBLE PRECISION R
      DOUBLE COMPLEX GAMB, ERFUNCNS1
      DOUBLE PRECISION RINV
      INCLUDE 'header.f'

      RINV = 1.0d0 / R

      GAMB = NGAMNS(K,0) / BETA0(NF)

      ERFUNCNS1 = ( 1.0d0 - RINV )

      RETURN
      END
C     *****

C     ****s* erfunc.f/ERFUNCF
C  NAME
C     ERFUNCF  --   scale dependent part of the singlet evolution operator
C  DESCRIPTION
C    calculates scale dependent part of the singlet evolution operator
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
      INTEGER I, J
      DOUBLE PRECISION RINV
      DOUBLE COMPLEX LAM(2)
      INCLUDE 'header.f'


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
