C     ****h* gepard/erfunc.f
C  FILE DESCRIPTION
C    calculation of scale dependent part of the evolution operator
C
C    $Id$
C     *******

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
C  SYNOPSIS
C     SUBROUTINE ERFUNCF (R, LAM, ERFUNC1, ERFUNC2)
C
C     DOUBLE PRECISION R
C     DOUBLE COMPLEX LAMB(2), ERFUNC1(2,2), ERFUNC2(2,2)
C  INPUTS
C           R -- ratio of astrong(mu)/astrong(mu0)
C         LAM -- eigenvalues of LO evolution operator
C  OUTPUT
C     ERFUNC1 -- matrix {\mathcal R}(mu, mu0 | 1}_{ab}
C     ERFUNC2 -- matrix {\mathcal R}(mu, mu0 | 2}_{ab}
C  PARENTS
C      EVOLF
C  CHILDREN
C      LAMBDAF
C  SOURCE
C

      SUBROUTINE ERFUNCF (R, LAMN, LAMK, ERFUNC1, ERFUNC2)

      IMPLICIT NONE
      DOUBLE PRECISION R
      DOUBLE COMPLEX LAMN(2), LAMK(2), ERFUNC1(2,2), ERFUNC2(2,2)
      INTEGER A, B
      DOUBLE PRECISION RINV
      DOUBLE COMPLEX LAMNB(2), LAMKB(2)
      INCLUDE 'header.f'


      RINV = 1.0d0 / R

      DO 10 A = 1, 2
        LAMNB(A) = LAMN(A) / BETA0(NF)
        LAMKB(A) = LAMK(A) / BETA0(NF)
 10   CONTINUE

      DO 20 A = 1, 2
      DO 20 B = 1, 2
      ERFUNC1(A,B)=( 1.0d0 - RINV**(1.0d0 + LAMNB(A) - LAMKB(B)) )/
     &      (1.0d0 + LAMNB(A) - LAMKB(B))
      IF (P. GE. 2) THEN
      ERFUNC2(A,B)=( 1.0d0 - RINV**(2.0d0 + LAMNB(A) - LAMKB(B)) )/
     &      (2.0d0 + LAMNB(A) - LAMKB(B))
      END IF
 20   CONTINUE

      RETURN
      END
C     *****
