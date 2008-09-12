C     ****h* gepard/erfunc.F
C  FILE DESCRIPTION
C    calculation of scale dependent part of the evolution operator
C
C    $Id$
C     *******

C     ****s* erfunc.F/ERFUNCF
C  NAME
C     ERFUNCF  --   scale dependent part of the singlet evolution operator
C  DESCRIPTION
C    calculates scale dependent part of the singlet evolution operator
C    i.e.  {\cal R}^{(n)}_{jk} - functions, see formula in 
C      Eq. (126) in [Kumericki:2007sa]
C  SYNOPSIS

      SUBROUTINE ERFUNCF (R, LAMN, LAMK, ERFUNC1, ERFUNC2)

      IMPLICIT NONE
      DOUBLE PRECISION R
      DOUBLE COMPLEX LAMN(2), LAMK(2), ERFUNC1(2,2), ERFUNC2(2,2)

C  INPUTS
C         R -- ratio of astrong(mu)/astrong(mu0)
C   LAMN(flav) -- eigenvalues of LO evolution operator at conformal moment n
C   LAMK(flav) -- eigenvalues of LO evolution operator at conformal moment k
C         
C  OUTPUT
C     ERFUNC1(flav,flav) -- matrix {\mathcal R}(mu, mu0 | 1}
C     ERFUNC2(flav,flav) -- matrix {\mathcal R}(mu, mu0 | 2}
C  PARENTS
C      EVOLF, CB1F
C  SOURCE
C
      INTEGER A, B
      DOUBLE PRECISION RINV
      DOUBLE COMPLEX LAMNB(2), LAMKB(2)
#include "header.F"


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
