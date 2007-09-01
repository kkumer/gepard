C     ****h* gepard/lambda.f
C  FILE DESCRIPTION
C     eigenvalues of the LO singlet anomalous dimensions matrix
C
C    $Id$
C     *******


C     ****s* lambda.f/LAMBDAF
C  NAME
C     LAMBDAF  --  eigenvalues of the LO singlet anomalous dimensions matrix
C  DESCRIPTION
C    calculates eigenvalues of the LO singlet anomalous dimensions
C    matrix using Dieter's trick for treating the square-root cut
C    in complex plane
C  SYNOPSIS
C     SUBROUTINE LAMBDAF (K, LAM)
C
C     INTEGER K
C     DOUBLE COMPLEX LAM(2)
C  INPUTS
C           K -- Mellin-Barnes integration point index
C  OUTPUT
C         LAM -- eigenvalues (LAM(1) is "+" and LAM(2) 
C                 is "-" eigenvalue)
C  PARENTS
C      PROJECTORSF, ERFUNCF
C  SOURCE
C

      SUBROUTINE LAMBDAF (K, LAM)

      IMPLICIT NONE
      INTEGER K
      DOUBLE COMPLEX LAM(2)
      DOUBLE COMPLEX AUX
      INCLUDE 'header.f'


*  Eigenvalues. Expression is adjusted to avoid the crossing
*    of the SQRT cut on the negative real axis

      AUX = (GAM(K,0,1,1) - GAM(K,0,2,2)) * SQRT (1.0d0 + 
     &     4.0d0 * GAM(K,0,1,2)*GAM(K,0,2,1)/(GAM(K,0,1,1) 
     &          - GAM(K,0,2,2))**2)
      LAM(1) = 0.5d0 * (GAM(K,0,1,1) + GAM(K,0,2,2) - AUX)
      LAM(2) = LAM(1) + AUX

      RETURN
      END
C     *****


C     ****s* lambda.f/LAMBDANDF
C  NAME
C     LAMBDANDF  --  eigenvalues of the LO singlet anomalous dimensions matrix
C  DESCRIPTION
C    calculates eigenvalues of the LO singlet anomalous dimensions
C    matrix using Dieter's trick for treating the square-root cut
C    in complex plane - version needed for non-diagonal evolution only
C  SYNOPSIS
C     SUBROUTINE LAMBDAF (K, LAM)
C
C     INTEGER K
C     DOUBLE COMPLEX LAM(2)
C  INPUTS
C           K -- Mellin-Barnes integration point index
C  OUTPUT
C         LAM -- eigenvalues (LAM(1) is "+" and LAM(2) 
C                 is "-" eigenvalue)
C  PARENTS
C      PROJECTORSF, ERFUNCF
C  SOURCE
C

      SUBROUTINE LAMBDANDF (GAMN, GAMK, LAMN, LAMK)

      IMPLICIT NONE
      INTEGER K1
      DOUBLE COMPLEX GAMN(2,2), GAMK(2,2), LAMN(2), LAMK(2)
      DOUBLE COMPLEX AUX
      INCLUDE 'header.f'

*  Eigenvalues. Expression is adjusted to avoid the crossing
*    of the SQRT cut on the negative real axis

      AUX = (GAMN(1,1) - GAMN(2,2)) * SQRT (1.0d0 + 
     &          4.0d0 * GAMN(1,2)*GAMN(2,1)/(GAMN(1,1) 
     &          - GAMN(2,2))**2)
      LAMN(1) = 0.5d0 * (GAMN(1,1) + GAMN(2,2) - AUX)
      LAMN(2) = LAMN(1) + AUX

      AUX = (GAMK(1,1) - GAMK(2,2)) * SQRT (1.0d0 + 
     &          4.0d0 * GAMK(1,2)*GAMK(2,1)/(GAMK(1,1) 
     &          - GAMK(2,2))**2)
      LAMK(1) = 0.5d0 * (GAMK(1,1) + GAMK(2,2) - AUX)
      LAMK(2) = LAMK(1) + AUX

      RETURN
      END
C     *****
