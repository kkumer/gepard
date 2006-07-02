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
C     SUBROUTINE LAMBDAF (NF, LAM)
C
C     INTEGER NF
C     DOUBLE COMPLEX LAM(2)
C  INPUTS
C          NF -- number of active flavours
C  OUTPUT
C         LAM -- eigenvalues (LAM(1) is "+" and LAM(2) 
C                 is "-" eigenvalue)
C  PARENTS
C      PROJECTORSF, ERFUNCF
C  NOTES
C      One needs to call COMMONF first, to initialize
C      common block WGAMMA
C  SOURCE
C

      SUBROUTINE LAMBDAF (NF, LAM)

      IMPLICIT NONE
      INTEGER NF
      DOUBLE COMPLEX LAM(2)
      DOUBLE COMPLEX GAM0(2,2), GAM1(2,2), GAM2(2,2)
      DOUBLE COMPLEX AUX

*   Input common-block

      COMMON / WGAMMA  /  GAM0, GAM1, GAM2

*  Eigenvalues. Expression is adjusted to avoid the crossing
*    of the SQRT cut on the negative real axis

!      AUX = SQRT ((GAM0(1,1) - GAM0(2,2))**2  + 
!     &          4.0d0 * GAM0(1,2)*GAM0(2,1))
      AUX = (GAM0(1,1) - GAM0(2,2)) * SQRT (1.0d0 + 
     &          4.0d0 * GAM0(1,2)*GAM0(2,1)/(GAM0(1,1) - GAM0(2,2))**2)
      LAM(1) = 0.5d0 * (GAM0(1,1) + GAM0(2,2) - AUX)
      LAM(2) = LAM(1) + AUX

      RETURN
      END
C     *****
