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
      INTEGER SPEED, ACC, P, NF
      INTEGER NPTSMAX
      PARAMETER (NPTSMAX = 512)
      DOUBLE COMPLEX NGAM(NPTSMAX,0:2,2,2)
      DOUBLE COMPLEX AUX

*   Input common-blocks

      COMMON / PARINT /  SPEED, ACC, P, NF

      COMMON / NGAM     /  NGAM

*  Eigenvalues. Expression is adjusted to avoid the crossing
*    of the SQRT cut on the negative real axis

      AUX = (NGAM(K,0,1,1) - NGAM(K,0,2,2)) * SQRT (1.0d0 + 
     &          4.0d0 * NGAM(K,0,1,2)*NGAM(K,0,2,1)/(NGAM(K,0,1,1) 
     &          - NGAM(K,0,2,2))**2)
      LAM(1) = 0.5d0 * (NGAM(K,0,1,1) + NGAM(K,0,2,2) - AUX)
      LAM(2) = LAM(1) + AUX

      RETURN
      END
C     *****
