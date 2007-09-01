C     ****h* gepard/msbar.f
C  FILE DESCRIPTION
C    calculation of  Wilson coefficients for DVCS in MSbar scheme
C    according to KMKPS06 paper
C
C    $Id$
C     *******


C     ****s* msbar.f/MSBARF
C  NAME
C     MSBARF  --   Wilson coefficients (in MSbar scheme)
C  DESCRIPTION
C    calculates Wilson coefficients for DVCS in MSbar scheme
C    according to KMKPS06 paper
C  SYNOPSIS
C     SUBROUTINE MSBARF (J, BIGC0, BIGC1)
C
C     DOUBLE COMPLEX J, BIGC0(2), BIGC1(2)
C  INPUTS
C           J -- conformal moment
C  OUTPUT
C       BIGC0 -- vector of two C^{(0)} (quark and gluon) Wilson 
C                coefficients (trivially equal to (1,0))
C       BIGC1 -- vector  C^{(1)} 
C  IDENTIFIERS
C       WGAMMA -- common block with moments of anomalous dimensions of DIS
C  PARENTS
C      INIT
C  CHILDREN
C      HS1
C  SOURCE
C

      SUBROUTINE MSBARF ( K, BIGCTMP )

      IMPLICIT NONE
      INTEGER  K
      DOUBLE COMPLEX BIGCTMP(0:2,2)
      DOUBLE PRECISION LRF2
      DOUBLE COMPLEX J, HS1, BIGC0(2), BIGC1(2)
      INCLUDE 'header.f'

      J = N(K) - 1

      LRF2 = LOG(RF2)

      BIGCTMP(0, 1) = (1.0d0,0.0d0)
      BIGCTMP(0, 2) = (0.0d0,0.0d0)

      BIGCTMP(1, 1) =  CF * (2.0d0 * HS1(J+1)**2 - (9.0d0/2.0d0)
     &      + (5.0d0 - 4.0d0 * HS1(J+1)) / (2.0d0 * (J+1) * (J+2))
     &   + 1.0d0 / ( (J+1)**2 * (J+2)**2 ) ) - GAM(K,0,1,1)*LRF2/ 2.0d0

      BIGCTMP(1, 2) =  - NF * ( (4.0d0 + 3.0d0 * J + J*J) * (HS1(J)
     &      + HS1(J+2)) + 2.0d0 + 3.0d0 * J + J*J ) / ((J+1) * (J+2)
     &      * (J+3)) - GAM(K,0,1,2)*LRF2/ 2.0d0

      RETURN
      END
C     *****
