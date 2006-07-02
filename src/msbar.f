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
C     SUBROUTINE MSBARF (NF, J, RF2, RR2, BIGC0, BIGC1)
C
C     INTEGER NF
C     DOUBLE PRECISION RF2, RR2
C     DOUBLE COMPLEX J, BIGC0(2), BIGC1(2)
C  INPUTS
C          NF -- number of active flavours
C           J -- conformal moment
C         RF2 -- ratio {\cal Q}^2/{\mu_{f}^2} of photon
C                virtuality over factorization scale
C         RR2 -- ratio {\cal Q}^2/{\mu_{r}^2} of photon
C                virtuality over renormalization scale
C  OUTPUT
C       BIGC0 -- vector of two C^{(0)} (quark and gluon) Wilson 
C                coefficients (trivially equal to (1,0))
C       BIGC1 -- vector  C^{(1)} 
C  IDENTIFIERS
C       WGAMMA -- common block with moments of anomalous dimensions of DIS
C  PARENTS
C      PARWAVF
C  CHILDREN
C      HS1
C  NOTES
C      One needs to call COMMONF first, to initialize common block WGAMMA
C  SOURCE
C

      SUBROUTINE MSBARF (NF, J, RF2, RR2, BIGC0, BIGC1)

      IMPLICIT NONE
      INTEGER NF, NFMIN, NFMAX, K
      DOUBLE PRECISION RF2, LRF2, RR2, CF
      DOUBLE COMPLEX J, GAM0(2,2), GAM1(2,2), GAM2(2,2)
      DOUBLE COMPLEX BIGC0(2), BIGC1(2)
      DOUBLE COMPLEX HS1
      PARAMETER (NFMIN = 3, NFMAX = 6)

*   Input common-block

      COMMON / WGAMMA  /  GAM0, GAM1, GAM2

      LRF2 = LOG(RF2)
      CF = 4.0d0/3.0d0

      BIGC0(1) = (1.0d0,0.0d0)
      BIGC0(2) = (0.0d0,0.0d0)

      BIGC1(1) =  CF * (2.0d0 * HS1(J+1)**2 - (9.0d0/2.0d0)
     &      + (5.0d0 - 4.0d0 * HS1(J+1)) / (2.0d0 * (J+1) * (J+2))
     &      + 1.0d0 / ( (J+1)**2 * (J+2)**2 ) ) - GAM0(1,1)*LRF2/ 2.0d0

      BIGC1(2) =  - NF * ( (4.0d0 + 3.0d0 * J + J*J) * (HS1(J)
     &      + HS1(J+2)) + 2.0d0 + 3.0d0 * J + J*J ) / ((J+1) * (J+2)
     &      * (J+3)) - GAM0(1,2)*LRF2/ 2.0d0

      RETURN
      END
C     *****
