C     ****h* gepard/wcmsbar.F
C  FILE DESCRIPTION
C    calculation of  Wilson coefficients for DVCS in MSbar scheme
C
C    $Id$
C     *******


C     ****s* wcmsbar.F/WCMSBAR
C  NAME
C     WCMSBAR  --   DVCS Wilson coefficients (in MSbar scheme)
C  DESCRIPTION
C    calculates Wilson coefficients for DVCS in MSbar scheme
C    according to formula Eqs. (127) and (129) of [Kumericki:2007sa]
C  SYNOPSIS

      SUBROUTINE WCMSBAR (SEC, K, WCTMP)

      IMPLICIT NONE
      INTEGER  SEC, K
      DOUBLE COMPLEX WCTMP(0:2,4)

C  INPUTS
C         SEC -- SO(3) PW
C           K -- conformal moment index
C  OUTPUT
C       WCTMP(P, flavour) -- Wilson coefficients
C  PARENTS
C      INIT
C  CHILDREN
C      HS1
C  SOURCE
C

      DOUBLE PRECISION LRF2
      DOUBLE COMPLEX J, HS1
#include "header.F"

      J = N(SEC, K) - 1

      LRF2 = LOG(RF2)

*     --- LO ---
      WCTMP(0, 1) = (1.0d0,0.0d0)
      WCTMP(0, 2) = (0.0d0,0.0d0)
      WCTMP(0, 3) = (1.0d0,0.0d0)
      WCTMP(0, 4) = (1.0d0,0.0d0)

*     --- NLO ---
      WCTMP(1, 1) =  CF * (2.0d0 * HS1(J+1)**2 - (9.0d0/2.0d0)
     &      + (5.0d0 - 4.0d0 * HS1(J+1)) / (2.0d0 * (J+1) * (J+2))
     &      + 1.0d0 / ( (J+1)**2 * (J+2)**2 ) ) - 
     &      GAM(SEC,K,0,1,1)*LRF2/ 2.0d0

      WCTMP(1, 2) =  - NF * ( (4.0d0 + 3.0d0 * J + J*J) * (HS1(J)
     &      + HS1(J+2)) + 2.0d0 + 3.0d0 * J + J*J ) / ((J+1) * (J+2)
     &      * (J+3)) - GAM(SEC,K,0,1,2)*LRF2/ 2.0d0

*     TODO: Check that this is always equal to WCTMP(1,1) and use that
      WCTMP(1, 3) =  CF * (2.0d0 * HS1(J+1)**2 - (9.0d0/2.0d0)
     &    + (5.0d0 - 4.0d0 * HS1(J+1)) / (2.0d0 * (J+1) * (J+2))
     &    + 1.0d0 / ( (J+1)**2 * (J+2)**2 ) ) 
     &    - GAMNS(SEC,K,0,0)*LRF2/2.0d0
*     TODO: implement NS-

      RETURN
      END
C     *****
