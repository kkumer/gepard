C     ****h* gepard/wcdvmp.F
C  FILE DESCRIPTION
C    calculation of  Wilson coefficients for DVMP
C
C    $Id$
C     *******


C     ****s* wcdvmp.F/WCDVMP
C  NAME
C     MSBARF  --   singlet Wilson coefficients (in MSbar scheme)
C  DESCRIPTION
C    calculates singlet Wilson coefficients for DVCS in MSbar scheme
C    according to formula Eqs. (127) and (129) of [Kumericki:2007sa]
C  SYNOPSIS

      SUBROUTINE WCDVMP (SEC, K, WCTMP )

      IMPLICIT NONE
      INTEGER  SEC, K
      DOUBLE COMPLEX WCTMP(0:2,4), MCU(3), KMOM

C  INPUTS
C         SEC -- SO(3) PW
C           K -- conformal moment index
C  OUTPUT
C       WCTMP(P, flavour) -- Wilson coefficients
C  PARENTS
C      INIT
C  CHILDREN
C      C1DVMP
C  SOURCE
C

      DOUBLE PRECISION LRF2
      DOUBLE COMPLEX J, HS1
#include "header.F"

c     WRITE (*,*)  "BEGIN doing WCDVMP SEC, K = ", SEC, K
      J = N(SEC, K) - 1

      LRF2 = LOG(RF2)


      WCTMP(0, 1) = (1.0d0,0.0d0)/NF
      WCTMP(0, 2) = (1.0d0,0.0d0)
*     FIXME: check division with NF???
      WCTMP(0, 3) = (1.0d0,0.0d0)/NF
      WCTMP(0, 4) = (1.0d0,0.0d0)/NF


*     FIXME: Approximate with asymptotic DA only: K=0
      KMOM = (0.0d0, 0.0d0)
      IF (P .GE. 1) THEN
        CALL C1DVMP (1, J, KMOM, MCU)
        WCTMP(1, 1) = MCU(1)/NF + MCU(2)
        WCTMP(1, 2) = MCU(3)
*       FIXME: Check NS+ and implement NS-
        WCTMP(1, 3) = MCU(1)/NF
      ENDIF

c     WRITE (*,*)  "END doing WCDVMP SEC, K = ", SEC, K
      RETURN
      END
C     *****
