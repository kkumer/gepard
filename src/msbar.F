C     ****h* gepard/msbar.F
C  FILE DESCRIPTION
C    calculation of  singlet Wilson coefficients for DVCS in MSbar scheme
C
C    $Id$
C     *******


C     ****s* msbar.F/MSBARF
C  NAME
C     MSBARF  --   singlet Wilson coefficients (in MSbar scheme)
C  DESCRIPTION
C    calculates singlet Wilson coefficients for DVCS in MSbar scheme
C    according to formula Eqs. (127) and (129) of [Kumericki:2007sa]
C  SYNOPSIS

      SUBROUTINE MSBARF (SEC, K, BIGCTMP )

      IMPLICIT NONE
      INTEGER  SEC, K
      DOUBLE COMPLEX BIGCTMP(0:2,2), MCU(3), KMOM

C  INPUTS
C           K -- conformal moment index
C  OUTPUT
C       BIGCTMP(P, flavour) -- Wilson coefficients
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

      IF ( PROCESS(:4) .EQ. 'DVMP' ) THEN

        BIGCTMP(0, 1) = (1.0d0,0.0d0)/NF
        BIGCTMP(0, 2) = (1.0d0,0.0d0)

        IF (P .GE. 1) THEN

*         FIXME: Approximate with asymptotic DA only: K=0
          KMOM = (0.0d0, 0.0d0)
          CALL CDVEMF (J, KMOM, MCU)
          BIGCTMP(1, 1) = MCU(1)/NF + MCU(2)
          BIGCTMP(1, 2) = MCU(3)
        END IF

      ELSE
*     --- DVCS or DIS ---

        BIGCTMP(0, 1) = (1.0d0,0.0d0)
        BIGCTMP(0, 2) = (0.0d0,0.0d0)

        BIGCTMP(1, 1) =  CF * (2.0d0 * HS1(J+1)**2 - (9.0d0/2.0d0)
     &      + (5.0d0 - 4.0d0 * HS1(J+1)) / (2.0d0 * (J+1) * (J+2))
     &   + 1.0d0 / ( (J+1)**2 * (J+2)**2 ) ) - 
     &   GAM(SEC,K,0,1,1)*LRF2/ 2.0d0

        BIGCTMP(1, 2) =  - NF * ( (4.0d0 + 3.0d0 * J + J*J) * (HS1(J)
     &      + HS1(J+2)) + 2.0d0 + 3.0d0 * J + J*J ) / ((J+1) * (J+2)
     &      * (J+3)) - GAM(SEC,K,0,1,2)*LRF2/ 2.0d0

      END IF

      RETURN
      END
C     *****
