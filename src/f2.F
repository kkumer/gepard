C     ****h* gepard/f2.F
C  FILE DESCRIPTION
C    calculation of F2 DIS form factor
C
C    $Id$
C     *******


C     ****s* f2.F/F2F
C  NAME
C       F2F  --   singlet F2 DIS form factor
C  DESCRIPTION
C    calculates singlet F2 DIS form factor by numerical Mellin-Barnes 
C    integration 
C    Unlike CFF, this is defined with the singlet charge factor Q_{S}^2
C    included
C      NB!! - only one of three F2(3) array elements, F2(P)
C      is changed by single call to F2F!
C  SYNOPSIS

      SUBROUTINE F2F

C  PARENTS
C     PROCDATA
C  CHILDREN
C     LOOKUPQ
C  SOURCE
C

      IMPLICIT NONE
      INTEGER K, QIND
      DOUBLE PRECISION RES, F2IMAG
      DOUBLE COMPLEX EPH, FPW, J, FCM(2)
#include "header.F"


      CALL LOOKUPQ(QIND)

*     Fixing forward kinematics:
      MTIND = 0

      EPH = EXP ( DBLCMPLX(0.0d0, PHI) )
      RES = 0.0d0

      DO 123 K = 1, NPTS
      J = N(1, K) - 1
      FCM(1) = HGRID(MTIND, K, 1)
      FCM(2) = HGRID(MTIND, K, 2)
      FPW = CGRIDDIS(QIND,K,1) * FCM(1) + CGRIDDIS(QIND,K,2) * FCM(2)
      F2IMAG = IMAGPART( EPH / XI**J * FPW )
      RES = RES + WG(K)*F2IMAG
 123  CONTINUE

      F2(P) = RES * CHARGEFAC / PI

      RETURN
      END
C     ***
