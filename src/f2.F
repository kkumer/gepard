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
C     GETQIND
C  SOURCE
C

      IMPLICIT NONE
      INTEGER K, QIND, GETQIND
      DOUBLE PRECISION RES, F2IMAG
      DOUBLE COMPLEX EPH, FPW, J
#include "header.F"


      PID = 0
      QIND = GETQIND()
      CALL GETMBGPD

*     Fixing forward kinematics:
      MTIND = 0

      EPH = EXP ( DBLCMPLX(0.0d0, PHI) )
      RES = 0.0d0

      DO 123 K = 1, NPTS
      J = N(1, K) - 1
      FPW = WCE(0,1,QIND,K,1) * GPDH(K,1) + 
     &      WCE(0,1,QIND,K,2) * GPDH(K,2)
      F2IMAG = IMAGPART( EPH / XI**J * FPW )
      RES = RES + WG(K)*F2IMAG
 123  CONTINUE

      IF (NF .EQ. 3) THEN
         CHARGEFAC =  2.0D0 / 9.0D0
      ELSE IF (NF .EQ. 4) THEN
         CHARGEFAC =  5.0D0 / 18.0D0
      ELSE
         STOP
      END IF

      F2(P) = RES * CHARGEFAC / PI

      RETURN
      END
C     ***
