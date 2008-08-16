C     ****h* gepard/dcomplex.f
C  FILE DESCRIPTION
C    Tangent and argument of double complex number Z
C
C    $Id$
C
C  NOTES
C    Adapted from SLATEC routines CTAN and CARG
C     *******


C     ****f* dcomplex.f/DCTAN
C  NAME
C    DCTAN  --  tangent function
C  SYNOPSIS

      DOUBLE COMPLEX FUNCTION DCTAN (Z)

      IMPLICIT NONE
      DOUBLE COMPLEX Z
#include "header.F"

C  PARENTS
C     CFFF, HJ
C  SOURCE
C

      DOUBLE PRECISION SMALL, X2, Y2, SN2X, DEN
      DATA SMALL /1.0D-8/


      X2 = 2.0d0 * DBLE(Z)
      Y2 = 2.0d0 * DIMAG(Z)

C     FIXME: This is a kludge
      IF ( ABS(Y2) .GT. 3.0D2 ) THEN
         DCTAN = DBLCMPLX(0.0d0, SIGN(1.0d0, Y2))
         RETURN
      END IF

      SN2X = SIN (X2)


      DEN = COS(X2) + COSH(Y2)
      IF (DEN .EQ. 0.0d0) CALL ERROR ('GeParD', 'DCTAN',
     &   'TAN IS SINGULAR FOR INPUT Z (X IS PI/2 OR 3*PI/2 AND Y IS 0)',
     &   2, 2)

      IF (ABS(DEN).GT.MAX(ABS(X2),1.)*SMALL) GO TO 10
      CALL ERROR ('GeParD', 'DCTAN',
     &   'ABS(X) TOO BIG OR X TOO NEAR PI/2 OR 3*PI/2                 ',
     &   1, 1)

 10   DCTAN = DBLCMPLX(SN2X/DEN, SINH(Y2)/DEN)

      RETURN
      END
C     ***

C     ****f* dcomplex.f/DCARG
C  NAME
C     DCARG  --   argument of a complex number
C  SYNOPSIS

      DOUBLE PRECISION FUNCTION DCARG (Z)

      IMPLICIT NONE
      DOUBLE COMPLEX Z

C  SOURCE
C

      DCARG = 0.0d0
      IF (DBLE(Z).NE.0.0d0 .OR. DIMAG(Z).NE.0.0d0) DCARG =
     &  ATAN2 (DIMAG(Z), DBLE(Z))

      RETURN
      END
C     ***
