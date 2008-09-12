C     ****h* gepard/dcomplex.F
C  FILE DESCRIPTION
C    Tangent and argument of double complex number Z
C
C    $Id$
C
C  NOTES
C    Adapted from SLATEC routines CTAN and CARG
C     *******


C     ****f* dcomplex.F/G77TAN
C  NAME
C    G77TAN  --  tangent function
C  SYNOPSIS

      DOUBLE COMPLEX FUNCTION G77TAN (Z)

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
      Y2 = 2.0d0 * IMAGPART(Z)

C     FIXME: This is a kludge
      IF ( ABS(Y2) .GT. 3.0D2 ) THEN
         G77TAN = DBLCMPLX(0.0d0, SIGN(1.0d0, Y2))
         RETURN
      END IF

      SN2X = SIN (X2)


      DEN = COS(X2) + COSH(Y2)
      IF (DEN .EQ. 0.0d0) CALL ERROR ('GeParD', 'G77TAN',
     &   'TAN IS SINGULAR FOR INPUT Z (X IS PI/2 OR 3*PI/2 AND Y IS 0)',
     &   2, 2)

      IF (ABS(DEN).GT.MAX(ABS(X2),1.0d0)*SMALL) GO TO 10
      CALL ERROR ('GeParD', 'G77TAN',
     &   'ABS(X) TOO BIG OR X TOO NEAR PI/2 OR 3*PI/2                 ',
     &   1, 1)

 10   G77TAN = DBLCMPLX(SN2X/DEN, SINH(Y2)/DEN)

      RETURN
      END
C     ***

C     ****f* dcomplex.F/DCARG
C  NAME
C     DCARG  --   argument of a complex number
C  SYNOPSIS

      DOUBLE PRECISION FUNCTION DCARG (Z)

      IMPLICIT NONE
      DOUBLE COMPLEX Z

C  SOURCE
C

      DCARG = 0.0d0
      IF (DBLE(Z).NE.0.0d0 .OR. IMAGPART(Z).NE.0.0d0) DCARG =
     &  ATAN2 (IMAGPART(Z), DBLE(Z))

      RETURN
      END
C     ***
