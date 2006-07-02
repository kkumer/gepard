C     ****h* gepard/dcomplex.f
C  FILE DESCRIPTION
C    Tangent and argument of double complex number Z
C
C    $Id: dcomplex.f,v 1.1 2006-04-12 15:47:13+02 kkumer Exp kkumer $
C
C  NOTES
C    Adapted from SLATEC routines CTAN and CARG
C     *******


C     ****f* dcomplex.f/DCTAN
C  NAME
C    DCTAN  --  tangent function
C  SYNOPSIS
C    DOUBLE COMPLEX FUNCTION DCTAN (Z)
C
C    DOUBLE COMPLEX Z
C  PARENTS
C     FREAL
C  CHILDREN
C     D1MACH, XERCLR, XERMSG  -- from SLATEC library
C  SOURCE
C

      DOUBLE COMPLEX FUNCTION DCTAN (Z)

      IMPLICIT NONE
      DOUBLE COMPLEX Z
      DOUBLE PRECISION D1MACH, SQEPS, X2, Y2, SN2X, DEN
      SAVE SQEPS
      DATA SQEPS /0.0d0/
      IF (SQEPS.EQ.0.0d0) SQEPS = SQRT (D1MACH(4))

      X2 = 2.0d0*DBLE(Z)
      Y2 = 2.0d0*DIMAG(Z)

      SN2X = SIN (X2)
      CALL XERCLR

      DEN = COS(X2) + COSH(Y2)
      IF (DEN .EQ. 0.0d0) CALL XERMSG ('SLATEC', 'DCTAN',
     &   'TAN IS SINGULAR FOR INPUT Z (X IS PI/2 OR 3*PI/2 AND Y IS 0)',
     &   2, 2)

      IF (ABS(DEN).GT.MAX(ABS(X2),1.)*SQEPS) GO TO 10
      CALL XERCLR
      CALL XERMSG ('SLATEC', 'DCTAN',
     &   'ANSWER LT HALF PRECISION, ABS(X) TOO BIG OR X TOO NEAR ' //
     &   'PI/2 OR 3*PI/2', 1, 1)

 10   DCTAN = COMPLEX (SN2X/DEN, SINH(Y2)/DEN)

      RETURN
      END
C     ***

C     ****f* dcomplex.f/DCARG
C  NAME
C     DCARG  --   argument of a complex number
C  SYNOPSIS
C     DOUBLE PRECISION FUNCTION DCARG (Z)
C
C     DOUBLE COMPLEX Z
C  SOURCE
C

      DOUBLE PRECISION FUNCTION DCARG (Z)

      IMPLICIT NONE
      DOUBLE COMPLEX Z
      DCARG = 0.0d0
      IF (DBLE(Z).NE.0.0d0 .OR. DIMAG(Z).NE.0.0d0) DCARG =
     &  ATAN2 (DIMAG(Z), DBLE(Z))

      RETURN
      END
C     ***
