C     ****h* gepard/deriv.f
C  FILE DESCRIPTION
C    calculation of derivative of a function
C
C    $Id$
C  NOTES
C    Adapted from Sect. 5.7 of Press et al., "Numerical Recipes".
C     *******


C     ****s* deriv.f/DERIV
C  NAME
C     DERIV  --   derivative of a function
C  DESCRIPTION
C    calculates derivative using Ridders-Neville's algorithm
C  SYNOPSIS
C     SUBROUTINE DERIV(FUNC, X, H, NEVALS, DRV, ERR)
C
C     INTEGER NEVALS
C     DOUBLE PRECISION FUNC, X, H, DRV, ERR
C     EXTERNAL FUNC
C  INPUTS
C        FUNC --  function to be differentiated. Should be
C                 declared EXTERNAL in the calling program
C           X --  differentiation point
C           H --  initial scale for f(x+h)-f(x-h)/(2h). Shouldn't
C                 be very small. (FUNC should change substantially
C                 over H.)
C      NEVALS --  number of steps (FUNC is evaluated 2*NEVALS times)
C  OUTPUT
C         DRV -- value of derivative FUNC'(X)
C         ERR -- estimated error
C  IDENTIFIERS
C         CON -- how much is scale H decreased each step
C        SAFE -- return when error is SAFE worse than the best so far
C  PARENTS
C       SCALEDEP
C  SOURCE
C

      SUBROUTINE DERIV(FUNC, X, H, NEVALS, DRV, ERR)

      IMPLICIT NONE
      INTEGER NEVALS
      INTEGER I, J
      DOUBLE PRECISION FUNC, X, H, DRV, ERR
      DOUBLE PRECISION CON, CON2, BIG, SAFE
      DOUBLE PRECISION ERRT, FAC, HH, A(NEVALS, NEVALS)
      PARAMETER (CON=1.4d0, CON2=CON*CON, BIG=1.d30, SAFE=2.0d0)
      EXTERNAL FUNC

      IF (H .EQ. 0.0d0) PAUSE 
      HH = H
      A(1,1) = (FUNC(X+HH)-FUNC(X-HH))/(2.0d0*HH)
      ERR=BIG
      DO 20 I=2,NEVALS
          HH=HH/CON
          A(1,I)=(FUNC(X+HH)-FUNC(X-HH))/(2.0d0*HH) 
          FAC=CON2
          DO 10 J=2,I
          A(J,I)=(A(J-1,I)*FAC-A(J-1,I-1))/(FAC-1.0d0)
          FAC=CON2*FAC
          ERRT=MAX(ABS(A(J,I)-A(J-1,I)),ABS(A(J,I)-A(J-1,I-1)))
          IF (ERRT.LE.ERR) THEN 
              ERR=ERRT
              DRV=A(J,I)
          END IF
 10   CONTINUE
      IF (ABS(A(I,I)-A(I-1,I-1)) .GE. SAFE*ERR) RETURN
 20   CONTINUE

      RETURN 
      END
C     *****
