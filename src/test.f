C     ****h* gepard/test.f
C  FILE DESCRIPTION
C    Simple program for testing.
C
C    $Id$
C     *******

C     ****p* plotsigma.f/PLOTSIGMA
C  NAME
C    TEST  --  Tests routines
C            
C  DESCRIPTION
C             Calculates various stuff for fixed 
C             paramterers.
C
C  CHILDREN
C      SIGMA
C
C  SOURCE
C


      PROGRAM TEST

      IMPLICIT NONE
      INTEGER N, NMAX, IFLAG, I
      DOUBLE PRECISION FITPAR(10)
      DOUBLE PRECISION F2(0:2)
      DOUBLE PRECISION W2, XI, DEL2, Q2, Q02
      DOUBLE PRECISION NG, NSEA, MG, MSEA, ALPHA0G, ALPHA0SEA

      DOUBLE PRECISION SIGMA

      COMMON / KINEMATICS /  XI, DEL2, Q2, Q02
      COMMON / FITPARAMS  /  FITPAR
      COMMON / F2P        /  F2


      NG = 0.4
      NSEA = 0.2
      MG = 1.2
      MSEA = 1.1
      ALPHA0G = 1.2
      ALPHA0SEA = 1.1
      Q02 = 1.0


      FITPAR(1) = NG 
      FITPAR(2) = NSEA 
      FITPAR(3) = MG 
      FITPAR(4) = MSEA 
      FITPAR(5) = ALPHA0G 
      FITPAR(6) = ALPHA0SEA 
      FITPAR(7) = Q02 




      WRITE (*, *) " ------  Test 1 ----- "
      W2 = 82.0d0**2
      Q2 = 3.0d0
      XI = Q2 / ( 2.0d0 * W2 + Q2)
      WRITE (*, *) " 407.598953 was old result"
      WRITE (*, *) SIGMA()

      WRITE (*, *) " ------  Test 2 ----- "
      XI = 0.002d0
      Q2 = 3.0d0
      WRITE (*, *) " 0.376324871 was old result"
      CALL F2F
      WRITE (*, *) F2(1)


      STOP
      END
C     ****
