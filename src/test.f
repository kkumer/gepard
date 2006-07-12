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
      INTEGER SPEED, P, NF
      DOUBLE PRECISION NG, NSEA, MG, MSEA, ALPHA0G, ALPHA0SEA
      DOUBLE PRECISION FITPAR(10)
      DOUBLE PRECISION W2, XI, DEL2, Q2, Q02
      DOUBLE PRECISION F2(0:2)
      DOUBLE PRECISION PARSIGMA, SIGMA

      COMMON / PARINT     /  SPEED, P, NF

      COMMON / FITPARAMS  /  FITPAR
      COMMON / KINEMATICS /  XI, DEL2, Q2, Q02
      COMMON / F2P        /  F2

      NG = 0.4
      NSEA = 2./3. - NG
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

      CALL READPAR
      CALL INIT

      P = 0

      WRITE (*, *) " ---  Test 1: Naive parton model F2 ---- "
      XI = 0.002d0
      Q2 = 1.d0
      WRITE (*, *) " 0.65783876286 is correct result"
      CALL F2F
      WRITE (*, *) F2(0)

      WRITE (*, *) " ---  Test 2: LO partial sigma, no evolution ---- "
      W2 = 82.0d0**2
      Q2 = 1.0d0
      XI = Q2 / ( 2.0d0 * W2 + Q2)
      WRITE (*, *) " 113570.49343 is result from gepard_devel.nb"
      WRITE (*, *) PARSIGMA(0.d0)

      P = 1

      WRITE (*, *) " ---  Test 3: NLO partial sigma, 1->3 GeV^2   ---- "
      W2 = 82.0d0**2
      Q2 = 3.0d0
      XI = Q2 / ( 2.0d0 * W2 + Q2)
      WRITE (*, *) " 6676.08735 is t=0 result from gepard_devel.nb"
      WRITE (*, *) PARSIGMA(0.d0)
      WRITE (*, *) " 97.346871 is t=-0.5 result from gepard_devel.nb"
      WRITE (*, *) PARSIGMA(0.5d0)

      WRITE (*, *) " ---  Test 4: NLO sigma, 1->3 GeV^2, xi~1e-4  --- "
      W2 = 82.0d0**2
      Q2 = 3.0d0
      XI = Q2 / ( 2.0d0 * W2 + Q2)
      WRITE (*, *) " 745.02480514 is result from gepard_devel.nb"
      WRITE (*, *) SIGMA()

      WRITE (*, *) " ---  Test 5: NLO sigma, 1->3 GeV^2, xi~1e-1  --- "
      W2 = 3.5d0**2
      Q2 = 3.d0
      XI = Q2 / ( 2.0d0 * W2 + Q2)
      WRITE (*, *) " 138.94079 is result from gepard_devel.nb"
      WRITE (*, *) SIGMA()


      STOP
      END
C     ****
