C     ****h* gepard/test.f
C  FILE DESCRIPTION
C    Simple program for testing.
C
C    $Id$
C     *******

C     ****p* test.f/TEST
C  NAME
C    TEST  --  Testing gepard routines
C            
C  DESCRIPTION
C             Calculates various stuff for fixed 
C             parameters and prints out together
C             with confirmed results from old
C             Mathematica version of gepard
C
C  CHILDREN
C      READPAR, INIT, SIGMA, PARSIGMA, F2F
C
C  SOURCE
C


      PROGRAM TEST

      IMPLICIT NONE
      INTEGER SPEED, ACC, P, NF
      DOUBLE PRECISION NG, NSEA, MG, MSEA, ALPHA0G, ALPHA0SEA
      DOUBLE PRECISION FITPAR(10)
      DOUBLE PRECISION W2, XI, DEL2, Q2, Q02, AUX
      DOUBLE PRECISION F2(0:2)
      DOUBLE PRECISION PARSIGMA, SIGMA

      COMMON / PARINT     /  SPEED, ACC, P, NF

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
      WRITE (*, *) "For values of test parameters see src/test.f"
      WRITE (*, *) 

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
      WRITE (*, *) " 5608.4194288225 is result from gepard_devel.nb"
      AUX = PARSIGMA(0.d0)
      WRITE (*, *) AUX

      P = 1

      WRITE (*, *) " ---  Test 3: NLO partial sigma, 1->3 GeV^2   ---- "
      W2 = 82.0d0**2
      Q2 = 3.0d0
      XI = Q2 / ( 2.0d0 * W2 + Q2)
      WRITE (*, *) " 329.68332610735 is t=0 result from gepard_devel.nb"
      AUX = PARSIGMA(0.d0)
      WRITE (*, *) AUX
      WRITE (*, *) " 4.807252903 is t=-0.5 result from gepard_devel.nb"
      AUX = PARSIGMA(0.5d0)
      WRITE (*, *) AUX

      WRITE (*, *) " ---  Test 4: NLO sigma, 1->3 GeV^2, xi~1e-4  --- "
      W2 = 82.0d0**2
      Q2 = 3.0d0
      XI = Q2 / ( 2.0d0 * W2 + Q2)
      WRITE (*, *) " 36.7913484 is result from gepard_devel.nb"
      AUX = SIGMA ()
      WRITE (*, *) AUX

      WRITE (*, *) " ---  Test 5: NLO sigma, 1->3 GeV^2, xi~1e-1  --- "
      W2 = 3.5d0**2
      Q2 = 3.d0
      XI = Q2 / ( 2.0d0 * W2 + Q2)
      WRITE (*, *) " 6.8612738 is result from gepard_devel.nb"
      AUX = SIGMA ()
      WRITE (*, *) AUX


      STOP
      END
C     ****
