C     ****h* gepard/test.F
C  FILE DESCRIPTION
C    Simple program for testing.
C
C    $Id$
C     *******

C     ****p* test.F/TEST
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
      INTEGER K
      DOUBLE PRECISION AUX
      DOUBLE PRECISION PARSIGMA, SIGMA
#include "header.F"

      FFTYPE  = 'SINGLET'

      CALL READPAR

*   GPD ansatz:

      PAR(21) = 0.4d0
      PAR(11) = 2./3. - PAR(21)
      PAR(12) = 1.1d0
      PAR(13) = 0.25d0
      PAR(14) = 1.1d0
      PAR(22) = 1.2d0
      PAR(23) = 0.25d0
      PAR(24) = 1.2d0

*   turn off xi-dependence
      PAR(19) = 0.0d0
      PAR(29) = 0.0d0


      P = 0
      PROCESS = 'DIS'
      DEL2 = 0.0d0
      CALL INIT
      CALL GETMBGPD
      DO 11  K = 1, NPTS
        HGRID(0, K, 1) = MBGPD(K, 1) 
        HGRID(0, K, 2) = MBGPD(K, 2)
 11   CONTINUE
      WRITE (*, *) "For values of test parameters see src/test.F"
      WRITE (*, *) 

      WRITE (*, *) " ---  Test 1: Naive parton model F2 ---- "
      XI = 0.002d0
      Q2 = 1.d0
      NQSDIS = 1
      QSDIS(1) = Q2
      CALL EVOLC(1)
      WRITE (*, *) " 0.65783876286 is correct result"
      CALL F2F
      WRITE (*, *) F2(0)

      P = 1
      PROCESS = 'DVCS'
      CALL INIT

      WRITE (*, *) " ---  Test 2: LO partial sigma, no evolution ---- "
      P = 0
      W2 = 82.0d0**2
      Q2 = 1.0d0
      NQS = 1
      QS(1) = Q2
      CALL EVOLC(1)
      XI = Q2 / ( 2.0d0 * W2 + Q2)
      WRITE (*, *) " 5608.4194288225 is result from gepard_devel.nb"
      MTIND = 0
      AUX = PARSIGMA()
      WRITE (*, *) AUX

      P = 1

      WRITE (*, *) " ---  Test 3: NLO partial sigma, 1->3 GeV^2   ---- "
      W2 = 82.0d0**2
      Q2 = 3.0d0
      NQS = 2
      QS(2) = Q2
      CALL EVOLC(2)
      XI = Q2 / ( 2.0d0 * W2 + Q2)
      WRITE (*, *) " 329.68332610735 is t=0 result from gepard_devel.nb"
      MTIND = 0
      AUX = PARSIGMA()
      WRITE (*, *) AUX
      WRITE (*, *) " 4.807252903 is t=-0.5 result from gepard_devel.nb"
      MTIND = NMTS + 3
      AUX = PARSIGMA()
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

      WRITE (*, *) " ---  Test 6: NLO sigma, 1->4 GeV^2, xi=1e-6  --- "
      Q2 = 4.0D0
      NQS = 3
      QS(3) = Q2
      CALL EVOLC(3)
      XI = 1.0D-6
      WRITE (*, *) " 99.8381042 is old gepard result"
      AUX = SIGMA ()
      WRITE (*, *) AUX

      STOP
      END
C     ****
