C     ****h* gepard/auxtest.F
C  FILE DESCRIPTION
C    Simple program for intermediate testing.
C
C    $Id: test.F 12 2006-07-12 08:36:24Z kuk05260 $
C     *******

C     ****p* auxtest.F/AUXTEST
C  NAME
C    AUXTEST  --  Test routine
C            
C  DESCRIPTION
C             Calculates LO and NLO CFFs corresponding
C             to testpoints A and B from letter-num.nb
C             and gepard_devel.nb
C  SYNOPSIS

      PROGRAM AUXTEST

C
C  CHILDREN
C      READPAR, INIT, CFF
C
C  SOURCE
C
      IMPLICIT NONE
#include "header.F"

      PID = 1

      CALL READPAR
*     Scales 

      Q02 = 2.5d0
      Q2 = 2.5d0


*     1. Point "A"
      WRITE (*,*) " --- Point A  (MSBAR) ----- "

      XI = 0.01
      DEL2 = -0.5
      ANSATZ = 'HARD'
      SCHEME = 'MSBDI'



      WRITE (*,*) "Expect NLO: 0.248594238 + 102.323245 I"
      P = 1
      CALL INIT
      CALL CFFF 
      WRITE (*,*) CFF(1) / CHARGEFAC
      WRITE (*,*) 
      P = 0
      CALL INIT
      CALL CFFF 
      WRITE (*,*) "Expect LO: 13.919751 + 189.278946 I "
      WRITE (*,*) CFF(0) / CHARGEFAC
      WRITE (*,*) 
      

*     2. Point "B"
      WRITE (*,*) " --- Point B  (CSBAR) ----- "

      XI = 0.00001
      DEL2 = 0.0d0
      ANSATZ = 'HARD'
      SCHEME = 'CSBAR'



      WRITE (*,*) "Expect NLO: -25716.505846 + 1134612.30 I"
      P = 1
      CALL INIT
      CALL CFFF 
      WRITE (*,*) CFF(1) / CHARGEFAC
      WRITE (*,*) 
      P = 0
      CALL INIT
      CALL CFFF 
      WRITE (*,*) "Expect LO: 399725.2294640 + 2520521.6 I "
      WRITE (*,*) CFF(0) / CHARGEFAC
      WRITE (*,*) 



      STOP
      END
C     ****
