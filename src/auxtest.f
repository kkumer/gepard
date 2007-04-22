C     ****h* gepard/auxtest.f
C  FILE DESCRIPTION
C    Simple program for intermediate testing.
C
C    $Id: test.f 12 2006-07-12 08:36:24Z kuk05260 $
C     *******

C     ****p* auxtest.f/AUXTEST
C  NAME
C    AUXTEST  --  Test routine
C            
C  DESCRIPTION
C             Calculates LO and NLO CFFs corresponding
C             to testpoints A and B from letter-num.nb
C             and gepard_devel.nb
C
C  CHILDREN
C      READPAR, INIT, CFF
C
C  SOURCE
C


      PROGRAM AUXTEST

      IMPLICIT NONE
      INTEGER PT, NPOINTS, LN, NDEL
      DOUBLE PRECISION LOGXI, LOGXISTART, LOGXIEND, LOGXISTEP
      PARAMETER ( NPOINTS = 2 )
      DOUBLE PRECISION POINTS(4, 6, NPOINTS)
      DOUBLE PRECISION XIS(NPOINTS)
      PARAMETER ( LOGXISTART = -5.0d0, LOGXIEND = -1.30103d0,
     &       LOGXISTEP = (LOGXIEND - LOGXISTART) / (NPOINTS - 1)  )
      INCLUDE 'header.f'

      PROCESS = 'DVCS'
      FFTYPE = 'SINGLET'

      CALL READPAR
*     Scales 

      Q2 = 2.5d0

      NQS = 1
      QS(1) = Q2

      PAR(1) = 2.5d0
      PAR(2) = 0.05d0
      PAR(3) = 2.5d0

      PAR(50) = 0.0d0
      PAR(51) = 0.0d0

*     1. Point "A"
      WRITE (*,*) " --- Point A  (MSBAR) ----- "

      DEL2 = -0.5d0
      XI = 0.01
      ANSATZ = 'HARD'
      P = 1
      SCHEME = 'MSBAR'

      CALL INIT

      CALL EVOLC(1, 1)
      CALL EVOLC(2, 1)

      CALL CFFF 
      P = P - 1
      CALL CFFF 

      WRITE (*,*) "Expect NLO: 0.8470317 + 345.3313465 I"
      WRITE (*,*) CFF(1)
      WRITE (*,*) 
      WRITE (*,*) "Expect LO: 46.980589576 + 638.808718636 I "
      WRITE (*,*) CFF(0)
      WRITE (*,*) 

*     2. Point "B"
      WRITE (*,*) " --- Point B  (CSBAR) ----- "

      DEL2 = 0.0d0
      XI = 0.00001
      ANSATZ = 'HARD'
      P = 1
      SCHEME = 'CSBAR'

      CALL INIT

      CALL CFFF 
      P = P - 1
      CALL CFFF 

      WRITE (*,*) "Expect NLO: -25716.505846 + 1134612.30 I"
      WRITE (*,*) CFF(1)
      WRITE (*,*) 
      WRITE (*,*) "Expect LO: 399725.2294640 + 2520521.6 I "
      WRITE (*,*) CFF(0)
      WRITE (*,*) 



      STOP
      END
C     ****
