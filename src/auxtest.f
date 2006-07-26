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
      INTEGER SPEED, ACC, P, NF
      CHARACTER SCHEME*5, ANSATZ*6
      INTEGER PT, NPOINTS, LN, NDEL
      DOUBLE PRECISION XI, DEL2, Q2, Q02
      DOUBLE COMPLEX CFF(0:2)
      DOUBLE PRECISION LOGXI, LOGXISTART, LOGXIEND, LOGXISTEP
      PARAMETER ( NPOINTS = 2 )
      DOUBLE PRECISION POINTS(4, 6, NPOINTS)
      DOUBLE PRECISION XIS(NPOINTS)
      PARAMETER ( LOGXISTART = -5.0d0, LOGXIEND = -1.30103d0,
     &       LOGXISTEP = (LOGXIEND - LOGXISTART) / (NPOINTS - 1)  )

*     Input common-blocks 

      COMMON / PARINT /  SPEED, ACC, P, NF
      COMMON / PARCHR /  SCHEME, ANSATZ

*     Output common-blocks 

      COMMON / KINEMATICS /  XI, DEL2, Q2, Q02
      COMMON / CFF        /  CFF

      CALL READPAR
      ACC = 6
*     Scales 

      Q2 = 2.5d0
      Q02 = 2.5d0

*     1. Point "A"
      WRITE (*,*) " --- Point A  (MSBAR) ----- "

      DEL2 = -0.5d0
      XI = 0.01
      ANSATZ = 'HARD'
      P = 1
      SCHEME = 'MSBAR'

      CALL INIT

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
