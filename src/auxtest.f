
      PROGRAM AUXTEST

      IMPLICIT NONE
      INTEGER P, PMAX, PT, NPOINTS, LN, NDEL
      DOUBLE PRECISION XI, DEL2, Q2, Q02
      DOUBLE COMPLEX CFF(0:2)
      CHARACTER SCHEME*5, ANSATZ*6
      DOUBLE PRECISION LOGXI, LOGXISTART, LOGXIEND, LOGXISTEP
      DOUBLE PRECISION DCARG
      PARAMETER ( NPOINTS = 2 )
      DOUBLE PRECISION POINTS(4, 6, NPOINTS)
      DOUBLE PRECISION XIS(NPOINTS)
      PARAMETER ( LOGXISTART = -5.0d0, LOGXIEND = -1.30103d0,
     &       LOGXISTEP = (LOGXIEND - LOGXISTART) / (NPOINTS - 1)  )

*     Output common-blocks 

      COMMON / INITPAR    /  PMAX
      COMMON / KINEMATICS /  XI, DEL2, Q2, Q02
      COMMON / APPROX     /  P
      COMMON / LABELS     /  SCHEME, ANSATZ
      COMMON / CFF        /  CFF

      PMAX = 1

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
