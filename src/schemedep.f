C     ****h* gepard/schemedep.f
C  FILE DESCRIPTION
C    Investigation of scheme dependence
C
C    $Id: test.f 12 2006-07-12 08:36:24Z kuk05260 $
C     *******

C     ****p* schemedep.f/SCHEMEDEP
C  NAME
C    SCHEMEDEP  --  Scheme dependence
C            
C  DESCRIPTION
C             Calculates LO and NLO CFFs in both MSBAR
C             and CSBAR schemes, with various choices
C             for evolution.
C
C  CHILDREN
C      READPAR, INIT, CFF
C
C  SOURCE
C


      PROGRAM SCHEMEDEP

      IMPLICIT NONE
      INTEGER PT, NPOINTS, LN, NDEL
      DOUBLE COMPLEX TMP
      DOUBLE PRECISION LOGXI, LOGXISTART, LOGXIEND, LOGXISTEP
      PARAMETER ( NPOINTS = 2 )
      DOUBLE PRECISION POINTS(4, 6, NPOINTS)
      DOUBLE PRECISION XIS(NPOINTS)
      PARAMETER ( LOGXISTART = -5.0d0, LOGXIEND = -1.30103d0,
     &       LOGXISTEP = (LOGXIEND - LOGXISTART) / (NPOINTS - 1)  )
      INCLUDE 'header.f'


      CALL READPAR
*     Scales 

      Q2 = 10.0d0
      PAR(1) = 10.0d0

*     alpha_strong 
      PAR(2) = 0.05d0
      PAR(3) = 2.5d0

*   Values below are for XI = 0.001
      XI = 1.0e-3
      DEL2 = 0.0d0
      SCHEME = 'MSBAR'
      ANSATZ = 'NSTOY'

      P = 0
      CALL INIT
      CALL CFFF 
      WRITE (*,*) "Expect LO (noevol): -116.731 + 130.145 I"
      WRITE (*,*) "MSBAR LO (noevol)"
      WRITE (*,*) CFF(0)
      WRITE (*,*) 

      WRITE (*,*) "Expect LO (evol): -133.895 + 160.002 I"
      WRITE (*,*) "MSBAR LO (evol):"
      PAR(1) = 1.0d0
      CALL INIT
      CALL CFFF 
      WRITE (*,*) CFF(0)
      WRITE (*,*) 

      P = 1
      CZERO = 0
      PAR(1) = 10.0d0
      CALL INIT
      CALL CFFF 
      WRITE (*,*) "Expect NLO (noevol): 14.3975 - 8.8217 I"
      WRITE (*,*) "MSBAR NLO (noevol):"
      WRITE (*,*) CFF(1)
      WRITE (*,*) 

      SCHEME = 'MSBLO'
      PAR(1) = 1.0d0
      WRITE (*,*) "Expect NLO (LO evol): 18.0721 - 13.6474 I"
      WRITE (*,*) "MSBAR NLO (LO evol)"
      CALL INIT
      CALL CFFF 
      WRITE (*,*) CFF(1)
      TMP = CFF(1)
      WRITE (*,*) 


      SCHEME = 'MSBAR'
      WRITE (*,*) "Expect NLO (evolD): -1.59391 + 6.58042 I"
      WRITE (*,*) "MSBAR NLO (evolD):"
      CALL INIT
      CALL CFFF 
      WRITE (*,*) CFF(1) - TMP
      TMP = CFF(1)
      WRITE (*,*) 

      SCHEME = 'MSBND'
      WRITE (*,*) "Expect NLO (evolND): -0.503399 + 0.883417 I"
      WRITE (*,*) "MSBAR NLO (evolND):"
      CALL INIT
      CALL CFFF 
      WRITE (*,*) CFF(1) - TMP
      WRITE (*,*) 

      SCHEME = 'CSBAR'

      P = 1
      CZERO = 0
      Q2 = 2.5d0
      PAR(1) = 2.5d0
      CALL INIT
      CALL CFFF 
      WRITE (*,*) "CSBAR NLO (noevol):"
      WRITE (*,*) CFF(1)
      WRITE (*,*) 

      PAR(1) = 1.0d0
      CALL INIT
      CALL CFFF 
      WRITE (*,*) "CSBAR NLO (evol):"
      WRITE (*,*) CFF(1)
      WRITE (*,*) 

      STOP
      END
C     ****
