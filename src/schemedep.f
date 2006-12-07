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


      NF = 4
      ANSATZ = 'NSFIT'

*       1  Q02       
        PAR(1) =   2.5d0    
*       2  AS0       
        PAR(2) =   0.05d0
*       3  MU02      
        PAR(3) =   2.5d0 

* ------------ ANSATZ  -------------
* ----  11 NS   --------------------
!        PAR(11) =  see below
*       12 AL0S      
        PAR(12) =  1.1d0 
*       13 ALPS      
        PAR(13) =  0.15d0
*       14 M02S      
        PAR(14) =  (2.0d0 * MP)**2
*       15 DELM2S    
        PAR(15) =  MP**2
*       16 PS        
        PAR(16) =  3.0d0
* ----  21 NG   --------------------
!        PAR(21) =  0.5d0
*       22 AL0G      
!        PAR(22) =  1.0d0 
*       23 ALPG      
        PAR(23) =  0.15d0
*       24 M02G      
        PAR(24) =  (2.0d0 * MP)**2
*       25 DELM2G    
        PAR(25) =  MP**2
*       26 PG        
        PAR(26) =  2.0d0
* ----  31 NU   --------------------
        PAR(31) =  2.0d0 
*       32 AL0U      
        PAR(32) =  0.5d0 
*       33 ALPU      
        PAR(33) =  1.0d0 
*       34 M02U      
        PAR(34) =  (2.0d0 * MP)**2
*       35 DELM2U    
        PAR(35) =  MP**2
*       36 PU        
        PAR(36) =  1.0d0 
* ----  41 ND   --------------------
        PAR(41) =  1.0d0 
*       42 AL0D      
        PAR(42) =  0.5d0 
*       43 ALPD      
        PAR(43) =  1.0d0 
*       44 M02D      
        PAR(44) =  (2.0d0 * MP)**2
*       45 DELM2D    
        PAR(45) =  MP**2
*       46 PD        
        PAR(46) =  1.0d0    


*     Scales 

      Q2 = 10.0d0

*   Values below are for XI = 0.001
      XI = 1.0e-3
      DEL2 = 0.0d0
      SCHEME = 'MSBAR'
*   sea or not?
      PAR(11) = 0.0d0

      WRITE (*,*) "Ansatz = ", ANSATZ
      WRITE (*,*) "Delta^2 =  ", DEL2
      WRITE (*,*) "XI = ", XI

      PAR(1) = 10.0d0
      P = 0
      CALL INIT
      CALL CFFF 
      WRITE (*,*) "Expect LO (noevol): -116.731 + 130.145 I"
      WRITE (*,*) "MSBAR LO (noevol)"
      WRITE (*,*) CFF(0)
      WRITE (*,*) 

      WRITE (*,*) "Expect LO (evol): -133.895 + 160.002 I"
      WRITE (*,*) "MSBAR LO (evol):"
      PAR(1) = 2.5d0
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
      PAR(1) = 2.5d0
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
      PAR(1) = 10.0d0
      CALL INIT
      CALL CFFF 
      WRITE (*,*) "CSBAR NLO (noevol):"
      WRITE (*,*) CFF(1)
      WRITE (*,*) 

      PAR(1) = 2.5d0
      CALL INIT
      CALL CFFF 
      WRITE (*,*) "CSBAR NLO (evol):"
      WRITE (*,*) CFF(1)
      WRITE (*,*) 

      STOP
      END
C     ****
