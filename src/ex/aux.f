C     ****h* gepard/aux.f
C  FILE DESCRIPTION
C    Auxilliary program for printing out non-singlet CFF's
C
C    $Id$
C     *******

C     ****p* aux.f/AUX
C  NAME
C    AUX  --  Auxilliary program for printing out singlet CFF's
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


      PROGRAM AUX

      IMPLICIT NONE
      INTEGER CZ
      DOUBLE COMPLEX TMP
      DOUBLE PRECISION Q02, Q2EXP, MP
      PARAMETER (MP = 0.938272d0 )
      CHARACTER SUBANSATZ*4

      INCLUDE '../header.f'

      CALL READPAR
      NF = 4
      ANSATZ = 'FIT'
      CZ = 0

*       1  Q02       
        PAR(1) =   2.5d0    
*       2  AS0       
        PAR(2) =   0.05d0
*       3  MU02      
        PAR(3) =   2.5d0 

* ------------ ANSATZ  -------------
* ----  11 NS   --------------------
!        PAR(11) =  0.2d0 
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

*       1  Q02       
!       PAR(1) =   2.5d0    
*       2  AS0       
        PAR(2) =   0.05d0
*       3  MU02      
        PAR(3) =   2.5d0 


      Q02 = 2.5d0
      Q2EXP = 10.0d0
      XI = 0.1d0
      DEL2 = -0.5d0
      SCHEME = 'MSBAR'

      SUBANSATZ = 'HARD'

      IF ( SUBANSATZ .EQ.  'HARD' ) THEN
!          'HARD'
            PAR(21) = 0.4d0
            PAR(22) = PAR(12) + 0.05d0
      ELSE
!           'SOFT'
            PAR(21) = 0.3d0
            PAR(22) = PAR(12) - 0.2d0
      END IF
      PAR(11) = (2.0d0/3.0d0) - PAR(21)


      WRITE (*,901) ANSATZ, DEL2, XI
      WRITE (*,903) Q02, Q2EXP

!      DO 10 CND = -0.9, 1.3, 0.2

      WRITE (*,904) CND

      Q2 = Q2EXP
      PAR(1) = Q2EXP
      P = 0
      CALL INIT
      CALL CFFF 
      WRITE (*,902) "LO (no evol)", CFF(0)

      PAR(1) = Q02
      CALL INIT
      CALL CFFF 
      WRITE (*,902) "LO (LO evol)", CFF(0)

      P = 1
      CZERO = CZ
      PAR(1) = Q2EXP
      CALL INIT
      CALL CFFF 
      WRITE (*,902) "MS NLO (no evol)", CFF(1)

      SCHEME = 'MSBLO'
      PAR(1) = Q02
      CALL INIT
      CALL CFFF 
      WRITE (*,902) "MS NLO (LO evol)", CFF(1)
      TMP = CFF(1)

      SCHEME = 'MSBAR'
      CALL INIT
      CALL CFFF 
      WRITE (*,902) "MS NLO (evol D)", CFF(1) - TMP
      TMP = CFF(1)

      SCHEME = 'MSBND'
      CALL INIT
      CALL CFFF 
      WRITE (*,902) "MS NLO (evol ND)", CFF(1) - TMP

      SCHEME = 'CSBAR'
      PAR(1) = Q2EXP
      CALL INIT
      CALL CFFF 
      WRITE (*,902) "CS NLO (no evol)", CFF(1)

      PAR(1) = Q02
      CALL INIT
      CALL CFFF 
      WRITE (*,902) "CS NLO (evol)", CFF(1)

 10   CONTINUE

 901  FORMAT (1X, /, 1X, 'ANSATZ = ', A, 4X, 'DELTA^2 = ', F5.2,
     &    4X,'XI = ', F8.5, /)
 902  FORMAT (1X, A16, ' H = ', F11.2, ' + ', F11.2, ' I')
 903  FORMAT (1X, 'evolution : ', F4.1, ' --> ', F5.1, ' GeV^2', /)
 904  FORMAT (1X, 'CND = ', F6.3)
      STOP
      END
C     ****
