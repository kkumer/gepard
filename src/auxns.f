C     ****h* gepard/auxns.f
C  FILE DESCRIPTION
C    Auxilliary program for printing out non-singlet CFF's
C
C    $Id$
C     *******

C     ****p* auxns.f/AUXNS
C  NAME
C    AUXNS  --  Auxilliary program for printing out non-singlet CFF's
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


      PROGRAM AUXNS

      IMPLICIT NONE
      INTEGER CZ
      DOUBLE COMPLEX TMP
      DOUBLE PRECISION Q02, Q2EXP, MP
      PARAMETER (MP = 0.938272d0 )

      INCLUDE 'header.f'

      CALL READPAR
      NF = 4
      ANSATZ = 'NSFIT'
      CZ = 1

*       1  Q02       
!       PAR(1) =   2.5d0    
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
* ----  21 NG  (irrelevant for NS !)-----
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
* ----  31 NU  ( irrelevant, see D)-
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
* ----  41 ND   (fakes whole val.) -
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


      Q02 = 2.5d0
      Q2EXP = 4.0d0
      XI = 0.04d0
      DEL2 = -0.25d0
      SCHEME = 'MSBAR'
*   sea or not?
      PAR(11) = 0.0d0
      PAR(11) = 4.0d0 / 15.0d0

      WRITE (*,901) ANSATZ, DEL2, XI
      WRITE (*,903) Q02, Q2EXP

      PAR(1) = Q02
      Q2 = Q02
      P = 0
      CALL INIT
      CALL CFFF 
      WRITE (*,902) "LO (no evol)", CFF(0)

      Q2 = Q2EXP
      CALL INIT
      CALL CFFF 
      WRITE (*,902) "LO (LO evol)", CFF(0)

      P = 1
      CZERO = CZ
      Q2 = Q02
      CALL INIT
      CALL CFFF 
      WRITE (*,902) "MS NLO (no evol)", CFF(1)

      SCHEME = 'MSBLO'
      Q2 = Q2EXP
      CALL INIT
      CALL CFFF 
      WRITE (*,902) "MS NLO (LO evol)", CFF(1)
*      TMP = CFF(1)

      SCHEME = 'MSBAR'
      CALL INIT
      CALL CFFF 
*     WRITE (*,902) "MS NLO (evol D)", CFF(1) - TMP
      WRITE (*,902) "MS NLO (evol D)", CFF(1)

      SCHEME = 'MSBND'
      CALL INIT
      CALL CFFF 
*     WRITE (*,902) "MS NLO (evol ND)", CFF(1) - TMP
      WRITE (*,902) "MS NLO (evol ND)", CFF(1)

      SCHEME = 'CSBAR'
      Q2 = Q02
      CALL INIT
      CALL CFFF 
      WRITE (*,902) "CS NLO (no evol)", CFF(1)

      Q2 = Q2EXP
      CALL INIT
      CALL CFFF 
      WRITE (*,902) "CS NLO (evol)", CFF(1)

 901  FORMAT (1X, /, 1X, 'ANSATZ = ', A, 4X, 'DELTA^2 = ', F5.2,
     &    4X,'XI = ', F8.5, /)
 902  FORMAT (1X, A16, ' H = ', F8.2, ' + ', F8.2, ' I')
 903  FORMAT (1X, 'evolution : ', F4.1, ' --> ', F5.1, ' GeV^2', /)
      STOP
      END
C     ****
