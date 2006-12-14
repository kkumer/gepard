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

      INCLUDE '../header.f'

      PROCESS = 'DVCS'
      FFTYPE = 'NONSINGLET'

      CALL READPAR
      CZ = 1

      INCLUDE 'ansatz.f'
      ANSATZ = 'NSFIT'


      Q02 = 2.5d0
      Q2EXP = 2.5d0
      XI = 0.4d0
      DEL2 = -0.25d0
      SCHEME = 'CSBAR'
*   sea or not?
      PAR(11) = 0.0d0
      PAR(11) = 4.0d0 / 15.0d0

      WRITE (*,901) ANSATZ, DEL2, XI
      WRITE (*,903) Q02, Q2EXP

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

 901  FORMAT (1X, /, 1X, 'ANSATZ = ', A, 4X, 'DELTA^2 = ', F5.2,
     &    4X,'XI = ', F8.5, /)
 902  FORMAT (1X, A16, ' H = ', F8.2, ' + ', F8.2, ' I')
 903  FORMAT (1X, 'evolution : ', F4.1, ' --> ', F5.1, ' GeV^2', /)
      STOP
      END
C     ****
