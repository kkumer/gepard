C     ****h* gepard/auxns.f
C  FILE DESCRIPTION
C    Auxilliary program for printing out non-singlet CFF's
C
C    $Id$
C     *******

C     ****p* auxns.f/AUXNS
C  NAME
C    AUXNS  --  Auxilliary program for printing out non-singlet CFF's
C  DESCRIPTION
C             Calculates LO and NLO CFFs in both MSBAR
C             and CSBAR schemes, with various choices
C             for evolution.
C
C   Expected result:
C      ANSATZ = NSFIT     DELTA^2 = -0.25    XI =  0.40000
C
C        evolution :  2.5 -->   4.0 GeV^2
C
C        LO (no evol) H =    -0.11 +     0.27 I
C        LO (LO evol) H =    -0.09 +     0.26 I
C    MS NLO (no evol) H =    -0.17 +     0.24 I
C    MS NLO (LO evol) H =    -0.14 +     0.22 I
C     MS NLO (evol D) H =     0.00 +     0.00 I
C    MS NLO (evol ND) H =     0.00 +     0.00 I
C    CS NLO (no evol) H =    -0.17 +     0.24 I
C       CS NLO (evol) H =    -0.15 +     0.24 I
C
C  SYNOPSIS

      PROGRAM AUXNS

C
C  CHILDREN
C      READPAR, INIT, CFF
C
C  SOURCE
C
      IMPLICIT NONE
      INTEGER CZ
      DOUBLE COMPLEX TMP
      DOUBLE PRECISION Q2EXP, Q02IN, MP
      PARAMETER (MP = 0.938272d0 )

#include "../header.F"

      PID = 1
      FFTYPE = 'NONSINGLET'

      CALL READPAR
      CZ = 1

#include "ansatz.F"
      ANSATZ = 'NSFIT'


      Q02IN = 2.5d0
      Q2EXP = 4.0d0
      XI = 0.4d0
      DEL2 = -0.25d0
      SCHEME = 'CSBAR'

*   sea or not?
      PAR(11) = 0.0d0
      PAR(11) = 4.0d0 / 15.0d0

      WRITE (*,901) ANSATZ, DEL2, XI
      WRITE (*,903) Q02, Q2EXP

      Q2 = Q2EXP
      Q02 = Q2EXP
      P = 0
      CALL INIT
      CALL CFFF 
      WRITE (*,902) "LO (no evol)", CFF(0)

      Q02 = Q02IN
      CALL INIT
      CALL CFFF 
      WRITE (*,902) "LO (LO evol)", CFF(0)

      P = 1
      CZERO = CZ
      Q02 = Q2EXP
      CALL INIT
      CALL CFFF 
      WRITE (*,902) "MS NLO (no evol)", CFF(1)

      SCHEME = 'MSBLO'
      Q02 = Q02IN
      CALL INIT
      CALL CFFF 
      WRITE (*,902) "MS NLO (LO evol)", CFF(1)
      TMP = CFF(1)

      SCHEME = 'MSBDI'
      CALL INIT
      CALL CFFF 
      WRITE (*,902) "MS NLO (evol D)", CFF(1) - TMP
      TMP = CFF(1)

      SCHEME = 'MSBAR'
      CALL INIT
      CALL CFFF 
      WRITE (*,902) "MS NLO (evol ND)", CFF(1) - TMP

      SCHEME = 'CSBAR'
      Q02 = Q2EXP
      CALL INIT
      CALL CFFF 
      WRITE (*,902) "CS NLO (no evol)", CFF(1)

      Q02 = Q02IN
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
