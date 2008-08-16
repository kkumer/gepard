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

      PROCESS = 'DVCS'
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
      NQS = 1

*   sea or not?
      PAR(11) = 0.0d0
      PAR(11) = 4.0d0 / 15.0d0

      WRITE (*,901) ANSATZ, DEL2, XI
      WRITE (*,903) Q02, Q2EXP

      Q2 = Q2EXP
      Q02 = Q2EXP
      P = 0
      CALL INIT
      QS(1) = Q2
      CALL EVOLC(1)
      CALL CFFF 
      WRITE (*,902) "LO (no evol)", CFF(0)

      Q02 = Q02IN
      CALL INIT
      QS(1) = Q2
      CALL EVOLC(1)
      CALL CFFF 
      WRITE (*,902) "LO (LO evol)", CFF(0)

      P = 1
      CZERO = CZ
      Q02 = Q2EXP
      CALL INIT
      QS(1) = Q2
      CALL EVOLC(1)
      CALL CFFF 
      WRITE (*,902) "MS NLO (no evol)", CFF(1)

      SCHEME = 'MSBLO'
      Q02 = Q02IN
      CALL INIT
      QS(1) = Q2
      CALL EVOLC(1)
      CALL CFFF 
      WRITE (*,902) "MS NLO (LO evol)", CFF(1)
      TMP = CFF(1)

      SCHEME = 'MSBAR'
      CALL INIT
      QS(1) = Q2
      CALL EVOLC(1)
      CALL CFFF 
      WRITE (*,902) "MS NLO (evol D)", CFF(1) - TMP
      TMP = CFF(1)

      SCHEME = 'MSBND'
      CALL INIT
      QS(1) = Q2
      CALL EVOLC(1)
      CALL CFFF 
      WRITE (*,902) "MS NLO (evol ND)", CFF(1) - TMP

      SCHEME = 'CSBAR'
      Q02 = Q2EXP
      CALL INIT
      QS(1) = Q2
      CALL EVOLC(1)
      CALL CFFF 
      WRITE (*,902) "CS NLO (no evol)", CFF(1)

      Q02 = Q02IN
      CALL INIT
      QS(1) = Q2
      CALL EVOLC(1)
      CALL CFFF 
      WRITE (*,902) "CS NLO (evol)", CFF(1)

 901  FORMAT (1X, /, 1X, 'ANSATZ = ', A, 4X, 'DELTA^2 = ', F5.2,
     &    4X,'XI = ', F8.5, /)
 902  FORMAT (1X, A16, ' H = ', F8.2, ' + ', F8.2, ' I')
 903  FORMAT (1X, 'evolution : ', F4.1, ' --> ', F5.1, ' GeV^2', /)
      STOP
      END
C     ****
