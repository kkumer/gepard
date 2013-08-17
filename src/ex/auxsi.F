C     ****h* gepard/auxsi.F
C  FILE DESCRIPTION
C    Auxilliary program for printing out non-singlet CFF's
C
C    $Id$
C     *******

C     ****p* auxsi.F/AUXSI
C  NAME
C    AUXSI  --  Auxilliary program for printing out singlet CFF's
C  DESCRIPTION
C             Calculates LO and NLO CFFs in both MSBAR
C             and CSBAR schemes, with various choices
C             for evolution.
C
C    Expected result:
C
C    ANSATZ = FITBP     DELTA^2 = -0.25    XI =  0.00100
C   
C    evolution :  2.5 -->   4.0 GeV^2
C   
C        LO (no evol) H =      297.88 +     3027.22 I
C        LO (LO evol) H =      382.21 +     3314.32 I
C    MS NLO (no evol) H =       37.54 +     1366.25 I
C    MS NLO (LO evol) H =     -221.85 +     1105.77 I
C     MS NLO (evol D) H =       22.76 +     1510.30 I
C    MS NLO (evol ND) H =       84.80 +     1632.84 I
C    CS NLO (no evol) H =       96.47 +     1775.06 I
C       CS NLO (evol) H =      141.67 +     2013.85 I
C   
C  SYNOPSIS

      PROGRAM AUXSI

C  CHILDREN
C      READPAR, INIT, CFF
C  SOURCE
C

      IMPLICIT NONE
      INTEGER CZ
      DOUBLE COMPLEX TMP
      DOUBLE PRECISION Q2EXP, MP
      PARAMETER (MP = 0.938272d0 )
      CHARACTER SUBANSATZ*4
#include "../header.F"

      PROCESS = 'DVCS'
      FFTYPE = 'SINGLET'

      CALL READPAR

#include "ansatz.F"
      ANSATZ = 'FITBP'


      Q02 = 2.5d0
      Q2EXP = 4.0d0
      XI = 1.0d-3
      DEL2 = -0.25d0
      SCHEME = 'MSBDI'
      CZ = 1
      NQS = 1

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

!     WRITE (*,904) CND

      Q2 = Q02
      P = 0
      CALL INIT
      QS(1) = Q2
      CALL EVOLC(1,1)
      CALL CFFF 
      WRITE (*,902) "LO (no evol)", CFF(0)

      Q2 = Q2EXP
      CALL INIT
      QS(1) = Q2
      CALL EVOLC(1,1)
      CALL CFFF 
      WRITE (*,902) "LO (LO evol)", CFF(0)

      P = 1
      CZERO = CZ
      Q2 = Q02
      CALL INIT
      QS(1) = Q2
      CALL EVOLC(1,1)
      CALL CFFF 
      WRITE (*,902) "MS NLO (no evol)", CFF(1)

      SCHEME = 'MSBLO'
      Q2 = Q2EXP
      CALL INIT
      QS(1) = Q2
      CALL EVOLC(1,1)
      CALL CFFF 
      WRITE (*,902) "MS NLO (LO evol)", CFF(1)
      TMP = CFF(1)

      SCHEME = 'MSBDI'
      CALL INIT
      QS(1) = Q2
      CALL EVOLC(1,1)
      CALL CFFF 
!     WRITE (*,902) "MS NLO (evol D)", CFF(1) - TMP
      WRITE (*,902) "MS NLO (evol D)", CFF(1)
      TMP = CFF(1)

      SCHEME = 'MSBAR'
      CALL INIT
      QS(1) = Q2
      CALL EVOLC(1,1)
      CALL CFFF 
!     WRITE (*,902) "MS NLO (evol ND)", CFF(1) - TMP
      WRITE (*,902) "MS NLO (evol ND)", CFF(1)

      SCHEME = 'CSBAR'
      Q2 = Q02
      CALL INIT
      QS(1) = Q2
      CALL EVOLC(1,1)
      CALL CFFF 
      WRITE (*,902) "CS NLO (no evol)", CFF(1)

      Q2 = Q2EXP
      CALL INIT
      QS(1) = Q2
      CALL EVOLC(1,1)
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
