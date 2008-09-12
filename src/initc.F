C     ****h* gepard/initc.F
C  FILE DESCRIPTION
C    Initialization of "evolved" Wilson coefficients
C
C    $Id$
C     *******

C     ****s* initc.F/INITC
C  NAME
C     INITC  --  initialize "evolved" Wilson coefficients
C            
C  DESCRIPTION
C    Reads experimental datasets and calculates evolved
C    Wilson coefficients for all needed Q2 values i.e.
C    Wilson coefficients multiplied by evolution operator
C  SYNOPSIS

      SUBROUTINE INITC

C  CHILDREN
C      LOOKUPQ
C  PARENTS
C      FIT, MNINTERFACE
C  SOURCE
C

      IMPLICIT NONE
      INTEGER NN, K, QIND, XPANELS, YPANELS, LOC
      CHARACTER*30 DATAFNAME, OUTFORMAT
      CHARACTER YOBS*8
      DOUBLE PRECISION X, YY, STAT, SYS
      DOUBLE PRECISION W, WIN, Q2IN, XBJIN
#include "header.F"


      NQS = 0
      NQSDIS = 0
      QS(1) = 0.0d0
      QSDIS(1) = 0.0d0

      OPEN (UNIT = 11, FILE = DATFILE, STATUS = 'OLD')
      READ (11, *) OUTFORMAT, XPANELS, YPANELS

*     Process only files specified in FIT.INI between
*      'START' and 'STOP'

 15   READ (11, *) DATAFNAME
      IF (DATAFNAME(:5) .NE. 'START') GOTO 15
 18   READ (11, *) DATAFNAME
      IF (DATAFNAME(:4) .EQ. 'STOP') THEN
        GOTO 20
      ELSE
      IF (DATAFNAME(:5) .EQ. 'CHISQ') THEN
*       Skip partial CHISQ
        READ (11, *) DATAFNAME
        IF (DATAFNAME(:4) .EQ. 'STOP') GOTO 20
      END IF
      IF (DATAFNAME(:6) .EQ. 'NEWFIG') THEN
*       Skip PGPLOT definitions
        READ (11, *) DATAFNAME
        READ (11, *) DATAFNAME
      END IF
C       BEGIN count datapoints from this file
      OPEN (UNIT = 12, FILE = DATAFNAME, STATUS = 'OLD')

      READ (12,*) YOBS

      IF (YOBS(:2) .EQ. 'PA') THEN

        PROCESS = 'DVCS'

        READ (12,*) W
        READ (12,*) Q2
        READ (12,*) NN

        
        CALL LOOKUPQ(LOC)
        IF ( LOC .EQ. 0 ) THEN
          NQS = NQS + 1
          QS(NQS) = Q2
        END IF

      ELSE IF (YOBS(:2) .EQ. 'SI') THEN

        PROCESS = 'DVCS'

        READ (12,*) WIN
        READ (12,*) Q2IN
        READ (12,*) NN

*   Which of W or Q2 is on x-axis?
        IF (WIN .LT. 0.) THEN
          Q2 = Q2IN
          CALL LOOKUPQ(LOC)
          IF ( LOC .EQ. 0 ) THEN
            NQS = NQS + 1
            QS(NQS) = Q2
          END IF
        ELSE IF (Q2IN .LT. 0) THEN
          DO 120 K = 1, NN
            READ (12, *) X, YY, STAT, SYS
            Q2 = X
            CALL LOOKUPQ(LOC)
            IF ( LOC .EQ. 0 ) THEN
              NQS = NQS + 1
              QS(NQS) = Q2
            END IF
 120      CONTINUE
        ELSE
          CALL ERROR ('GeParD', 'PROCDATA',
     &    'Either W or Q2 in ' // DATAFNAME // ' should be negative!',
     &    4, 2)
        END IF

      ELSE IF (YOBS(:2) .EQ. 'F2') THEN

        PROCESS = 'DIS'

        READ (12,*) XBJIN
        READ (12,*) Q2IN
        READ (12,*) NN

*   Which of X_BJ or Q2 is on x-axis?
        IF (XBJIN .LT. 0.) THEN
          Q2 = Q2IN
          CALL LOOKUPQ(LOC)
          IF ( LOC .EQ. 0 ) THEN
            NQSDIS = NQSDIS + 1
            QSDIS(NQSDIS) = Q2
          END IF
        ELSE IF (Q2IN .LT. 0) THEN
          DO 130 K = 1, NN
          READ (12, *) X, YY, STAT, SYS
            Q2 = X
            CALL LOOKUPQ(LOC)
            IF ( LOC .EQ. 0 ) THEN
              NQSDIS = NQSDIS + 1
              QSDIS(NQSDIS) = Q2
            END IF
 130      CONTINUE
        ELSE
          CALL ERROR ('GeParD', 'PROCDATA',
     &    'Either XBJ or Q2 in ' // DATAFNAME // ' should be negative!',
     &    4, 2)
        END IF
      ELSE

        CALL ERROR ('GeParD', 'PROCDATA',
     &  'Record ' // YOBS // ' in ' // DATAFNAME // ' unrecognized',
     &  3, 2)

      END IF

      CLOSE (12)

C       END count datapoints from this file
      END IF
      GOTO 18
 20   CLOSE(11)

C     Initialize grids with values of evolved C's
      PROCESS = 'DVCS'
      DO 310 QIND = 1, NQS
        Q2 = QS(QIND)
        CALL EVOLC(QIND)
 310  CONTINUE
      PROCESS = 'DIS'
      DO 320 QIND = 1, NQSDIS
        Q2 = QSDIS(QIND)
        CALL EVOLC(QIND)
 320  CONTINUE
        

      RETURN
      END
C     ****


      SUBROUTINE LOOKUPQ(LOC)

      INTEGER K, LOC
#include "header.F"

      LOC = 0
      IF ( PROCESS .EQ. 'DVCS' ) THEN
        DO 160 K = 1, NQS
           IF ( QS(K) .EQ. Q2 ) LOC = K
 160    CONTINUE
      ELSE
        DO 170 K = 1, NQSDIS
           IF ( QSDIS(K) .EQ. Q2 ) LOC = K
 170    CONTINUE
      END IF


      RETURN
      END
