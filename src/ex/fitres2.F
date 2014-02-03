C     ****h* gepard/fitres2.f
C  FILE DESCRIPTION
C    calculation of Q^2 dependence of {\cal H} singlet
C    DVCS form factor
C    Produces data files for Fig. 1. of  smallx paper
C
C    $Id$
C     *******


C     ****p* fitres2.f/FITRES2
C  NAME
C     FITRES2  --  Ploting final fit results
C  SYNOPSIS

      PROGRAM FITRES2

C  DESCRIPTION
C    Produces data files for Fig. 1 of smallx paper
C  OUTPUT
C       fitres2[0-3].dat  --  4 files (one for each panel of Figure
C                            fitres2) with 1-3 sets of point coordinates
C                            corresponding to 1-3 lines on a graph
C  IDENTIFIERS
C          LNPOINTS -- number of points on each line
C              XS  -- array(LNPOINTS), values of x-coordinates (-tt, Q2, W)
C
C                P -- approximation order N^{P}LO P=0,1,2
C
C  CHILDREN
C      READPAR, INIT, CFFF, DCARG
C  SOURCE
C
      IMPLICIT NONE
      INTEGER PT, LNPOINTS, LN, PANEL, K
      DOUBLE PRECISION MT, MP, LSTART, LEND, LSTEP
      PARAMETER ( LNPOINTS = 40 )
      DOUBLE PRECISION W, PARSIGMA, SIGMA
      CHARACTER FNAME*30
*     Proton mass:
      PARAMETER (MP = 0.938272d0 )

#include "../header.F"

      PROCESS = PID

      CALL READPAR

      ANSATZ = 'FIT'
      Q02    =    4.0d0
    
#include "ansatz.F"

* ----  "RESULTING FIT" ----------------


*    NNLO 
*      P = 2
*      SCHEME = 'CSBAR'
*
*      PAR(11)   =    0.1673
*      PAR(12)   =    1.14
*      PAR(14)   =    1.167
*      PAR(15)   =    0.0
*
*      PAR(21)   =    0.571
*      PAR(22)   =    1.068
*      PAR(24)   =    1.391
*      PAR(25)   =    0.0

*    NLO 
*      P = 1
*      SCHEME = 'CSBAR'
*
*      PAR(11)   =    0.167
*      PAR(12)   =     1.14
*      PAR(14)   =     1.34
*      PAR(15)   =     0.0
*
*      PAR(21)   =    0.535
*      PAR(22)   =     1.09
*      PAR(24)   =     1.59
*      PAR(25)   =     0.0

*    LO 
      P = 0
      SCHEME = 'CSBAR'

      PAR(11)   =    0.1520
      PAR(12)   =    1.157506
      PAR(13)   =    0.15
      PAR(14)   =    0.35612
      PAR(15)   =    0.0
      PAR(16)   =    2.0
      PAR(18)   =    0.0
      PAR(19)   =   -9.8177889

      PAR(22)   =    1.2473
      PAR(23)   =    0.15
      PAR(24)   =    0.7
      PAR(25)   =    0.0
      PAR(26)   =    2.0
      PAR(28)   =    0.0
      PAR(29)   =  -32.401


*     Files that will hold results

      OPEN (UNIT = 10, FILE = "fitres20.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 11, FILE = "fitres21.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 12, FILE = "fitres22.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 13, FILE = "fitres23.dat", STATUS = "UNKNOWN")

      DO 5 PANEL = 0, 3
  5         WRITE (10 + PANEL, *) '# Output of fitres2.f. See prolog of 
     & that program'


*   ------  PANEL 0  -------------------
*     Looping over three choices for Q2
      DO 20 LN =1, 3

        IF (LN .EQ. 1) THEN
          Q2 = 8.0d0 
          W = 82.d0
          XI = Q2 / ( 2.0d0 * W**2 + Q2)
          FNAME = 'DATA/DVCS_H1_0506_06.DAT'
        ELSE IF (LN .EQ. 2) THEN
          Q2 = 15.5d0 
          W = 82.d0
          XI = Q2 / ( 2.0d0 * W**2 + Q2)
          FNAME = 'DATA/DVCS_H1_0506_07.DAT'
        ELSE
          Q2 = 25.0d0 
          W = 82.d0
          XI = Q2 / ( 2.0d0 * W**2 + Q2)
          FNAME = 'DATA/DVCS_H1_0506_08.DAT'
        END IF
        CALL PRINTDATA(10, FNAME)

        CALL INIT
        
*     Looping over MT
        LSTART = 0.05
        LEND = 0.9
        LSTEP = (LEND - LSTART) / (LNPOINTS - 1)

        DO 10 PT = 1, LNPOINTS

          MT = LSTART + (PT - 1) * LSTEP
          MTIND = NMTS + NMTSEXP + 1
          MTS(MTIND) = MT
          WRITE (UNIT=10,FMT=998) MT, PARSIGMA(), 0.0

 10     CONTINUE
        WRITE (UNIT=10, FMT=999)
 20   CONTINUE

*   ------  PANEL 1  -------------------
*     Looping over choices for W
      DO 30 LN =1, 3

        IF (LN .EQ. 1) THEN
          W = 82.d0
          FNAME = 'DATA/DVCS_H1_3.DAT'
        ELSE IF (LN .EQ. 2) THEN
          W = 82.d0
          FNAME = 'DATA/DVCS_H1_0506_01.DAT'
        ELSE
          W = 89.d0
          FNAME = 'DATA/DVCS_ZEUS_1.DAT'
        END IF
        CALL PRINTDATA(11, FNAME)


        CALL INIT
        
*     Looping over Q2
        LSTART = 3.0d0
        LEND = 90.0d0
        LSTEP = (LEND - LSTART) / (LNPOINTS - 1)

        DO 40 PT = 1, LNPOINTS

          Q2 = LSTART + (PT - 1) * LSTEP
          XI = Q2 / ( 2.0d0 * W**2 + Q2)
          WRITE (UNIT=11,FMT=998) Q2, SIGMA(), 0.0

 40     CONTINUE
        WRITE (UNIT=11, FMT=999)
 30   CONTINUE

*   ------  PANEL 2  -------------------
*     Looping over three choices for Q2
      DO 50 LN =1, 3

        IF (LN .EQ. 1) THEN
          Q2 = 8.d0
          FNAME = 'DATA/DVCS_H1_0506_03.DAT'
        ELSE IF (LN .EQ. 2) THEN
          Q2 = 15.5d0
          FNAME = 'DATA/DVCS_H1_0506_04.DAT'
        ELSE
          Q2 = 25.d0
          FNAME = 'DATA/DVCS_H1_0506_05.DAT'
        END IF
        CALL PRINTDATA(12, FNAME)

        CALL INIT
        
*     Looping over W
        LSTART = 30.d0
        LEND = 140.d0
        LSTEP = (LEND - LSTART) / (LNPOINTS - 1)

        DO 60 PT = 1, LNPOINTS

          W = LSTART + (PT - 1) * LSTEP
          XI = Q2 / ( 2.0d0 * W**2 + Q2)
          WRITE (UNIT=12,FMT=998) W, SIGMA(), 0.0

 60     CONTINUE
        WRITE (UNIT=12, FMT=999)
 50   CONTINUE


*   ------  PANEL 3  -------------------
*     Looping over choices for X_BJ
      PROCESS = 'DIS'
      DO 70 LN =1, 4
        IF (LN .EQ. 1) THEN
          XI = 0.008d0
          FNAME = 'DATA/DIS2H1.DAT'
        ELSE IF (LN .EQ. 2) THEN
          XI = 0.0032d0
          FNAME = 'DATA/DIS4H1.DAT'
        ELSE IF (LN .EQ. 3) THEN
          XI = 0.0013d0
          FNAME = 'DATA/DIS6H1.DAT'
        ELSE
          XI = 0.0005d0
          FNAME = 'DATA/DIS8H1.DAT'
        END IF
        CALL PRINTDATA(13, FNAME)

        DEL2 = 0.0d0
        CALL INIT
        CALL GETMBGPD
        
*     Looping over Q2
        LSTART = 2.d0
        LEND = 95.d0
        LSTEP = (LEND - LSTART) / (LNPOINTS - 1)

        DO 80 PT = 1, LNPOINTS

          Q2 = LSTART + (PT - 1) * LSTEP
          CALL F2F
          WRITE (UNIT=13,FMT=998) Q2, F2(P), 0.0

 80     CONTINUE
        WRITE (UNIT=13, FMT=999)
 70   CONTINUE


998   FORMAT (F8.3,5X,F12.7,5X,F3.1)
999   FORMAT (1X)

      STOP
      END
C     ****


      SUBROUTINE PRINTDATA (UNT, FNAME)

      IMPLICIT NONE
      INTEGER UNT
      CHARACTER  FNAME*30
      INTEGER K, NN
      DOUBLE PRECISION X(20), YY(20), DY(20)
      DOUBLE PRECISION STAT, SYS, AUX
      CHARACTER  YOBS*8


      OPEN (UNIT = 19, FILE = FNAME, STATUS = 'OLD')

      READ (19,*) YOBS
      READ (19,*) AUX
      READ (19,*) AUX
      READ (19,*) NN

      DO 10 K = 1, NN
        READ (19, *) X(K), YY(K), STAT, SYS
        DY(K) = SQRT( STAT*STAT + SYS*SYS )
        WRITE (UNIT=UNT,FMT=898) X(K), YY(K), DY(K)
 10   CONTINUE
      WRITE (UNIT=UNT,FMT=899)

898   FORMAT (1X,3F9.2)
899   FORMAT (1X)
      RETURN
      END
