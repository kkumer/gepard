C     ****h* gepard/fitres.f
C  FILE DESCRIPTION
C    calculation of Q^2 dependence of {\cal H} singlet
C    DVCS form factor
C    Produces data files for Fig. 8. of  KMPKS06b paper
C
C    $Id$
C     *******


C     ****p* fitres.f/FITRES
C  NAME
C     FITRES  --  Ploting final fit results
C  DESCRIPTION
C    Produces data files for Fig. ? of  KMPKS06b
C  OUTPUT
C       fitres[0-3].dat  --  4 files (one for each panel of Figure
C                            fitres) with 1-3 sets of point coordinates
C                            corresponding to 1-3 lines on a graph
C  IDENTIFIERS
C          NPOINTS -- number of points on each line
C              XS  -- array(NPOINTS), values of x-coordinates (-tt, Q2, W)
C           POINTS -- array(4, 3, NPOINTS), holds results for
C                 four panels (Fig ?, in order  a, b, c, d)
C
C                P -- approximation order N^{P}LO P=0,1,2
C
C  CHILDREN
C      READPAR, INIT, CFFF, DCARG
C  SOURCE
C

      PROGRAM FITRES

      IMPLICIT NONE
      INTEGER PT, NPOINTS, LN, PANEL
      DOUBLE PRECISION MT, MP, LSTART, LEND, LSTEP
      PARAMETER ( NPOINTS = 10 )
      DOUBLE PRECISION POINTS(0:3, 4, NPOINTS)
      DOUBLE PRECISION W, PARSIGMA, SIGMA
      CHARACTER FNAME*30
*     Proton mass:
      PARAMETER (MP = 0.938272d0 )

      INCLUDE '../header.f'

      PROCESS = 'DVCS'
      FFTYPE = 'SINGLET'

      CALL READPAR
    
      INCLUDE 'ansatz.f'

* ----  "RESULTING FIT" ----------------

      PAR(1)    =    7.0d0
      PAR(2)    =    0.0536d0


      PAR(11)   =    0.183
      PAR(21)   =    0.531
      PAR(12)   =     1.18
      PAR(22)   =     1.20
      PAR(14)   =     1.39
      PAR(24)   =     1.21

*     NLO_DVCS

      PAR(11)   =    0.136
      PAR(21)   =    0.350
      PAR(12)   =     1.26
      PAR(22)   =     1.31
      PAR(14)   =     0.993
      PAR(24)   =     0.599

*   Removing valence quarks
      PAR(31)   =     0.0d0
      PAR(41)   =     0.0d0

*     Files that will hold results

      OPEN (UNIT = 10, FILE = "fitres0.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 11, FILE = "fitres1.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 12, FILE = "fitres2.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 13, FILE = "fitres3.dat", STATUS = "UNKNOWN")

      DO 5 PANEL = 0, 3
  5         WRITE (10 + PANEL, *) '# Output of fitres.f. See prolog of 
     & that program'


*   ------  PANEL 0  -------------------
*     Looping over two choices for Q2
      DO 20 LN =1, 2

        IF (LN .EQ. 1) THEN
          Q2 = 4.0d0 
          W = 71.d0
          XI = Q2 / ( 2.0d0 * W**2 + Q2)
          FNAME = 'DATA/DVCS_H1_1.DAT'
        ELSE
          Q2 = 8.0d0 
          W = 82.d0
          XI = Q2 / ( 2.0d0 * W**2 + Q2)
          FNAME = 'DATA/DVCS_H1_2.DAT'
        END IF
        CALL PRINTDATA(10, FNAME)

        CALL INIT
        
*     Looping over MT
        LSTART = 0.05
        LEND = 0.9
        LSTEP = (LEND - LSTART) / (NPOINTS - 1)

        DO 10 PT = 1, NPOINTS

          MT = LSTART + (PT - 1) * LSTEP
          WRITE (UNIT=10,FMT=998) MT, PARSIGMA(MT), 0.0

 10     CONTINUE
        WRITE (UNIT=10, FMT=999)
 20   CONTINUE

*   ------  PANEL 1  -------------------
*     Looping over two choices for W
      DO 30 LN =1, 2

        IF (LN .EQ. 1) THEN
          W = 82.d0
          FNAME = 'DATA/DVCS_H1_3.DAT'
        ELSE
          W = 89.d0
          FNAME = 'DATA/DVCS_ZEUS_1.DAT'
        END IF
        CALL PRINTDATA(11, FNAME)

        CALL INIT
        
*     Looping over Q2
        LSTART = 3.0d0
        LEND = 90.0d0
        LSTEP = (LEND - LSTART) / (NPOINTS - 1)

        DO 40 PT = 1, NPOINTS

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
          Q2 = 4.d0
          FNAME = 'DATA/DVCS_H1_5.DAT'
        ELSE IF (LN .EQ. 2) THEN
          Q2 = 8.d0
          FNAME = 'DATA/DVCS_H1_6.DAT'
        ELSE
          Q2 = 9.6d0
          FNAME = 'DATA/DVCS_ZEUS_2.DAT'
        END IF
        CALL PRINTDATA(12, FNAME)

        CALL INIT
        
*     Looping over W
        LSTART = 30.d0
        LEND = 140.d0
        LSTEP = (LEND - LSTART) / (NPOINTS - 1)

        DO 60 PT = 1, NPOINTS

          W = LSTART + (PT - 1) * LSTEP
          XI = Q2 / ( 2.0d0 * W**2 + Q2)
          WRITE (UNIT=12,FMT=998) W, SIGMA(), 0.0

 60     CONTINUE
        WRITE (UNIT=12, FMT=999)
 50   CONTINUE
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
