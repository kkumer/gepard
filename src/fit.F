C     ****h* gepard/fit.f
C  FILE DESCRIPTION
C    Fitting of DVCS (and DIS in future) experimental data.
C
C    $Id$
C     *******

C     ****p* fit.f/FIT
C  NAME
C    FIT   --  Determines parameters of GPDs (and low-energy
C              input point Q0), by fitting to DVCS and DIS
C              experimental data.
C            
C  DESCRIPTION
C             Just a driver. Calls minuit subroutine
C             which does everything else.
C
C  INPUTS
C          FIT.CMD -- file with specification of fit parameters
C                     and with minuit commands.
C
C  OUTPUT
C          Goes to standard output. Good way to both observe
C          progress and save results is to use:
C            time ./fit | tee -a savefile
C
C  CHILDREN
C      MINUIT, FCN (via MINUIT)
C
C  SOURCE
C


      PROGRAM FIT

      IMPLICIT NONE
      INTEGER ISINTER
      CHARACTER SCHEME*5, ANSATZ*6
      CHARACTER COMFILE*10
      EXTERNAL FCN
      DATA COMFILE /'MINUIT.CMD'/

      COMMON /INTERACTIVE/ ISINTER
      COMMON / LABELS    /  SCHEME, ANSATZ

*   Reading from first line do we want interactive Minuit session.
      OPEN (UNIT = 17, FILE = COMFILE, STATUS = 'OLD')
      READ (17, *) ISINTER
      CLOSE (17)

      SCHEME = 'CSBAR'
      CALL READPAR
      CALL INIT

      IF (ISINTER .NE. 1) THEN
        OPEN (UNIT = 5, FILE = COMFILE, STATUS = 'OLD')
      ELSE
        WRITE (*, 801) COMFILE
      END IF

      CALL MINUIT(FCN,0)  

801   FORMAT (1X, / 1X, 72('-') / 1X, "Interactive mode" / 1X,
     &  "Answer 'SET INPUT 18 ", A, "' to the next querry." / 1X,
     &  "(Rewinding the unit is irrelevant, I guess.)" / 1X, 72('-') //)
      STOP
      END
C     ****

C     ****s* fit.f/FCN
C  NAME
C     FCN  --  minimization subroutine called by minuit
C            
C  DESCRIPTION
C    Reads experimental datasets and calculates chi-square by
C    comparing to theoretical prediction functions.
C    Writes out theoretical line points for fit plotting.
C  SYNOPSIS
C     SUBROUTINE FCN(NPAR,GIN,F,A,IFLAG,FUTIL)
C
C     IMPLICIT NONE
C     INTEGER NPAR, IFLAG
C     DOUBLE PRECISION F, FUTIL
C     DOUBLE PRECISION A(*),GIN(*)
C    
C  INPUTS 
C            DATASET?  --  ? = 1, ..., 9, 0.  (name of 8 chars exactly!!)
C                          files with experimental data formatted
C                          like:   x  y  dy_stat  dy_syst
C  OUTPUT
C             FIT.PLT  --  File with fit results formatted like:
C                          x  y(x)_theor
C                          (with empty lines separating graphs)
C  IDENTIFIERS
C          N  --  Number of points in a given data set
C          X  --  Array with x-values of experimental points
C          Y  --  Array with y-values of experimental points
C         DY  --  Array with errors of y-values
C          A  --  Array with fit parameters
C     FITPAR  --  Equal to A
C          W2 -- (gamma*-proton) invariant mass squared
C          XI -- DVCS scaling parameter or Bjorken x
C        DEL2 -- DVCS asymmetry parameter (P2-P1)^2
C          Q2 -- photon virtuality squared
C
C  CHILDREN
C      PARSIGMA, SIGMA, F2F
C
C  PARENTS
C      MINUIT, FCN (via MINUIT)
C
C  SOURCE
C


      SUBROUTINE FCN(NPAR,GIN,CHISQ,A,IFLAG,FUTIL)

      IMPLICIT NONE
      INTEGER NPAR, IFLAG, NPARMAX
      PARAMETER (NPARMAX = 10)
      DOUBLE PRECISION CHISQ, FUTIL
      DOUBLE PRECISION A(NPARMAX),GIN(NPARMAX)
      INTEGER K, IER, XPANELS, YPANELS
      DOUBLE PRECISION FITPAR(NPARMAX)
      DOUBLE PRECISION CHISQPART
      CHARACTER DATAFNAME*15
      CHARACTER OUTFILE*15
#ifndef NOPGPLOT
      INTEGER PGBEG
#endif

      COMMON / FITPARAMS  /  FITPAR

*     Clear output buffer

      CALL FLUSHOUT()
      
*     Transfer parameters to common block

      DO 10 K = 1, NPARMAX
 10     FITPAR(K) = A(K)

*     Initialization of printing the results

      IF (IFLAG .EQ. 3)  THEN
*       File for writing out points of final fit
        OPEN (UNIT = 21, FILE = 'FIT.OUT', STATUS = 'UNKNOWN')
      END IF

      CHISQ = 0.0d0

      OPEN (UNIT = 11, FILE = 'FIT.INI', STATUS = 'OLD')
      READ (11, *) OUTFILE
      READ (11, *) XPANELS, YPANELS
#ifndef NOPGPLOT
*    Initailization of plotting on 2x2 grid of panels
      IF (IFLAG .EQ. 3)  THEN
        IER = PGBEG(0, OUTFILE, XPANELS, YPANELS)
        IF (IER.NE.1) STOP
        CALL PGSCH(1.8)
      END IF
#endif

*     Process only files specified in FIT.INI between
*      'START' and 'STOP'

 15   READ (11, *) DATAFNAME
      IF (DATAFNAME(:5) .NE. 'START') GOTO 15
 18   READ (11, *) DATAFNAME
      IF (DATAFNAME(:4) .EQ. 'STOP') THEN
        GOTO 20
      ELSE
        CALL PROCDATA (DATAFNAME, IFLAG, CHISQPART)
        CHISQ = CHISQ + CHISQPART
        IF (IFLAG .EQ. 3) WRITE (21, *)
      END IF
      GOTO 18
 20   CLOSE(11)

      IF (IFLAG .EQ. 3)  THEN
         WRITE (*,*) 'CHISQ = ', CHISQ
         CLOSE (21)
#ifndef NOPGPLOT
         CALL PGEND
#endif
      END IF


      RETURN
      END
C     ****


C     ****s* fit.f/PROCDATA
C  NAME
C    PROCDATA   --  Process a single data set
C  DESCRIPTION
C    Reads experimental datasets and calculates chi-square by
C    comparing to theoretical prediction functions.
C    Writes out theoretical line points for fit plotting.
C  SYNOPSIS
C     SUBROUTINE PROCDATA (FNAME, IFLAG, CHISQPART)
C
C     CHARACTER  FNAME*15
C     INTEGER IFLAG
C     DOUBLE PRECISION CHISQPART
C  INPUTS
C      FNAME  --  Name of the file containing data
C      IFLAG  --  If .EQ. 3 then print results out, see minuit doc
C  OUTPUT
C   CHISQPART --  contribution to chi-square from processed data set
C  SOURCE
C
      SUBROUTINE PROCDATA (FNAME, IFLAG, CHISQPART)

      IMPLICIT NONE
      CHARACTER  FNAME*15
      INTEGER IFLAG
      DOUBLE PRECISION CHISQPART
      CHARACTER YOBS*8
      INTEGER N, K
      DOUBLE PRECISION X, Y, DY, THY, DIFSG
      DOUBLE PRECISION STAT, SYS
      DOUBLE PRECISION XI, DEL2, Q2, Q02,  W
      DOUBLE PRECISION WIN, Q2IN
      DOUBLE PRECISION XBJ, XBJIN
      DOUBLE PRECISION PARSIGMA, SIGMA
      DOUBLE PRECISION F2(0:2)
#ifndef NOPGPLOT
      LOGICAL DRAW
      REAL LASTX, LASTY
#endif

      COMMON / KINEMATICS /  XI, DEL2, Q2, Q02
      COMMON / F2P        /  F2

      CHISQPART = 0.d0

      IF (IFLAG .EQ. 3) THEN
        WRITE (21, 902) FNAME
#ifndef NOPGPLOT
        CALL PGPAGE
        CALL PGVSTD
        DRAW = .FALSE.
#endif
      END IF

      OPEN (UNIT = 12, FILE = FNAME, STATUS = 'OLD')

      READ (12,*) YOBS

      IF (YOBS(:2) .EQ. 'PA') THEN

        READ (12,*) W
        READ (12,*) Q2
        READ (12,*) N

        XI = Q2 / ( 2.0d0 * W**2 + Q2)
#ifndef NOPGPLOT
        IF (IFLAG .EQ. 3) THEN
          CALL PGSWIN (0., 1., -1.3, 2.)
          CALL PGBOX ('BCNST1', 0.0, 0, 'BCLNST', 0.0, 0)
          CALL PGLAB("-t", 'd\\gs/dt', FNAME)
        END IF
#endif

        DO 110 K = 1, N
        READ (12, *) X, Y, STAT, SYS
        DY = SQRT( STAT*STAT + SYS*SYS )
        THY = PARSIGMA(X)
        CHISQPART = CHISQPART + ( (THY - Y)**2 / DY**2 )
        IF (IFLAG .EQ. 3)  THEN
          DIFSG = (THY - Y) / DY
          WRITE (21, 901) X, THY, Y, DY, DIFSG
#ifndef NOPGPLOT
          CALL PGERRY(1, SNGL(X), LOG10(SNGL(Y-DY)), LOG10(SNGL(Y+DY)),
     &           3.0)
          CALL PGSCI(2)
          IF (DRAW) THEN
            CALL PGMOVE(LASTX, LASTY)
            CALL PGDRAW(SNGL(X), LOG10(SNGL(THY)))
          END IF
          CALL PGPT1(SNGL(X), LOG10(SNGL(THY)), 17)
          LASTX = SNGL(X)
          LASTY = LOG10(SNGL(THY)) 
          DRAW = .TRUE.
          CALL PGSCI(1)
#endif
        END IF
110     CONTINUE

      ELSE IF (YOBS(:2) .EQ. 'SI') THEN

        READ (12,*) WIN
        READ (12,*) Q2IN
        READ (12,*) N

        DO 120 K = 1, N
        READ (12, *) X, Y, STAT, SYS
        DY = SQRT( STAT*STAT + SYS*SYS )
*   Which of W or Q2 is on x-axis?
        IF (WIN .LT. 0.) THEN
          W = X
          Q2 = Q2IN
#ifndef NOPGPLOT
          IF (IFLAG .EQ. 3) THEN
            CALL PGSWIN (30., 140., 0., 12.)
            CALL PGBOX ('BCNST1', 0.0, 0, 'BCNST', 0.0, 0)
            CALL PGLAB("W", '\\gs', FNAME)
          END IF
#endif
        ELSE IF (Q2IN .LT. 0) THEN
          Q2 = X
          W = WIN
#ifndef NOPGPLOT
          IF (IFLAG .EQ. 3) THEN
            CALL PGSWIN (0., 90., -1.7, 1.4)
            CALL PGBOX ('BCNST1', 0.0, 0, 'BCLNST', 0.0, 0)
            CALL PGLAB("Q\\u2\\d", '\\gs', FNAME)
          END IF
#endif
        ELSE
          CALL ERROR ('GeParD', 'PROCDATA',
     &    'Either W or Q2 in ' // FNAME // ' should be negative!',
     &    4, 2)
        END IF
        XI = Q2 / (2.0d0 * W**2 + Q2)
        THY = SIGMA()
        CHISQPART = CHISQPART + ( (THY - Y)**2 / DY**2 )
        IF (IFLAG .EQ. 3)  THEN
          DIFSG = (THY - Y) / DY
          WRITE (21, 901) X, THY, Y, DY, DIFSG
#ifndef NOPGPLOT
          IF (WIN .LT. 0.) THEN
            CALL PGERRY(1, SNGL(X), SNGL(Y-DY), SNGL(Y+DY), 3.0)
            CALL PGSCI(2)
            IF (DRAW) THEN
              CALL PGMOVE(LASTX, LASTY)
              CALL PGDRAW(SNGL(X), SNGL(THY))
            END IF
            CALL PGPT1(SNGL(X), SNGL(THY), 17)
            LASTX = SNGL(X)
            LASTY = SNGL(THY)
            DRAW = .TRUE.
          ELSE
            CALL PGERRY(1,SNGL(X), LOG10(SNGL(Y-DY)), LOG10(SNGL(Y+DY)),
     &             3.0)
            CALL PGSCI(2)
            IF (DRAW) THEN
              CALL PGMOVE(LASTX, LASTY)
              CALL PGDRAW(SNGL(X), LOG10(SNGL(THY)))
            END IF
            CALL PGPT1(SNGL(X), LOG10(SNGL(THY)), 17)
            LASTX = SNGL(X)
            LASTY = LOG10(SNGL(THY)) 
            DRAW = .TRUE.
          END IF
          CALL PGSCI(1)
#endif
        END IF
120     CONTINUE

      ELSE IF (YOBS(:2) .EQ. 'F2') THEN

        READ (12,*) XBJIN
        READ (12,*) Q2IN
        READ (12,*) N

        DO 130 K = 1, N
        READ (12, *) X, Y, STAT, SYS
        DY = SQRT( STAT*STAT + SYS*SYS )
*   Which of X_BJ or Q2 is on x-axis?
        IF (XBJIN .LT. 0.) THEN
          XBJ = X
          Q2 = Q2IN
#ifndef NOPGPLOT
          IF (IFLAG .EQ. 3) THEN
            CALL PGSWIN (0., 0.01, 0., 1.5)
            CALL PGBOX ('BCNST1', 0.0, 0, 'BCNST', 0.0, 0)
            CALL PGLAB("x\\dBJ\\u", 'F\\d2\\u(x\\dBJ\\u)', FNAME)
          END IF
#endif
        ELSE IF (Q2IN .LT. 0) THEN
          Q2 = X
          XBJ = XBJIN
#ifndef NOPGPLOT
          IF (IFLAG .EQ. 3) THEN
            CALL PGSWIN (0., 70., 0., 1.5)
            CALL PGBOX ('BCNST1', 0.0, 0, 'BCNST', 0.0, 0)
            CALL PGLAB("Q\\u2\\d", 'F\\d2\\u(Q\\u2\\d)', FNAME)
          END IF
#endif
        ELSE
          CALL ERROR ('GeParD', 'PROCDATA',
     &    'Either XBJ or Q2 in ' // FNAME // ' should be negative!',
     &    4, 2)
        END IF
        XI = XBJ
        CALL F2F
        THY = F2(1)
        CHISQPART = CHISQPART + ( (THY - Y)**2 / DY**2 )
        IF (IFLAG .EQ. 3)  THEN
          DIFSG = (THY - Y) / DY
          WRITE (21, 901) X, THY, Y, DY, DIFSG
#ifndef NOPGPLOT
          CALL PGERR1(6, SNGL(X), SNGL(Y), SNGL(DY), 3.0)
          CALL PGSCI(2)
          IF (DRAW) THEN
            CALL PGMOVE(LASTX, LASTY)
            CALL PGDRAW(SNGL(X), SNGL(THY))
          END IF
          CALL PGPT1(SNGL(X), SNGL(THY), 17)
          LASTX = SNGL(X)
          LASTY = SNGL(THY) 
          DRAW = .TRUE.
          CALL PGSCI(1)
#endif
        END IF
130     CONTINUE
      ELSE

        CALL ERROR ('GeParD', 'PROCDATA',
     &  'Record ' // YOBS // ' in ' // FNAME // ' unrecognized',
     &  3, 2)

      END IF

      CLOSE (12)
901   FORMAT (1X, F10.6, 5X, F9.4, 5X, F9.4, 3X, F7.2, 7X, F5.1)
902   FORMAT (1X, 72('-') / 1X, 35(' '), A, 21(' '), / 
     &  6X, "X", 11X, "Y(X)_FIT", 6X, "Y(X)_EXP",
     &  3X, "DY_EXP", 4X, "(FIT-EXP)/DY_EXP" / 1X, 72('-'))

      RETURN
      END
C     ****

