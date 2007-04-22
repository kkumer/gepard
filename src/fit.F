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
C             Calls minuit subroutine for minimization and
C             prints results in various ways.
C
C  INPUTS
C          FIT.INI -- file with names of files with experimental
C                     data to be fitted to
C       MINUIT.CMD -- file with specification of fit parameters
C                     and with minuit commands.
C
C  OUTPUT
C          Goes to files with extensions .out (tabular representation
C          of fit), .ps (graphical representation of fit), and .min
C          (Minuit output)
C
C  CHILDREN
C      READPAR, INIT, MINUIT, FCN (via MINUIT)
C
C  SOURCE
C


      PROGRAM FIT

      IMPLICIT NONE
      INTEGER ISINTER
      CHARACTER COMFILE*10, OUTFILE*20
      EXTERNAL FCN
      DATA COMFILE /'MINUIT.CMD'/
      INCLUDE 'header.f'

      COMMON /INTERACTIVE/ ISINTER

*   Reading from first line whether we want interactive Minuit session.
      OPEN (UNIT = 17, FILE = COMFILE, STATUS = 'OLD')
      READ (17, *) ISINTER
      CLOSE (17)

      CALL READPAR
      FFTYPE = 'SINGLET'

      PROCESS = 'DVCS'
      CALL INIT
      PROCESS = 'DIS'
      CALL INIT

      CALL INITC

      IF (ISINTER .NE. 1) THEN
*   We work in batch mode
*     File where Minuit batch commands are
        OPEN (UNIT = 5, FILE = COMFILE, STATUS = 'OLD')
        OPEN (UNIT = 11, FILE = 'FIT.INI', STATUS = 'OLD')
        READ (11, *) OUTFILE
        CLOSE (11)
*     File for writing out Minuit output
        OPEN (UNIT = 6, FILE = 
     &    OUTFILE(1:MAX(1,INDEX(OUTFILE//' ',' ')-1)) // '.min',
     &    STATUS = 'UNKNOWN')
      ELSE
*   We work in interactive mode
        WRITE (*, 801) COMFILE
      END IF

      CALL MINUIT(FCN,0)  

801   FORMAT (1X, // 1X, 72('-') / 1X, " -- Interactive mode ! --" / 1X,
     & "You can answer 'SET INPUT 18 ", A, "' to the querry below" / 1X,
     & "(Second question about rewinding the unit is irrelevant.)" / 1X,
     & 72('-') //)
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
C     INTEGER NPAR, IFLAG, NPARMAX
C     PARAMETER (NPARMAX = 10)
C     DOUBLE PRECISION CHISQ, FUTIL
C     DOUBLE PRECISION A(NPARMAX),GIN(NPARMAX)
C    
C  INPUTS 
C                NPAR  --  Number of variable parameters
C                   A  --  Values of all parameters
C               IFLAG  --  If .EQ. 3 then print results out, see minuit doc
C  OUTPUT
C                   F  --  Value of chi-square
C  IDENTIFIERS
C     FITPAR  --  Equal to A
C  CHISQPART  --  chi-square produced by single data set
C
C  CHILDREN
C      PROCDATA, PGPLOT routines
C
C  PARENTS
C      MINUIT, FIT (via MINUIT)
C
C  SOURCE
C


      SUBROUTINE FCN(NPAR,GIN,CHISQ,A,IFLAG,FUTIL)

      IMPLICIT NONE
      INTEGER NPAR, IFLAG
      DOUBLE PRECISION CHISQ, FUTIL
      INTEGER NN, NTOT, K, IER, XPANELS, YPANELS
      DOUBLE PRECISION CHISQPART
      CHARACTER*30 DATAFNAME, OUTFILE, OUTFORMAT
      CHARACTER TXT*45
      DOUBLE COMPLEX CBETA
      CHARACTER*9 FITNM
      DIMENSION FITNM(8)
      INCLUDE 'header.f'
      DOUBLE PRECISION A(NPARMAX),GIN(NPARMAX)
      DATA FITNM /'NS', 'AL0S', 'ALPS', 'M02S',
     & 'NG', 'AL0G', 'ALPG', 'M02G'/

#ifndef NOPGPLOT
      INTEGER PGBEG
      LOGICAL DRAW
      REAL NORMG, NORMSEA
      REAL XB, LGXB, LASTLGXB, YG, YS, LASTYG, LASTYS
#endif


*     Clear output buffer

      CALL FLUSHOUT()
      
*     Transfer parameters to common block

      DO 10 K = 1, NPARMAX
        PAR(K) = A(K)
 10   CONTINUE

*     Initialize GPDs
      CALL INITGPD

      NTOT = 0
      CHISQ = 0.0d0

      OPEN (UNIT = 11, FILE = 'FIT.INI', STATUS = 'OLD')
      READ (11, *) OUTFILE
      READ (11, *) OUTFORMAT, XPANELS, YPANELS
#ifndef NOPGPLOT
*    Initailization of plotting on a grid of panels
      IF (IFLAG .EQ. 3)  THEN
        IF (OUTFORMAT(:2) .EQ. '/X') THEN
*         Ploting to /XSERVE
          IER = PGBEG(0, "/XSERVE", XPANELS, YPANELS)
        ELSE
*         Ploting to /PSDRIV
        IER = PGBEG(0, 
     &  OUTFILE(1:MAX(1,INDEX(OUTFILE//' ',' ')-1)) // '.ps' //
     &  OUTFORMAT(1:MAX(1,INDEX(OUTFORMAT//' ',' ')-1)),
     &  XPANELS, YPANELS)
        END IF
        IF (IER.NE.1) STOP
        CALL PGSCH(1.8)
      END IF
#endif

*     Initialization of printing of the results tables to .out

      IF (IFLAG .EQ. 3)  THEN
        OPEN (UNIT = 21, FILE = 
     &  OUTFILE(1:MAX(1,INDEX(OUTFILE//' ',' ')-1)) // '.out',
     &  STATUS = 'UNKNOWN')
        WRITE (21, 803) "Gepard parameters:"
        WRITE (21, 804) SPEED, P, NF
        WRITE (21, 806) SCHEME, ANSATZ
        IF ( ANSATZ(:3) .EQ. 'FIT' ) THEN
          WRITE (21, 803) "Final fit parameters:"
          WRITE (21, 807) 'Q02      ', PAR(1)
          DO 13 K = 1, 4
          WRITE (21, 807) FITNM(K), PAR(10+K)
 13       WRITE (21, 807) FITNM(4+K), PAR(20+K)
        END IF
        WRITE (21, *) 
*     Same goes to .min
        WRITE (6, 803) "Gepard parameters:"
        WRITE (6, 804) SPEED, P, NF
        WRITE (6, 806) SCHEME, ANSATZ
        WRITE (6, *) 
      END IF

*     Process only files specified in FIT.INI between
*      'START' and 'STOP'

 15   READ (11, *) DATAFNAME
      IF (DATAFNAME(:5) .NE. 'START') GOTO 15
 18   READ (11, *) DATAFNAME
      IF (DATAFNAME(:4) .EQ. 'STOP') THEN
        GOTO 20
      ELSE
        CALL PROCDATA (DATAFNAME, IFLAG, NN, CHISQPART)
        CHISQ = CHISQ + CHISQPART
        NTOT = NTOT + NN
        IF (IFLAG .EQ. 3) WRITE (21, *)
      END IF
      GOTO 18
 20   CLOSE(11)

      IF (IFLAG .EQ. 3)  THEN
         WRITE (21, *) 'FREE PARAMS = ', NPAR
         WRITE (21, *) 'D.O.F = ', NTOT - NPAR
         WRITE (21, *) 'CHISQ = ', CHISQ
         CLOSE (21)
#ifndef NOPGPLOT
*     Draw next-to-last plot containing just parameter values
        CALL PGSCH(1.3)
        CALL PGPAGE
        CALL PGVSTD
        CALL PGSWIN (0., 12., 0., 12.)
        CALL PGTEXT (-1., 12., "Gepard parameters:")
        WRITE(TXT, 804) SPEED, P, NF
        CALL PGTEXT (1., 11., TXT)
        WRITE(TXT, 806) SCHEME, ANSATZ
        CALL PGTEXT (1., 10., TXT)
        IF ( ANSATZ(:3) .EQ. 'FIT' ) THEN
          CALL PGTEXT (-1., 9., "Final fit parameters:")
          WRITE(TXT, 807) 'Q02      ', PAR(1)
          CALL PGTEXT (-1., 8., TXT)
          DO 30 K = 1, 4
          WRITE (TXT, 807) FITNM(K), PAR(10+K)
          CALL PGTEXT (-1., REAL(8-K), TXT)
          WRITE (TXT, 807) FITNM(4+K), PAR(20+K)
 30       CALL PGTEXT (-1., REAL(4-K), TXT)
          WRITE (TXT, 808) NTOT, NPAR, NTOT-NPAR, CHISQ
          CALL PGTEXT (-1., -1., TXT)
        END IF
*     Draw last plot containing PDFs
        CALL PGSCH(1.8)
        CALL PGPAGE
        CALL PGVSTD
        DRAW = .FALSE.
        IF (PAR(1) .LT. 3.0d0) THEN
            CALL PGSWIN (-4., -1., -2., 5.)
        ELSEIF (PAR(1) .LT. 8.0d0) THEN
            CALL PGSWIN (-4., -1., 0., 12.)
        ELSEIF (PAR(1) .LT. 22.0d0) THEN
            CALL PGSWIN (-4., -1., 0., 24.)
        ELSE
            CALL PGSWIN (-4., -1., 0., 50.)
        ENDIF
        CALL PGBOX ('BCLNST', 0.0, 0, 'BCNST1', 0.0, 0)
        CALL PGLAB('x\\dBJ\\u', 'x\\gS(x)(blue), xg(x)(red)', 
     &       'Resulting PDFs')
        NORMSEA = PAR(11) / 
     &              CBETA(COMPLEX(2.0D0-PAR(12),0.0d0), (8.0d0,0.0d0))
        NORMG = PAR(21) / 
     &              CBETA(COMPLEX(2.0D0-PAR(22),0.0d0), (6.0d0,0.0d0))
        DO 40 LGXB = -4., -1., 0.02
          XB = 10**(LGXB)
          YS = NORMSEA * XB**(1.0d0-PAR(12)) * (1.0d0-XB)**7.0
          YG = NORMG * XB**(1.0d0-PAR(22)) * (1.0d0-XB)**5.0
          IF (DRAW) THEN
            CALL PGSCI(2)
            CALL PGSLS(2)
            CALL PGMOVE(LASTLGXB, LASTYG)
            CALL PGDRAW(LGXB, YG)
            CALL PGSCI(4)
            CALL PGSLS(1)
            CALL PGMOVE(LASTLGXB, LASTYS)
            CALL PGDRAW(LGXB, YS)
          END IF
          LASTLGXB = LGXB
          LASTYG = YG
          LASTYS = YS
          DRAW = .TRUE.
 40     CONTINUE
        CALL PGEND
#endif
      END IF

803   FORMAT (1X, /, 1X, '-----    ', A, '    -----', /)
804   FORMAT (1X, 'SPEED = ', I1, 5X, 'P = ', I1, 5X, 'NF = ', I1)
806   FORMAT (1X, 'SCHEME = ', A, 3X, 'ANSATZ = ', A)
807   FORMAT (1X,  A15, ' = ', 3X, G9.3)
808   FORMAT (1X, 'd.o.f = ', I3, ' - ', I1, ' = ', I3, '   ',
     &            '\\gx\\u2\\d = ', F8.1)

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
C     SUBROUTINE PROCDATA (FNAME, IFLAG, N, CHISQPART)
C
C     CHARACTER  FNAME*15
C     INTEGER IFLAG, N
C     DOUBLE PRECISION CHISQPART
C  INPUTS
C      FNAME  --  Name of the file containing data
C      IFLAG  --  If .EQ. 3 then print results out, see minuit doc
C  OUTPUT
C   CHISQPART --  contribution to chi-square from processed data set
C          N  --  Number of points in a given data set
C  IDENTIFIERS
C          X  --  Array with x-values of experimental points
C          Y  --  Array with y-values of experimental points
C         DY  --  Array with errors of y-values
C         W2  -- (gamma*-proton) invariant mass squared
C         XI  -- DVCS scaling parameter or Bjorken x
C       DEL2  -- DVCS asymmetry parameter (P2-P1)^2
C         Q2  -- photon virtuality squared
C  PARENTS
C      FCN
C  CHILDREN
C      PARSIGMA, SIGMA, F2F, PGPLOT routines
C  SOURCE
C
      SUBROUTINE PROCDATA (FNAME, IFLAG, NN, CHISQPART)

      IMPLICIT NONE
      CHARACTER  FNAME*30
      INTEGER IFLAG, NN
      DOUBLE PRECISION CHISQPART
      CHARACTER YOBS*8
      INTEGER K
      DOUBLE PRECISION X, YY, DY, THY, DIFSG
      DOUBLE PRECISION STAT, SYS
      DOUBLE PRECISION W
      DOUBLE PRECISION WIN, Q2IN
      DOUBLE PRECISION XBJ, XBJIN
      DOUBLE PRECISION PARSIGMA, SIGMA
#ifndef NOPGPLOT
      LOGICAL DRAW
      REAL LASTX, LASTY
#endif
      INCLUDE 'header.f'

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

        PROCESS = 'DVCS'

        READ (12,*) W
        READ (12,*) Q2
        READ (12,*) NN

        XI = Q2 / ( 2.0d0 * W**2 + Q2)
#ifndef NOPGPLOT
        IF (IFLAG .EQ. 3) THEN
          CALL PGSWIN (0., 1., -1.3, 2.)
          CALL PGBOX ('BCNST1', 0.0, 0, 'BCLNST', 0.0, 0)
          CALL PGLAB("-t", 'd\\gs/dt', FNAME)
        END IF
#endif

        DO 110 K = 1, NN
        READ (12, *) X, YY, STAT, SYS
        DY = SQRT( STAT*STAT + SYS*SYS )
        MTIND = NMTS + K
        THY = PARSIGMA()
        CHISQPART = CHISQPART + ( (THY - YY)**2 / DY**2 )
        IF (IFLAG .EQ. 3)  THEN
          DIFSG = (THY - YY) / DY
          WRITE (21, 901) X, THY, YY, DY, DIFSG
#ifndef NOPGPLOT
          CALL PGERRY(1, SNGL(X), LOG10(SNGL(YY-DY)), 
     &                  LOG10(SNGL(YY+DY)), 3.0)
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

        PROCESS = 'DVCS'

        READ (12,*) WIN
        READ (12,*) Q2IN
        READ (12,*) NN

        DO 120 K = 1, NN
        READ (12, *) X, YY, STAT, SYS
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
        CHISQPART = CHISQPART + ( (THY - YY)**2 / DY**2 )
        IF (IFLAG .EQ. 3)  THEN
          DIFSG = (THY - YY) / DY
          WRITE (21, 901) X, THY, YY, DY, DIFSG
#ifndef NOPGPLOT
          IF (WIN .LT. 0.) THEN
            CALL PGERRY(1, SNGL(X), SNGL(YY-DY), SNGL(YY+DY), 3.0)
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
            CALL PGERRY(1,SNGL(X), LOG10(SNGL(YY-DY)), 
     &              LOG10(SNGL(YY+DY)), 3.0)
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

        PROCESS = 'DIS'

        READ (12,*) XBJIN
        READ (12,*) Q2IN
        READ (12,*) NN

        DO 130 K = 1, NN
        READ (12, *) X, YY, STAT, SYS
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
        THY = F2(P)
        CHISQPART = CHISQPART + ( (THY - YY)**2 / DY**2 )
        IF (IFLAG .EQ. 3)  THEN
          DIFSG = (THY - YY) / DY
          WRITE (21, 901) X, THY, YY, DY, DIFSG
#ifndef NOPGPLOT
          CALL PGERR1(6, SNGL(X), SNGL(YY), SNGL(DY), 3.0)
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
