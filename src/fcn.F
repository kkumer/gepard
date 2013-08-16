C     ****h* gepard/fcn.F
C  FILE DESCRIPTION
C    subroutine returning the chi-sqare
C
C    $Id$
C     *******

C     ****s* fcn.F/FCN
C  NAME
C     FCN  --  called by minuit, returns chi-square
C            
C  DESCRIPTION
C    Reads experimental datasets and calculates chi-square by
C    comparing to theoretical prediction functions.
C    Writes out theoretical line points for fit plotting.
C  SYNOPSIS
C
      SUBROUTINE FCN(NPAR,GIN,CHI2,A,IFLAG,FUTIL)

      IMPLICIT NONE
      INTEGER NPAR, IFLAG
      DOUBLE PRECISION CHI2, FUTIL
#include "header.F"
      DOUBLE PRECISION A(NPARMAX),GIN(NPARMAX)
C    
C  INPUTS 
C                NPAR  --  Number of variable parameters
C                   A  --  Values of all parameters
C               IFLAG  --  If .EQ. 3 then print results out, see minuit doc
C  OUTPUT
C               CHI2   --  Value of total chi-square = CHISQ(0)
C  IDENTIFIERS
C     FITPAR  --  Equal to A
C  CHISQPART  --  chi-square produced by a single data set
C
C  CHILDREN
C      GETMBGPD, GETMBGPDMMA, PROCDATA, PGPLOT routines, FLUSHOUT
C
C  PARENTS
C      MINUIT, FIT (via MINUIT)
C
C  SOURCE
C
      INTEGER NN, K, XPANELS, YPANELS, ICOL
      INTEGER PRPARMAX, NSUBTOT, IPANEL
      DOUBLE PRECISION CHISUBTOT, CHISQPART
      CHARACTER*30 DATAFNAME, OUTFORMAT
      CHARACTER*20 FIGLABEL, XLABEL, YLABEL, XTYPE, YTYPE
      CHARACTER TXT*50, CHTXT*15, CHTXTOT*60
      REAL XMIN, XMAX, YMIN, YMAX
      INTEGER IVAR
      CHARACTER*10 NAML, NAMR
      CHARACTER*1 VFLAGL, VFLAGR
      DOUBLE PRECISION VALL, VALR, ERL, ERR, BND1, BND2
      DOUBLE PRECISION BCA0
*      Parameters to be printed on the plot (only second digit)
      PARAMETER (PRPARMAX = 5)
      INTEGER PRPAR(PRPARMAX)
      DATA PRPAR / 1, 2, 4, 8, 9 /

#ifdef PGPLOT
      INTEGER PGBEG, IER, COL(6)
      LOGICAL DRAW
      DOUBLE PRECISION HX(2), BH(2)
      REAL XBJ, LGXBJ, DELLGXBJ, LASTLGXBJ, YG, YS, LASTYG, LASTYS
      REAL CHSIZE
      PARAMETER (DELLGXBJ = 0.02)
*      Colors used by PGPLOT:
      DATA COL /1, 2, 4, 8, 12, 13/
      DRAW = .FALSE.
#endif

      CHTXTOT = ''

*     Clear output buffer

      CALL FLUSHOUT()
      
*     Transfer parameters to common block

      DO 10 K = 1, NPARMAX
        PAR(K) = A(K)
 10   CONTINUE

*     Initialize (only forward!) GPDs
      DEL2 = 0.0d0
      PROCESS = "DIS"
#ifdef MMA
      IF (ANSATZ(:3) .EQ. 'MMA') THEN
        CALL GETMBGPDMMA
      ELSE
        CALL GETMBGPD
      END IF
#else
      IF (ANSATZ(:3) .EQ. 'MMA') THEN
        WRITE(*,*) 'ANSATZ=MMA unavailable. Recompile with -D MMA'
      ELSE
        CALL GETMBGPD
      END IF
#endif

      ICOL = 0
      IPANEL = 0

      NDATAPTS(0) = 0
      CHISQ(0) = 0.0d0
      CHIN = 0

      OPEN (UNIT = 11, FILE = DATFILE, STATUS = 'OLD')
      READ (11, *) OUTFORMAT, XPANELS, YPANELS
#ifdef PGPLOT
*    Initailization of plotting on a grid of panels
      IF (IFLAG .EQ. 3)  THEN
        IF (OUTFORMAT(:2) .EQ. '/X') THEN
*         Ploting to /XSERVE
          IER = PGBEG(0, "/XSERVE", XPANELS, YPANELS)
        ELSE
*         Ploting to /PSDRIV
        IER = PGBEG(0, 
     &  OUTFILE(1:INDEX(OUTFILE,'.out')-1) // '.ps' //
     &  OUTFORMAT(1:MAX(1,INDEX(OUTFORMAT//' ',' ')-1)),
     &  XPANELS, YPANELS)
        END IF
        IF (IER.NE.1) STOP
        CALL PGSCH(1.8)
      END IF
#endif

*     Initialization of printing of the results tables to .out

      IF (IFLAG .EQ. 3)  THEN
        OPEN (UNIT = 21, FILE = OUTFILE, STATUS = 'UNKNOWN')
        WRITE (21, 803) "Gepard parameters:"
        WRITE (21, 804) SPEED, P, NF
        WRITE (21, 806) SCHEME, Q02
        WRITE (21, *) 
*     Same goes to .min
        WRITE (6, 803) "Gepard parameters:"
        WRITE (6, 804) SPEED, P, NF
        WRITE (6, 806) SCHEME, Q02
        WRITE (6, *) 
      END IF

*     Process only files specified in DATFILE between
*      'START' and 'STOP'

 15   READ (11, *) DATAFNAME
      IF (DATAFNAME(:5) .NE. 'START') GOTO 15
 18   READ (11, *) DATAFNAME
      IF (DATAFNAME(:4) .EQ. 'STOP') THEN
        GOTO 20 
      ELSE
        IF (DATAFNAME(:6) .EQ. 'NEWFIG') THEN
*         New panel is desired.
*         Read next line with PGPLOT parameters.
          READ (11, *) FIGLABEL, XLABEL, YLABEL, 
     &      XMIN, XMAX, YMIN, YMAX, XTYPE, YTYPE
*         start plotting panel
#ifdef PGPLOT
          IF (IFLAG .EQ. 3) THEN
*         check if second X window is needed
            IPANEL = IPANEL + 1
            IF ( (IPANEL .EQ. XPANELS*YPANELS+1) .AND. 
     &           (OUTFORMAT(:2) .EQ. '/X') ) THEN
              CALL PGEND
              IER = PGBEG(0, "2/XSERVE", XPANELS, YPANELS)
              IF (IER.NE.1) STOP
              CALL PGSCH(1.8)
            END IF
            CALL PGPAGE
            CALL PGVSTD
            CALL PGSLW(2)
            CALL PGSWIN (XMIN, XMAX, YMIN, YMAX)
            CALL PGBOX (XTYPE, 0.0, 0, YTYPE, 0.0, 0)
            CALL PGLAB(XLABEL, YLABEL, FIGLABEL)
          END IF
#endif
        ELSE IF (DATAFNAME(:5) .EQ. 'CHISQ') THEN
          IF (IFLAG .EQ. 3)  THEN
*         Calculate partial chi-square and corresponding number
*         of data points used for fitting and write to the
*         next free space in arrays CHISQ and NDATAPTS ...
            NSUBTOT = 0
            CHISUBTOT = 0.0d0
            DO 19 K = 1, CHIN
              NSUBTOT = NSUBTOT + NDATAPTS(K)
              CHISUBTOT = CHISUBTOT + CHISQ(K)
 19         CONTINUE
            CHIN = CHIN + 1
            NDATAPTS(CHIN) = NDATAPTS(0) - NSUBTOT
            CHISQ(CHIN) = CHISQ(0) - CHISUBTOT
*         ... and also add to string CHXTOT for final printout on graph
            WRITE(CHTXT, 811) CHISQ(CHIN), NDATAPTS(CHIN)
            CHTXTOT = 
     &       CHTXTOT(1:MAX(1,INDEX(CHTXTOT//'   ','   ')-1)) // CHTXT
*         ... and to .out file
            WRITE(21, *) '====== SUBTOTAL ====='
            WRITE(21, 810) CHISQ(CHIN), NDATAPTS(CHIN)
*         ... and to title of last panel
            WRITE(TXT, 810) CHISQ(CHIN), NDATAPTS(CHIN)
#ifdef PGPLOT
            CALL PGQCH(CHSIZE)
            CALL PGSCH(1.4)
            CALL PGTEXT (XMIN, YMAX+(YMAX-YMIN)/30., TXT)
            CALL PGSCH(CHSIZE)
#endif
          END IF
        ELSE
#ifdef PGPLOT
          IF (IFLAG .EQ. 3)  THEN
            ICOL = ICOL + 1
            IF (ICOL .EQ. 5) ICOL = 1
            CALL PGSCI(COL(ICOL))
            CALL PGSLS(ICOL)
          END IF
#endif
          CALL PROCDATA (DATAFNAME, IFLAG, NN, CHISQPART)
          CHISQ(0) = CHISQ(0) + CHISQPART
          NDATAPTS(0) = NDATAPTS(0) + NN
#ifdef PGPLOT
          IF (IFLAG .EQ. 3)  THEN
            CALL PGSCI(1)
            CALL PGSLS(1)
          END IF
#endif
        END IF
        IF (IFLAG .EQ. 3) WRITE (21, *)
      END IF
      GOTO 18
 20   CLOSE(11)

      IF (IFLAG .EQ. 3)  THEN
*     Write parameter values, chi-squares and BCA to *.out
        IF ((ANSATZ(:3) .EQ. 'FIT') .OR. (ANSATZ(:3) .EQ. 'MMA')
     &   .OR. (ANSATZ(:1) .EQ. 'E')  ) THEN
          WRITE (21, 803) "Final fit parameters (*=fixed):"
          DO 25 K = 1, PRPARMAX
          CALL MNPOUT(10+PRPAR(K), NAML, VALL, ERL, BND1, BND2, IVAR)
          IF (IVAR .EQ. 0) THEN
            VFLAGL="*"
          ELSE
            VFLAGL=" "
          ENDIF
          NAML = VFLAGL//NAML
          CALL MNPOUT(20+PRPAR(K), NAMR, VALR, ERR, BND1, BND2, IVAR)
          IF (IVAR .EQ. 0) THEN
            VFLAGR="*"
          ELSE
            VFLAGR=" "
          ENDIF
          NAMR = VFLAGR//NAMR
          WRITE (21, 809) NAML, VALL, NAMR, VALR
 25       CONTINUE
        END IF
        WRITE (21, 803) '[partial chi-squares / datapoints]:'
        WRITE (21, *) CHTXTOT
        WRITE (21, 803) "Chi-square:"
        WRITE (21, *) 'FREE PARAMS = ', NPAR
        WRITE (21, *) 'D.O.F = ', NDATAPTS(0) - NPAR
        WRITE (21, *) 'CHISQ = ', CHISQ(0)
        WRITE (21, 803) 'Beam charge asymmetry'
*    FIXME: hardwired values for H1 here:
        CALL BCA(1, 82.0d0, 8.0d0, 0.0d0, BCA0)
        WRITE (21, *) 'BCA(PHI=0) = ', BCA0
        CLOSE (21)
#ifdef PGPLOT
*     Write parameter values, chi-squares and BCA to plot
        CALL PGSCH(1.3)
        CALL PGPAGE
        CALL PGVSTD
        CALL PGSWIN (0., 12., 0., 12.)
        CALL PGTEXT (-1., 13., "Gepard parameters:")
        WRITE(TXT, 804) SPEED, P, NF
        CALL PGTEXT (1., 12., TXT)
        WRITE(TXT, 806) SCHEME, Q02
        CALL PGTEXT (1., 11., TXT)
        IF ((ANSATZ(:3) .EQ. 'FIT') .OR. (ANSATZ(:3) .EQ. 'MMA')
     &      .OR. (ANSATZ(:1) .EQ. 'E')  ) THEN
          CALL PGTEXT (-1., 10., "Final fit parameters (*=fixed):")
          DO 30 K = 1, PRPARMAX
          CALL MNPOUT(10+PRPAR(K), NAML, VALL, ERL, BND1, BND2, IVAR)
          IF (IVAR .EQ. 0) THEN
            VFLAGL="*"
          ELSE
            VFLAGL=" "
          ENDIF
          NAML = VFLAGL//NAML
          CALL MNPOUT(20+PRPAR(K), NAMR, VALR, ERR, BND1, BND2, IVAR)
          IF (IVAR .EQ. 0) THEN
            VFLAGR="*"
          ELSE
            VFLAGR=" "
          ENDIF
          NAMR = VFLAGR//NAMR
          WRITE (TXT, 809) NAML, VALL, NAMR, VALR
          CALL PGTEXT (-1., REAL(10-K), TXT)
 30       CONTINUE
          CALL PGTEXT(-1.,4.,'[\\gx\\u2\\d\\dpartial\\u / datapoints]:')
          CALL PGTEXT (-1., 3., CHTXTOT)
          WRITE (TXT, 818) NDATAPTS(0), NPAR, NDATAPTS(0)-NPAR, CHISQ(0)
          CALL PGTEXT(-1.,2.,'Total \\gx\\u2\\d:')
          CALL PGTEXT (-1., 1., TXT)
          CALL PGTEXT(-1.,0.,'Beam charge asymmetry:')
          WRITE (TXT, 812) BCA0
          CALL PGTEXT (-1., -1., TXT)
        END IF
*     Draw plot of PDFs
        CALL PGSCH(1.8)
        CALL PGSCI(1)
        CALL PGPAGE
        CALL PGVSTD
        DRAW = .FALSE.
        IF (Q02 .LT. 3.0d0) THEN
            CALL PGSWIN (-4., -1., -2., 5.)
        ELSEIF (Q02 .LT. 8.0d0) THEN
            CALL PGSWIN (-4., -1., 0., 12.)
        ELSEIF (Q02 .LT. 22.0d0) THEN
            CALL PGSWIN (-4., -1., 0., 24.)
        ELSE
            CALL PGSWIN (-4., -1., 0., 50.)
        ENDIF
        CALL PGBOX ('BCLNST', 0.0, 0, 'BCNST1', 0.0, 0)
        CALL PGLAB('x', 'x\\gS(x)(blue), xg(x)(red)', 
     &       'Resulting PDFs')
        LGXBJ = -4.
        DO WHILE (LGXBJ .LE. -1.)
          XBJ = 10**(LGXBJ)
          PROCESS = 'DIS'
          CALL XSPACE(HX, DBLE(XBJ), 0.0d0)
          YS = HX(1)
          YG = HX(2)
          IF (DRAW) THEN
            CALL PGSCI(2)
            CALL PGSLS(2)
            CALL PGMOVE(LASTLGXBJ, LASTYG)
            CALL PGDRAW(LGXBJ, YG)
            CALL PGSCI(4)
            CALL PGSLS(1)
            CALL PGMOVE(LASTLGXBJ, LASTYS)
            CALL PGDRAW(LGXBJ, YS)
          END IF
          LASTLGXBJ = LGXBJ
          LASTYG = YG
          LASTYS = YS
          DRAW = .TRUE.
          LGXBJ = LGXBJ + DELLGXBJ
        END DO
*     Draw plot of GPD slopes B(x)
        CALL PGSCH(1.8)
        CALL PGSCI(1)
        CALL PGPAGE
        CALL PGVSTD
        DRAW = .FALSE.
        CALL PGSWIN (-4., -1., 0., 6.)
        CALL PGBOX ('BCLNST', 0.0, 0, 'BCNST1', 0.0, 0)
        CALL PGLAB('x', 'slope B(x) (blue=q, red=g)', 
     &       'GPD slopes')
        LGXBJ = -4.
        DO WHILE (LGXBJ .LE. -1.)
          XBJ = 10**(LGXBJ)
          CALL SLOPEH(BH, DBLE(XBJ))
          YS = BH(1)
          YG = BH(2)
          IF (DRAW) THEN
            CALL PGSCI(2)
            CALL PGSLS(2)
            CALL PGMOVE(LASTLGXBJ, LASTYG)
            CALL PGDRAW(LGXBJ, YG)
            CALL PGSCI(4)
            CALL PGSLS(1)
            CALL PGMOVE(LASTLGXBJ, LASTYS)
            CALL PGDRAW(LGXBJ, YS)
          END IF
          LASTLGXBJ = LGXBJ
          LASTYG = YG
          LASTYS = YS
          DRAW = .TRUE.
          LGXBJ = LGXBJ + DELLGXBJ
        END DO
        CALL PGIDEN
        CALL PGEND
#endif
*   set to actual number of d.o.f
        NDATAPTS(0) = NDATAPTS(0) - NPAR
      END IF
*   return total chi-square value to MINUIT:
      CHI2 = CHISQ(0)

803   FORMAT (1X, /, 1X, '-----    ', A, '    -----', /)
804   FORMAT (1X, 'SPEED = ', I1, 5X, 'P = ', I1, 5X, 'NF = ', I1)
806   FORMAT (1X, 'SCHEME = ', A, 3X, 'Q02 = ', F5.1)
809   FORMAT (3X,  A6, '=', 1X, G8.2, 5X, A6, '=', 1X, G8.2)
810   FORMAT (1X, 'CHISQPART = ', F5.1, 3X, 'PTS = ', I3)
811   FORMAT (1X, '[', F5.1, '/', I3, ']', 3X)
#ifdef PGPLOT
812   FORMAT (1X, 'BCA(\\gf=0) = ', F7.3)
818   FORMAT (1X, 'd.o.f = ', I3, ' - ', I2, ' = ', I3, '   ',
     &            'chi-square = ', F8.1)
#endif

      RETURN
      END
C     ****


C     ****s* fit.F/PROCDATA
C  NAME
C    PROCDATA   --  Process a single data set
C  DESCRIPTION
C    Reads experimental datasets and calculates chi-square by
C    comparing to theoretical prediction functions.
C    Writes out theoretical line points for fit plotting.
C  SYNOPSIS
C
      SUBROUTINE PROCDATA (FNAME, IFLAG, NN, CHISQPART)

      IMPLICIT NONE
      CHARACTER  FNAME*30
      INTEGER IFLAG, NN
      DOUBLE PRECISION CHISQPART
C
C  INPUTS
C      FNAME  --  Name of the file containing data
C      IFLAG  --  If .EQ. 3 then print results out, see minuit doc
C  OUTPUT
C   CHISQPART --  contribution to chi-square from processed data set
C         NN  --  Number of points in a given data set
C  IDENTIFIERS
C          X  --  Array with x-values of experimental points
C          Y  --  Array with y-values of experimental points
C         DY  --  Array with errors of y-values
C         W2  -- (gamma*-proton) invariant mass squared
C         XI  -- DVCS scaling parameter or Bjorken x
C       DEL2  -- DVCS asymmetry parameter (P2-P1)^2
C         Q2  -- photon virtuality squared
C         LS  -- line style
C  PARENTS
C      FCN
C  CHILDREN
C      PARSIGMA, SIGMA, F2F, PGPLOT routines
C  SOURCE
C
      CHARACTER YOBS*8
      INTEGER K
      DOUBLE PRECISION X, YY, DY, THY, DIFSG
      DOUBLE PRECISION STAT, SYS
      DOUBLE PRECISION W
      DOUBLE PRECISION WIN, Q2IN
      DOUBLE PRECISION XBJ, XBJIN
      DOUBLE PRECISION PARSIGMA, SIGMA
#include "header.F"
#ifdef PGPLOT
      INTEGER LS
      LOGICAL DRAW
      REAL LASTX, LASTY
      DRAW = .FALSE.
#endif

      CHISQPART = 0.d0

      IF (IFLAG .EQ. 3) THEN
        WRITE (21, 902) FNAME
#ifdef PGPLOT
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

        DO 110 K = 1, NN
        READ (12, *) X, YY, STAT, SYS
        DY = SQRT( STAT*STAT + SYS*SYS )
        MTIND = NMTS + K
        THY = PARSIGMA()
        CHISQPART = CHISQPART + ( (THY - YY)**2 / DY**2 )
        IF (IFLAG .EQ. 3)  THEN
          DIFSG = (THY - YY) / DY
          WRITE (21, 901) X, THY, YY, DY, DIFSG
#ifdef PGPLOT
          CALL PGQLS(LS)
          CALL PGSLS(1)
          CALL PGERRY(1, SNGL(X), LOG10(SNGL(YY-DY)), 
     &                  LOG10(SNGL(YY+DY)), 3.0)
          CALL PGSLS(LS)
          IF (DRAW) THEN
            CALL PGMOVE(LASTX, LASTY)
            CALL PGDRAW(SNGL(X), LOG10(SNGL(THY)))
          END IF
          LASTX = SNGL(X)
          LASTY = LOG10(SNGL(THY)) 
          DRAW = .TRUE.
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
        ELSE IF (Q2IN .LT. 0) THEN
          Q2 = X
          W = WIN
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
#ifdef PGPLOT
          IF (WIN .LT. 0.) THEN
            CALL PGQLS(LS)
            CALL PGSLS(1)
            CALL PGERRY(1, SNGL(X), SNGL(YY-DY), SNGL(YY+DY), 3.0)
            CALL PGSLS(LS)
            IF (DRAW) THEN
              CALL PGMOVE(LASTX, LASTY)
              CALL PGDRAW(SNGL(X), SNGL(THY))
            END IF
            LASTX = SNGL(X)
            LASTY = SNGL(THY)
            DRAW = .TRUE.
          ELSE
            CALL PGQLS(LS)
            CALL PGSLS(1)
            CALL PGERRY(1,SNGL(X), LOG10(SNGL(YY-DY)), 
     &              LOG10(SNGL(YY+DY)), 3.0)
            CALL PGSLS(LS)
            IF (DRAW) THEN
              CALL PGMOVE(LASTX, LASTY)
              CALL PGDRAW(SNGL(X), LOG10(SNGL(THY)))
            END IF
            LASTX = SNGL(X)
            LASTY = LOG10(SNGL(THY)) 
            DRAW = .TRUE.
          END IF
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
        ELSE IF (Q2IN .LT. 0) THEN
          Q2 = X
          XBJ = XBJIN
        ELSE
          CALL ERROR ('GeParD', 'PROCDATA',
     &    'Either XBJ or Q2 in ' // FNAME // ' should be negative!',
     &    4, 2)
          STOP
        END IF
        XI = XBJ
        CALL F2F
        THY = F2(P)
        CHISQPART = CHISQPART + ( (THY - YY)**2 / DY**2 )
        IF (IFLAG .EQ. 3)  THEN
          DIFSG = (THY - YY) / DY
          WRITE (21, 901) X, THY, YY, DY, DIFSG
#ifdef PGPLOT
          CALL PGQLS(LS)
          CALL PGSLS(1)
          CALL PGERR1(6, SNGL(X), SNGL(YY), SNGL(DY), 3.0)
          CALL PGSLS(LS)
          IF (DRAW) THEN
            CALL PGMOVE(LASTX, LASTY)
            CALL PGDRAW(SNGL(X), SNGL(THY))
          END IF
          LASTX = SNGL(X)
          LASTY = SNGL(THY) 
          DRAW = .TRUE.
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
