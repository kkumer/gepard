C     ****h* gepard/fit.f
C  FILE DESCRIPTION
C    Fitting of DVCS (and DIS in future) experimental data.
C
C    $Id: fit.F 76 2007-04-23 15:39:11Z kkumer $
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
