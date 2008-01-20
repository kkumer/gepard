C     ****h* gepard/fit.f
C  FILE DESCRIPTION
C    Fitting of GPD model to DVCS and DIS  experimental data.
C
C    $Id: fit.F 76 2007-04-23 15:39:11Z kkumer $
C     *******

C     ****p* fit.f/FIT
C  NAME
C    FIT   --  Determines parameters of GPDs by fitting to DVCS and DIS
C              experimental data.
C            
C  DESCRIPTION
C             Calls minuit subroutine for minimization and
C             prints results in various ways.
C    Output goes to files with extensions .out (tabular representation
C          of fit), .ps (graphical representation of fit), and .mnt
C          (Minuit output)
C  SYNOPSIS

      PROGRAM FIT

C
C
C  CHILDREN
C      READPAR, INIT, INITC, MINUIT, FCN (via MINUIT)
C
C  SOURCE
C

      IMPLICIT NONE
      INTEGER ISINTER, NARGS
      CHARACTER FNAME*20
      EXTERNAL FCN
      INCLUDE 'header.f'

      COMMON /INTERACTIVE/ ISINTER

      CALL READPAR

*   Overriding filenames from READPAR if they are specified
*   on commad line:

      NARGS = IARGC()
*   *.dat (1st argument) is file with datasets specification
      IF (NARGS .GE. 1) THEN
        CALL GETARG(1, FNAME)
        DATFILE = FNAME(1:MAX(1,INDEX(FNAME//' ',' ')-1)) // '.dat'
      END IF
*   *.out *.ps *.mnt (2nd argument) is where results go
      IF (NARGS .GE. 2) THEN
        CALL GETARG(2, FNAME)
        OUTFILE = FNAME(1:MAX(1,INDEX(FNAME//' ',' ')-1)) // '.out'
      ELSE 
*     default is to be same as DATFILE but .dat -> .out
        OUTFILE = DATFILE(1:INDEX(DATFILE,'.dat')-1) // '.out'
      END IF
*   *.cmd (3rd argument) is from where MINUIT commands are read
      IF (NARGS .GE. 3) THEN
        CALL GETARG(3, FNAME)
        CMDFILE = FNAME(1:MAX(1,INDEX(FNAME//' ',' ')-1)) // '.cmd'
      END IF
*    other arguments are ignored


*   Reading from first line of CMDFILE whether we want interactive Minuit session.
      OPEN (UNIT = 17, FILE = CMDFILE, STATUS = 'OLD')
      READ (17, *) ISINTER
      CLOSE (17)

      FFTYPE = 'SINGLET'
      PROCESS = 'DVCS'
      CALL INIT
      PROCESS = 'DIS'
      CALL INIT

      CALL INITC

      IF (ISINTER .NE. 1) THEN
*   We work in batch mode
*      File where Minuit batch commands are
        OPEN (UNIT = 5, FILE = CMDFILE, STATUS = 'OLD')
*      File wher MINUIT (and STDOUT) output is written
        OPEN (UNIT = 6, FILE = 'fit.mnt', STATUS = 'UNKNOWN')
      ELSE
*   We work in interactive mode
        WRITE (*, 801) CMDFILE
      END IF

      CALL MINUIT(FCN,0)  

801   FORMAT (1X, // 1X, 72('-') / 1X, " -- Interactive mode ! --" / 1X,
     & "You can answer 'SET INPUT 18 ", A, "' to the querry below" / 1X,
     & "(Second question about rewinding the unit is irrelevant.)" / 1X,
     & 72('-') //)
      STOP
      END
C     ****
