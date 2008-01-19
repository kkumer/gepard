C     ****h* gepard/mninterface.f
C  FILE DESCRIPTION
C    Mathematica - Minuit interface (Fortran)
C
C    $Id:$
C     *******

C     ****s* mninterface.f/FITINIT
C  NAME
C     FITINIT  --   Initialization of fitting procedure
C  DESCRIPTION
C     Initializes GeParD (MB contour, evolved Wilson coefs) and MINUIT
C  SYNOPSIS

      SUBROUTINE FITINIT

C  PARENTS
C    MinuitInit
C  CHILDREN
C    READPAR, INIT, INITC, MNINIT, MNSETI
C  SOURCE
C
      IMPLICIT NONE
      INCLUDE 'header.f'


      FFTYPE = 'SINGLET'
      PROCESS = 'DVCS'
      CALL INIT
      PROCESS = 'DIS'
      CALL INIT

      CALL INITC

*      File where MINUIT (and STDOUT) output is written
        OPEN (UNIT = 6, FILE = 
     &    OUTFILE(1:INDEX(OUTFILE,'.out')-1) // '.mnt',
     &    STATUS = 'UNKNOWN')


      CALL MNINIT(5, 6, 7)

      CALL MNSETI('Fitting to DVCS and DIS data.')  

      RETURN
      END
C     ****



C     ****s* mninterface.f/MPAR
C  NAME
C     MPAR  --   Define a Minuit parameter and assigns it value
C  SYNOPSIS

      SUBROUTINE MPAR(ID, SIZ, PNAM, VSTRT, STP, LO, HI)

      IMPLICIT NONE
      INTEGER ID, SIZ
      CHARACTER PNAM*10
      DOUBLE PRECISION VSTRT, STP, LO, HI

C  INPUTS 
C                  ID  --  parameter number
C                 SIZ  --  length of PNAM
C                PNAM  --  parameter name
C               VSTRT  --  starting value of parameter
C                 STP  --  starting step size or approx. parameter error
C                  LO  --  lower bound on parameter
C                  HI  --  upper bound on parameter
C                          (if LO=HI=0 parameter is considered unbounded)
C  PARENTS
C     MinuitSetParameter
C  CHILDREN
C     MNPARM
C  BUGS    
C     IERFLG should be returned to caller
C  SOURCE
C
      INTEGER IERFLG
      CHARACTER TPNAM*10

      TPNAM = PNAM(1:SIZ)

      CALL MNPARM(ID, TPNAM, VSTRT, STP, LO, HI, IERFLG)
        IF (IERFLG .NE. 0) THEN
          WRITE(6, '(A,I)') ' Unable to define parametr no. ', ID
          STOP
        END IF

      RETURN
      END
C     ****


C     ****s* mninterface.f/MCOM
C  NAME
C     MCOM  --  Executes Minuit command (as a character string)
C  SYNOPSIS

      SUBROUTINE MCOM(SIZ, CMD, IERFLG)

      IMPLICIT NONE
      INTEGER SIZ, IERFLG
      CHARACTER CMD*100 

C  INPUTS 
C                 SIZ  --  length of CMD
C                 CMD  --  Minuit command (as a character string)
C  OUTPUT 
C              IERFLG  --  error code, =0 if command executed normally
C  PARENTS
C     MinuitCommand
C  CHILDREN
C     MNCOMD
C  SOURCE
C
      CHARACTER TCMD*100
      EXTERNAL FCN
      INCLUDE 'header.f'

*   First flush the .mnt file
C      CLOSE(6)
C      OPEN(UNIT = 6, FILE = 
C     &  OUTFILE(1:INDEX(OUTFILE,'.out')-1) // '.mnt', STATUS = "OLD")
C 10   READ(6, *, END=20)
C      GOTO 10
C 20   BACKSPACE (11)

      TCMD = CMD(1:SIZ)
      CALL MNCOMD(FCN, TCMD, IERFLG, 0)

      RETURN
      END
C     ****


C     ****s* mninterface.f/GETPAR
C  NAME
C     GETPAR  --  Gets the current value of the parameter
C  SYNOPSIS

      SUBROUTINE GETPAR(ID, VAL, ERROR, IVARBL)

      IMPLICIT NONE
      INTEGER ID, IVARBL
      DOUBLE PRECISION VAL, ERROR

C  INPUTS 
C                  ID  --  parameter number
C  OUTPUT 
C                 VAL  --  current parameter value
C               ERROR  --  current estimate of parameter uncertainty
C              IVARBL  --  internal parmeter number (0 if constant,
C                          negative if undefined)
C  PARENTS
C     MinuitGetParameter
C  CHILDREN
C     MNPOUT
C  SOURCE
C
      DOUBLE PRECISION BND1, BND2
      CHARACTER CHNAM*10


      CALL MNPOUT(ID, CHNAM, VAL, ERROR, BND1, BND2, IVARBL)

      RETURN
      END
C     ****

C     ****s* mninterface.f/PARINIT
C  NAME
C     PARINIT  --   sets optional parameters
C  SYNOPSIS

      SUBROUTINE PARINIT(INSPEED, INP, INSCHEMESIZ, INSCHEME, 
     &   INANSATZSIZ, INANSATZ, 
     &   INDATFILESIZ, INDATFILE, 
     &   INOUTFILESIZ, INOUTFILE
     &   )

      IMPLICIT NONE
      INTEGER INSPEED, INP
      INTEGER INSCHEMESIZ, INANSATZSIZ 
      INTEGER INDATFILESIZ, INOUTFILESIZ
      CHARACTER INSCHEME*6, INANSATZ*7 
      CHARACTER INDATFILE*16, INOUTFILE*16
      INCLUDE 'header.f'

C  INPUTS 
C           IN*SIZ  --  length of IN*
C              IN*  --  input value of *
C  PARENTS
C     GepardInitInternal
C  CHILDREN
C     READPAR
C  SOURCE
C

*  read the defaults:

      CALL READPAR

*  override what's requested:

      IF ( INSPEED .GT. 0 ) SPEED = INSPEED

      IF ( INP .GE. 0 ) P = INP

      IF ( INSCHEME(1:INSCHEMESIZ) .NE. 'DFLT' )
     &         SCHEME = INSCHEME(1:INSCHEMESIZ)

      IF ( INANSATZ(1:INANSATZSIZ) .NE. 'DFLT' )
     &         ANSATZ = INANSATZ(1:INANSATZSIZ)
      
      IF ( INDATFILE(1:INDATFILESIZ) .NE. 'DFLT' )
     &        DATFILE = INDATFILE(1:INDATFILESIZ)//'.dat'

      IF ( INOUTFILE(1:INOUTFILESIZ) .NE. 'DFLT' ) THEN
              OUTFILE = INOUTFILE(1:INOUTFILESIZ)//'.out'
      ELSE
*       same rootname as .dat
        OUTFILE = DATFILE(1:INDEX(DATFILE,'.dat')-1) // '.out'
      END IF

      RETURN
      END
C     ****
