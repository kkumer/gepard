* This is interface to Fortran Minuit routies

      SUBROUTINE FITINIT

      IMPLICIT NONE
      CHARACTER OUTFILE*20
      INCLUDE 'header.f'


      CALL READPAR
      FFTYPE = 'SINGLET'

      PROCESS = 'DVCS'
      CALL INIT
      PROCESS = 'DIS'
      CALL INIT

      CALL INITC

*     File where Minuit batch commands are
      OPEN (UNIT = 11, FILE = 'FIT.INI', STATUS = 'OLD')
      READ (11, *) OUTFILE
      CLOSE (11)
*     File for writing out Minuit output
      OPEN (UNIT = 6, FILE = 
     &    OUTFILE(1:MAX(1,INDEX(OUTFILE//' ',' ')-1)) // '.min',
     &    STATUS = 'UNKNOWN')


      CALL MNINIT(5, 6, 7)

      CALL MNSETI('Fitting to DVCS and DIS data.')  

      RETURN
      END
C     ****



      SUBROUTINE MPAR(ID, SIZ, PNAM, VSTRT, STP, LO, HI)

      IMPLICIT NONE
      INTEGER SIZ, IERFLG, ID
      DOUBLE PRECISION VSTRT, STP, LO, HI
      CHARACTER PNAM*10, TPNAM*10

      TPNAM = PNAM(1:SIZ)

      CALL MNPARM(ID, TPNAM, VSTRT, STP, LO, HI, IERFLG)
        IF (IERFLG .NE. 0) THEN
          WRITE(6, '(A,I)') ' Unable to define parametr no. ', ID
          STOP
        END IF

      RETURN
      END
C     ****

      SUBROUTINE MCOM(SIZ, CMD, IERFLG)

      IMPLICIT NONE
      INTEGER SIZ, IERFLG
      CHARACTER CMD*100, TCMD*100
      EXTERNAL FCN

      TCMD = CMD(1:SIZ)

      CALL MNCOMD(FCN, TCMD, IERFLG, 0)

      RETURN
      END
C     ****


      SUBROUTINE GETPAR(NUM, VAL, ERROR, IVARBL)

      IMPLICIT NONE
      INTEGER NUM, IVARBL
      DOUBLE PRECISION VAL, ERROR, BND1, BND2
      CHARACTER CHNAM*10


      CALL MNPOUT(NUM, CHNAM, VAL, ERROR, BND1, BND2, IVARBL)

      RETURN
      END
C     ****
