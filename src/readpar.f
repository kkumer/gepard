C     ****h* gepard/readpar.f
C  FILE DESCRIPTION
C    initial reading of parameters
C
C    $Id: dvcs.f 11 2006-07-12 07:50:28Z kuk05260 $
C     *******

C     ****s* readpar.f/READPAR
C  NAME
C        READPAR -- reads parameters from file GEPARD.INI
C  SYNOPSIS
C     SUBROUTINE READPAR
C  IDENTIFIERS
C     (see gepard.pdf or provided self-documented GEPARD.INI)
C  PARENTS
C      AUXTEST, TEST, RADCORR, SCALEDEP, FIT
C  SOURCE
C

      SUBROUTINE READPAR

      IMPLICIT NONE
      INCLUDE 'header.f'


*   Read fixed initialization parameters
      OPEN (UNIT = 61, FILE = "GEPARD.INI", STATUS = "OLD")
      READ (61, *) SPEED
      READ (61, *) ACC
      READ (61, *) P
      READ (61, *) NF
      READ (61, *) CZERO
      READ (61, *) RF2
      READ (61, *) RR2
      READ (61, *) C
      READ (61, *) PHI
      READ (61, *) CND
      READ (61, *) PHIND
      READ (61, *) SCHEME
      READ (61, *) ANSATZ
      READ (61, *) PROCESS
      CLOSE (61)

      RETURN
      END
C     ***
