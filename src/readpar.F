C     ****h* gepard/readpar.F
C  FILE DESCRIPTION
C    initial reading of parameters
C
C    $Id$
C     *******

C     ****s* readpar.f/READPAR
C  NAME
C        READPAR -- reads parameters from file GEPARD.INI
C  DESCRIPTION
C     (See gepard.pdf or provided self-documented GEPARD.INI for meaning
C     of each parameter.)
C  SYNOPSIS

      SUBROUTINE READPAR

C  PARENTS
C      FIT
C  SOURCE
C

      IMPLICIT NONE
#include "header.F"


*   Read fixed initialization parameters
      OPEN (UNIT = 61, FILE = "GEPARD.INI", STATUS = "OLD")
      READ (61, *) SPEED
      READ (61, *) ACC
      READ (61, *) P
      READ (61, *) NF
      READ (61, *) CZERO
      READ (61, *) MU02
      READ (61, *) ASP(0)
      READ (61, *) ASP(1)
      READ (61, *) ASP(2)
      READ (61, *) Q02
      READ (61, *) RF2
      READ (61, *) RR2
      READ (61, *) C
      READ (61, *) PHI
      READ (61, *) CND
      READ (61, *) PHIND
      READ (61, *) SCHEME
      READ (61, *) ANSATZ
      CLOSE (61)

*   default filenames
      DATFILE = 'fit.dat'
      OUTFILE = 'fit.out'
      CMDFILE = 'fit.cmd'

      RETURN
      END
C     ***
