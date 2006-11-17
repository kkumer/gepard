C     ****h* gepard/error.f
C  FILE DESCRIPTION
C    treating errors
C
C    $Id: erfunc.f 11 2006-07-12 07:50:28Z kuk05260 $
C     *******


C     ****s* error.f/ERROR
C  NAME
C     ERROR  --   produces error messages
C  SYNOPSIS
C     SUBROUTINE ERROR (LIB, SUBR, MSG, NERR, LEVEL)
C
C     IMPLICIT NONE
C     INTEGER NERR, LEVEL
C     CHARACTER LIB*6, SUBR*5, MSG*60
C  INPUTS
C         LIB -- name of the library
C        SUBR -- name of the subroutine
C         MSG -- message itself
C        NERR -- unique error number
C       LEVEL -- severity level of error
C  PARENTS
C      PROCDATA, DCTAN
C  SOURCE
C
      SUBROUTINE ERROR (LIB, SUBR, MSG, NERR, LEVEL)

      IMPLICIT NONE
      INTEGER NERR, LEVEL
      CHARACTER LIB*6, SUBR*8, MSG*60

      WRITE (*,*) LIB, ": ", SUBR, " ERROR !!", MSG
      WRITE (*,*) " --> ERROR NUMBER: ", NERR
      WRITE (*,*) " --> ERROR LEVEL:  ", LEVEL

      STOP

      RETURN
      END
C     *****
