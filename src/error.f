      SUBROUTINE ERROR (LIB, SUBR, MSG, NERR, LEVEL)

      IMPLICIT NONE
      INTEGER NERR, LEVEL
      CHARACTER LIB*6, SUBR*5, MSG*60

      WRITE (*,*) LIB, ": ", SUBR, " ERROR !!", MSG
      WRITE (*,*) " --> ERROR NUMBER: ", NERR
      WRITE (*,*) " --> ERROR LEVEL:  ", LEVEL

      STOP

      RETURN
      END
