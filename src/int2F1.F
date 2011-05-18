
*  Just a wrapper for subroutine below. Subroutine is
*   used for direct ML call, and function from cdvem.F

      DOUBLE COMPLEX FUNCTION FINT(MU, J, NU, K, IND)

      IMPLICIT NONE
      REAL MU, NU, IND, ACC
      DOUBLE COMPLEX J, K

      ACC = 4.
      CALL INT2F1F(MU, J, NU, K, IND, ACC, FINT)

      RETURN
      END


      SUBROUTINE INT2F1F(MU, J, NU, K, IND, ACC, FINT)

      IMPLICIT NONE
      REAL MU, NU, IND, ACC
      DOUBLE COMPLEX J, K
      DOUBLE COMPLEX FINT

      INTEGER EFFACC, NBNDMAX
      PARAMETER (NBNDMAX = 1)
      INTEGER MTINDMAX, L, NMTS, LI
      PARAMETER (MTINDMAX = 100)
      DOUBLE PRECISION MTS(MTINDMAX)
      DOUBLE PRECISION MTWG(MTINDMAX)
      DOUBLE PRECISION BND(NBNDMAX+1)
      DOUBLE COMPLEX INNER, OUTER
      DOUBLE PRECISION X, Z, PREF
      DOUBLE COMPLEX  INTGR, HYP_2F1
      DOUBLE COMPLEX ac, bc, cc, aln, bln, cln, zln
      EXTERNAL HYP_2F1

*     Integration is done on the interval:
      DATA BND /0.0D0, 1.0D0/ 
*     Effective numerical accuracy used for this integration is:
      EFFACC = NINT(ACC)
*     which means that resulting number of integration points is:
      NMTS = 2**EFFACC
*     Calculating abscissas and weights
      CALL INTEGRAF(EFFACC, 1, BND, NBNDMAX, MTS, MTWG)

*     PRINT *, "ARGS = ", MU, J, NU, K, IND, EFFACC
      aln = 1.0d0
      bln = 1.0d0
      cln = J + 3.5d0 - MU
      ac = K + 2.5d0 - NU
      bc = K + 2.0d0
      cc = 2.0d0*K+4.0d0
*     WRITE(*,*) "IND, NINT(IND) = ", IND, NINT(IND)

      OUTER = 0.0d0
      DO 20 L = 1, NMTS
        Z = MTS(L)
        INNER = 0.0d0
        IF (NINT(IND) .EQ. 1) THEN
            PREF = 1.0d0
        ELSE IF (NINT(IND) .EQ. 2) THEN
            PREF = (2.0d0-Z)*Z
        ELSE IF (NINT(IND) .EQ. 3) THEN
            PREF = (1.0d0-Z)/Z
        ELSE IF (NINT(IND) .EQ. 4) THEN
            PREF = Z - 2.0d0
        ELSE IF (NINT(IND) .EQ. 5) THEN
            PREF = Z
        ELSE
            PREF = 0.0d0
        END IF
*       WRITE(*,*) 'PREF = ', NINT(IND), L, Z, PREF
        DO 10 LI = 1, NMTS
          X = MTS(LI)
          INTGR = PREF * X**(MU-0.5) * (1.-X)**(1.0d0 + J) *
     &             Z**(K+1.5-NU)
          zln = 1.0d0 - X*Z
          INTGR = INTGR *
     &      HYP_2F1(aln, bln, cln, zln) * 
     &      HYP_2F1(ac, bc, cc, Z)
          INNER = INNER + MTWG(LI) * INTGR
 10     CONTINUE
        OUTER = OUTER + MTWG(L) * INNER
 20   CONTINUE

      FINT = OUTER
*     WRITE (*, *)  'INTEGRAL = ', FINT


      RETURN
      END 

