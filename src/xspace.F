C     ****h* gepard/xspace.F
C  FILE DESCRIPTION
C    calculation of x-space GPDs and PDFs
C
C  BUGS
C     *******


C     ****s* xspace.F/XSPACE
C  NAME
C     XSPACE  --   x-space GPDs and PDFs
C  DESCRIPTION
C    calculates x-space GPD or PDF by numerical Mellin-Barnes 
C    integration of conformal GPD moment, cf. Eq. (195) from the
C    hep-ph/0703179
C  SYNOPSIS

      SUBROUTINE XSPACE (HX, X, T)

      IMPLICIT NONE
      DOUBLE PRECISION HX(2), X, T

C  INPUT
C           X  --  GPD parameter x
C           T  --  GPD parameter t=DEL2=\Delta^2
C  OUTPUT
C       HX(2)  --  x PDF(x)  (1: quark, 2: gluon)
C  CHILDREN
C    GETMBGPD
C  SOURCE
C

      INTEGER K, L
      DOUBLE COMPLEX EPH, J
#include "header.F"

      DEL2 = T

      EPH = EXP ( DBLCMPLX(0.0d0, PHI) )

      CALL GETMBGPD

*   Integration ...

      DO 20 L = 1, 2
      HX(L) = 0.0d0
      DO 10 K = 1, NPTS

      J = N(K) - 1

 10   HX(L) = HX(L) + WG(K) * IMAGPART(EPH * (1.0d0/X)**J * MBGPD(K,L))
 20   HX(L) = HX(L) / PI

      RETURN
      END
C     ***

C     ****s* xspace.F/SLOPEH
C  NAME
C     SLOPEH  --  Return slope B of GPD {\cal H}
C  DESCRIPTION
C     Calculates GPD slope, cf. Eq. (205) from hep-ph/0703179
C  INPUT
C      X -- GPD parameter x
C  OUTPUT
C      BH(2) -- GPD slope B, (1: quark, 2: gluon)
C  SYNOPSIS

      SUBROUTINE SLOPEH (BH, X)

      IMPLICIT NONE
      DOUBLE PRECISION BH(2), X
C
C  CHILDREN
C      LNHXQ, LNHXG, DERIV
C  SOURCE
C

      INTEGER NEVALS
      DOUBLE PRECISION LNHXQ, LNHXG, T, H, DRV, ERR
      EXTERNAL LNHXQ, LNHXG
      PARAMETER ( NEVALS = 5, H = 0.41d0 )
#include "header.F"

      PROCESS = 'DVCS'
      XI = X

      T = 0.0d0

      CALL DERIV(LNHXQ, T, H, NEVALS, DRV, ERR)
!      WRITE (30, *) XI, DRV, ERR
      BH(1) = DRV

      T = 0.0d0

      CALL DERIV(LNHXG, T, H, NEVALS, DRV, ERR)
!      WRITE (31, *) XI, DRV, ERR
      BH(2) = DRV

      RETURN
      END
C     ****



      DOUBLE PRECISION FUNCTION LNHXQ(T)

      IMPLICIT NONE
      DOUBLE PRECISION HX(2), T
#include "header.F"

      CALL XSPACE ( HX, XI, T)

      LNHXQ = LOG( HX(1) )

      RETURN
      END

      DOUBLE PRECISION FUNCTION LNHXG(T)

      IMPLICIT NONE
      DOUBLE PRECISION HX(2), T
#include "header.F"

      CALL XSPACE ( HX, XI, T)

      LNHXG = LOG( HX(2) )

      RETURN
      END

