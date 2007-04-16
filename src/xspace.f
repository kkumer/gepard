C     ****h* gepard/xspace.f
C  FILE DESCRIPTION
C    calculation of x-space GPDs and PDFs
C
C  BUGS
C     *******


C     ****s* xspace.f/XSPACE
C  NAME
C     XSPACE  --   x-space GPDs and PDFs
C  DESCRIPTION
C    calculates x-space GPD or PDF by numerical Mellin-Barnes 
C    integration of conformal GPD moment, cf. Eq. (195) from the
c    big paper
C  SYNOPSIS
C     SUBROUTINE XSPACE
C  OUTPUT
C       XPDF  --  x PDF(x)
C  IDENTIFIERS
C           K -- Mellin-Barnes integration point index
C           J -- conformal moment
C  CHILDREN
C    HJ
C  SOURCE
C

      SUBROUTINE XSPACE (HX, X, T)

      IMPLICIT NONE
      DOUBLE PRECISION HX(2), X, T
      INTEGER K, L
      DOUBLE COMPLEX EPH, J, FCM(2)
      DOUBLE PRECISION RESIMAG
      DOUBLE PRECISION FIMAG
      INCLUDE 'header.f'

      DEL2 = T

      EPH = EXP ( COMPLEX(0.0d0, PHI) )

*   Integration ...

      DO 20 L = 1, 2
      HX(L) = 0.0d0
      DO 10 K = 1, NPTS

      J = N(1,K) - 1
      CALL HJ(J, FCM)

 10   HX(L) = HX(L) + WG(K) * IMAGPART( EPH * (1.0d0/X)**J * FCM(L) )
 20   HX(L) = HX(L) / PI

      RETURN
      END
C     ***
