C     ****h* gepard/tff.F
C  FILE DESCRIPTION
C    calculation of transition form factors for DVMP
C
C    $Id$
C     *******


C     ****s* tff.F/TFFF
C  NAME
C     TFFF  --  calculates flavour components of Transition form factor \mathcal{H}_rho
C  DESCRIPTION
C    calculates singlet (sea and gluon) and non-singlet components of DVMP
C    transition form factors \mathcal{H}_M, \mathcal{E}_M, ... by numerical Mellin-Barnes 
C    integration according to formula Eq. (??) from [Mueller:2013??]
C      NB!! - only one of three TFFH(0:2) array elements, TFFH(P),
C             is changed by single call to TFFF!
C      NB - only CALH is done atm
C  SYNOPSIS

      SUBROUTINE TFFF

C  PARENTS
C    PARSIGMA, tffHInternal, tffEInternal
C  CHILDREN
C    GETMBGPD, LOOKUPQ, DCTAN
C  SOURCE
C


      IMPLICIT NONE
      INTEGER K, L, QIND, PIDK, SEC, NPWS
      DOUBLE COMPLEX HRHO(4), CH(4)
      DOUBLE COMPLEX EPH, J, CFACRE, CFACIM
      DOUBLE COMPLEX EFSIG, EFG, EFNS
      DOUBLE COMPLEX TFFHS, TFFHNS
      DOUBLE PRECISION RE(4), IM(4)
#include "header.F"
      DOUBLE PRECISION PWPAR(NPWMAX, 4)
      DOUBLE COMPLEX DCTAN
* ASMUR2 = alpha_s/(2 pi) at renormalization scale
      DOUBLE PRECISION ASMUR2, ASTRONG
      CALL AS2PF (ASMUR2, Q2/RR2, ASP(P), MU02)
      ASTRONG = 2.d0 * PI * ASMUR2

*     DVMP process ID = 3
      PIDK = 3
      EPH = EXP ( DBLCMPLX(0.0d0, PHI) )

      CALL LOOKUPQ(QIND)

*   Getting GPD values on MB contour ...

      CALL GETMBGPD

*   Setting the number and strength of SO(3) partial waves
      NPWS = 3
      PWPAR(1, 1) = 1.d0
      PWPAR(1, 2) = 1.d0
      PWPAR(2, 1) = PAR(17)
      PWPAR(2, 2) = PAR(27)
      PWPAR(3, 1) = PAR(37)
      PWPAR(3, 2) = PAR(47)

*   MB integration ...

      DO 5 L = 1, 4
        RE(L) = 0.d0
 5      IM(L) = 0.d0

      DO 30 K = 1, NPTS

        J = N(1,K) - 1
        CFACRE = EPH * DCTAN(PIHALF*J) / XI**(J+1.d0)
        CFACIM = EPH / XI**(J+1.d0)

*       --- rotation of flavour components ---

*       (SIG, G, NSP, NSM) = FROT * (sea, G, uval, dval)
c       FIXME: neglecting valence components atm
        HRHO(1) = GPDH(K,1) / SQRT2
        HRHO(2) = GPDH(K,2) / SQRT2
        HRHO(3) = GPDH(K,1) / SQRT2  * (3.d0/20.d0)
        HRHO(4) = 0.d0
        


c       FIXME: with NS(-), following loop has to go to 4
        DO 30 L = 1, 3
          CH(L) = (0.d0, 0.d0)
 10     CONTINUE
*         --- building evolved Wilson coef. summed over SO(3) partial waves ---
          DO 20 SEC = 1, NPWS
            CH(L) = CH(L) + PWPAR(SEC, L)*WCE(PIDK,SEC,QIND,K,L)
 20       CONTINUE
*         --- multiplying with GPD, to finalize that part of MB integrand
          CH(L) = CH(L) * HRHO(L)

          RE(L) = RE(L) + WG(K)*IMAGPART(CFACRE * CH(L))
          IM(L) = IM(L) + WG(K)*IMAGPART(CFACIM * CH(L))

 30   CONTINUE
      
*     Building up complex TFFs with right prefactors

      TFFH(P) = (0.d0, 0.d0)
*     FIXME: singlet-only 
      DO 40 L = 1, 2
        TFFH(P) = TFFH(P) + DBLCMPLX(RE(L), IM(L))
 40   CONTINUE

      TFFH(P) = CF * FRHO * ASTRONG / NC / sqrt(Q2) * TFFH(P)

      RETURN
      END
C     ***

