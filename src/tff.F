C     ****h* gepard/tff.F
C  FILE DESCRIPTION
C    calculation of transition form factors for DVMP
C
C    $Id$
C     *******


C     ****s* tff.F/TFFF
C  NAME
C     TFFF  --  calculates flavour components of Transition form factor \mathcal{F}_V
C  DESCRIPTION
C    calculates singlet (sea and gluon) and non-singlet components of DVMP
C    transition form factors \mathcal{H}_M, \mathcal{E}_M, ... by numerical Mellin-Barnes 
C    integration according to formula Eq. (??) from [Mueller:2013??]
C      NB!! - only one of three TFFH(0:2) array elements, TFFH(P),
C             is changed by single call to TFFF!
C  SYNOPSIS

      SUBROUTINE TFFF

C  PARENTS
C    PARSIGMA, tffHInternal, tffEInternal
C  CHILDREN
C    GETMBGPD, LOOKUPQ, DCTAN
C  SOURCE
C


      IMPLICIT NONE
      INTEGER K, QIND
      DOUBLE COMPLEX EPH, J, CFACRE, CFACIM
      DOUBLE COMPLEX EFSIG, EFG, EFNS
      DOUBLE COMPLEX TFFHS, TFFHNS
      DOUBLE PRECISION REEFSIG, REEFG, REEFNS
      DOUBLE PRECISION IMEFSIG, IMEFG, IMEFNS
#include "header.F"
      DOUBLE COMPLEX DCTAN
* ASMUR2 = alpha_s/(2 pi) at renormalization scale
      DOUBLE PRECISION ASMUR2, ASTRONG
      CALL AS2PF (ASMUR2, Q2/RR2, ASP(P), MU02)
      ASTRONG = 2.d0 * PI * ASMUR2

      EPH = EXP ( DBLCMPLX(0.0d0, PHI) )

      REEFSIG = 0.d0
      IMEFSIG = 0.d0
      REEFG = 0.d0
      IMEFG = 0.d0
      REEFNS = 0.d0
      IMEFNS = 0.d0

      CALL LOOKUPQ(QIND)
*   Getting GPD values on MB contour ...

      CALL GETMBGPD

*   MB integration ...

      DO 123 K = 1, NPTS

        J = N(1,K) - 1
        CFACRE = EPH * DCTAN(PIHALF*J) / XI**(J+1.d0)
        CFACIM = EPH / XI**(J+1.d0)


*       --- flavour components ---

        EFSIG = CGRID(1,QIND,K,1)*MBGPD(K,1)
        EFG   = CGRID(1,QIND,K,2)*MBGPD(K,2)
        EFNS  = CGRIDNS(QIND,K)*MBGPD(K,1)
        
*           --- ading 2nd SO(3) partial wave ---
C       EFSIG = EFSIG + PAR(17) * CGRID(2,QIND,K,1) * MBGPD(K,1)
C       EFG   = EFG   + PAR(27) * CGRID(2,QIND,K,2) * MBGPD(K,2)

*           --- ading 3nd SO(3) partial wave ---
C       EFSIG = EFSIG + PAR(37) * CGRID(3,QIND,K,1) * MBGPD(K,1)
C       EFG   = EFG   + PAR(47) * CGRID(3,QIND,K,2) * MBGPD(K,2)


        REEFSIG = REEFSIG + WG(K)*IMAGPART(CFACRE * EFSIG)
        IMEFSIG = IMEFSIG + WG(K)*IMAGPART(CFACIM * EFSIG)

        REEFG = REEFG + WG(K)*IMAGPART(CFACRE * EFG)
        IMEFG = IMEFG + WG(K)*IMAGPART(CFACIM * EFG)

        REEFNS = REEFNS + WG(K)*IMAGPART(CFACRE * EFNS)
        IMEFNS = IMEFNS + WG(K)*IMAGPART(CFACIM * EFNS)

 123  CONTINUE
      
*     Building up complex TFFs with right prefactors

      TFFHS     =  DBLCMPLX(REEFSIG, IMEFSIG) / NF 
     &    + DBLCMPLX(REEFG, IMEFG) * 2.d0/CF/(J+3.d0)

      TFFHNS    =  DBLCMPLX(REEFNS, IMEFNS)

*     Combining with charge and CG factors for rho-meson production
*     (using model as in Dieter&Mantas paper)

      TFFH(P) = 3.0d0 * CF * FRHO * ASTRONG / NC / sqrt(Q2)
     & * ( TFFHS + 3*TFFHNS/20.d0 ) / SQRT(2.0d0)


      RETURN
      END
C     ***

