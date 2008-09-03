C     ****h* gepard/cff.f
C  FILE DESCRIPTION
C    calculation of Compton form factors
C
C    $Id$
C     *******


C     ****s* cff.f/CFFF
C  NAME
C     CFFF  --   Compton form factor \mathcal{H}
C  DESCRIPTION
C    calculates singlet Compton form factor H by numerical Mellin-Barnes 
C    integration according to formula Eq. (58) from [Kumericki:2007sa]
C      NB!! - only one of three CFF(0:2) array elements, CFF(P)
C             is changed by single call to CFFF!
C  SYNOPSIS

      SUBROUTINE CFFF

C  PARENTS
C    PARSIGMA, cffHInternal, cffEInternal
C  CHILDREN
C    GETMBGPD, GETMBGPDMMA, LOOKUPQ, DCTAN
C  SOURCE
C


      IMPLICIT NONE
      INTEGER K, QIND
      DOUBLE COMPLEX EPH, J, CFAC, FPWH, FPWE
      DOUBLE COMPLEX DCTAN
      DOUBLE PRECISION RESREALH, RESIMAGH
      DOUBLE PRECISION FREALH, FIMAGH
      DOUBLE PRECISION RESREALE, RESIMAGE
      DOUBLE PRECISION FREALE, FIMAGE
      DOUBLE PRECISION XIMEM
#include "header.F"

      EPH = EXP ( DBLCMPLX(0.0d0, PHI) )

      RESREALH = 0.0d0
      RESIMAGH = 0.0d0
      RESREALE = 0.0d0
      RESIMAGE = 0.0d0


      XIMEM = XI
*     -- when calculating just GPDs at eta=0 --
      IF ( SCHEME(:4) .EQ. 'ZERO' ) XI = 0

      CALL LOOKUPQ(QIND)
#ifdef MMA
      IF (ANSATZ(:3) .EQ. 'MMA') THEN
        CALL GETMBGPDMMA
      ELSE
        CALL GETMBGPD
      END IF
#else
      IF (ANSATZ(:3) .EQ. 'MMA') THEN
        WRITE(*,*) 'ANSATZ=MMA unavailable. Recompile with -D MMA'
      ELSE
        CALL GETMBGPD
      END IF
#endif

      XI = XIMEM

*   Integration ...

      DO 123 K = 1, NPTS

      J = N(K) - 1
      CFAC = EPH / XI**(J+1.d0)

      IF ( FFTYPE(:7) .EQ. 'NONSING' ) THEN
        FPWH = CGRIDNS(QIND,K)*MBGPD(K,1)        
      ELSE
        FPWH = CGRID(QIND,K,1)*MBGPD(K,1) + CGRID(QIND,K,2)*MBGPD(K,2)
        FPWE = CGRID(QIND,K,1)*PAR(18)*MBGPD(K,1) 
     &       + CGRID(QIND,K,2)*PAR(28)*MBGPD(K,2)
      END IF

      FREALH = IMAGPART(CFAC * DCTAN(PIHALF*J) * FPWH)
      FREALE = IMAGPART(CFAC * DCTAN(PIHALF*J) * FPWE)

      FIMAGH = IMAGPART(CFAC * FPWH)
      FIMAGE = IMAGPART(CFAC * FPWE)

      RESREALH = RESREALH + WG(K)*FREALH
      RESIMAGH = RESIMAGH + WG(K)*FIMAGH
      
      RESREALE = RESREALE + WG(K)*FREALE
      RESIMAGE = RESIMAGE + WG(K)*FIMAGE

 123  CONTINUE
      
*     Multiplying with charge factor

      CFF(P) = CHARGEFAC * DBLCMPLX(RESREALH, RESIMAGH)
      CFFE(P) = CHARGEFAC * DBLCMPLX(RESREALE, RESIMAGE) 

      RETURN
      END
C     ***

