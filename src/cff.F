C     ****h* gepard/cff.F
C  FILE DESCRIPTION
C    calculation of Compton form factors
C
C    $Id$
C     *******


C     ****s* cff.F/CFFF
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
      DOUBLE PRECISION RESREALH, RESIMAGH
      DOUBLE PRECISION FREALH, FIMAGH
      DOUBLE PRECISION RESREALE, RESIMAGE
      DOUBLE PRECISION FREALE, FIMAGE
      DOUBLE PRECISION XIMEM
#include "header.F"
      DOUBLE COMPLEX DCTAN

      EPH = EXP ( DBLCMPLX(0.0d0, PHI) )

      RESREALH = 0.0d0
      RESIMAGH = 0.0d0
      RESREALE = 0.0d0
      RESIMAGE = 0.0d0


      XIMEM = XI
*     -- when calculating just GPDs at eta=0 --
      IF ( PROCESS(5:5) .EQ. 'Z' ) XI = 0

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

      J = N(1,K) - 1
      CFAC = EPH / XI**(J+1.d0)

      IF ( FFTYPE(:7) .EQ. 'NONSING' ) THEN
        FPWH = CGRIDNS(QIND,K)*MBGPD(K,1)        
*        FIXME: Following line could be wrong parameter no. It's ther
*               just to avoid compiler warning ... NONSING is maybe
*               broken!
        FPWE = CGRIDNS(QIND,K)*PAR(18)*MBGPD(K,3)        
      ELSE
*       --- singlet case ---
*            --- CFF H ---
        FPWH = CGRID(1,QIND,K,1)*MBGPD(K,1) + 
     &         CGRID(1,QIND,K,2)*MBGPD(K,2)
        FPWH = FPWH +  
     &              PAR(17) * CGRID(2,QIND,K,1) * MBGPD(K,1) + 
     &              PAR(27) * CGRID(2,QIND,K,2) * MBGPD(K,2)
        FPWH = FPWH +  
     &              PAR(37) * CGRID(3,QIND,K,1) * MBGPD(K,1) + 
     &              PAR(47) * CGRID(3,QIND,K,2) * MBGPD(K,2)
*            --- CFF E ---
        IF (.NOT. ( (ANSATZ .EQ. 'EFL') .OR. (ANSATZ .EQ. 'EFLEXP')
     &           )) THEN
*         we set strength of subleading PWs for E equal to those of H
*         to get full proportionality  E \propto H
          PAR(117) = PAR(17)
          PAR(127) = PAR(27)
          PAR(137) = PAR(37)
          PAR(147) = PAR(47)
        END IF
        IF (ANSATZ(:1) .EQ. 'E') THEN
*         we implement BS + BG = 0 sum rule
          PAR(28) = - PAR(18)*PAR(11) / (0.6d0 - PAR(11))
        END IF

        FPWE = CGRID(1,QIND,K,1)*PAR(18)*MBGPD(K,3) 
     &       + CGRID(1,QIND,K,2)*PAR(28)*MBGPD(K,4)
        FPWE = FPWE +  
     &            PAR(18) * PAR(117) * CGRID(2,QIND,K,1) * MBGPD(K,3) + 
     &            PAR(28) * PAR(127) * CGRID(2,QIND,K,2) * MBGPD(K,4)
        FPWE = FPWE +  
     &            PAR(18) * PAR(137) * CGRID(3,QIND,K,1) * MBGPD(K,3) + 
     &            PAR(28) * PAR(147) * CGRID(3,QIND,K,2) * MBGPD(K,4)
      END IF

      IF (ANSATZ(:1) .EQ. 'E') THEN
*     PAR(38)=Delta-B for CFF-E as proposed by Markus
        FPWE = FPWE * EXP(PAR(48)*DEL2)
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

