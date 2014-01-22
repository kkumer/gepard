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
C      NB!! - only one of three CFF(0:2) array elements, CFF(P),
C             is changed by single call to CFFF!
C  SYNOPSIS

      SUBROUTINE CFFF

C  PARENTS
C    PARSIGMA, cffHInternal, cffEInternal
C  CHILDREN
C    GETMBGPD, GETMBGPDMMA, GETQIND, DCTAN
C  SOURCE
C


      IMPLICIT NONE
      INTEGER K, QIND, GETQIND
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
      IF ((PID .EQ. -1) .OR. (PID .EQ. -2)) XI = 0

      QIND = GETQIND()
      CALL GETMBGPD

      XI = XIMEM

*   Integration ...

      DO 123 K = 1, NPTS

      J = N(1,K) - 1
      CFAC = EPH / XI**(J+1.d0)

      IF ( FFTYPE(:7) .EQ. 'NONSING' ) THEN
        IF (ANSATZ .EQ. 'NSHOUC') THEN
*              NS(-)
          FPWH = WCE(PID,1,QIND,K,4)*GPDH(K,1)        
        ELSE
*              NS(+)
          FPWH = WCE(PID,1,QIND,K,3)*GPDH(K,1)        
        END IF
*        FIXME: Following line could be wrong parameter no. It's ther
*               just to avoid compiler warning ... NONSING is maybe
*               broken!
        FPWE = WCE(PID,1,QIND,K,3)*PAR(17)*GPDE(K,1)        
      ELSE
*       --- singlet case ---
*            --- CFF H ---
        FPWH = WCE(PID,1,QIND,K,1)*GPDH(K,1) + 
     &         WCE(PID,1,QIND,K,2)*GPDH(K,2)
        FPWH = FPWH +  
     &              PAR(17) * WCE(PID,2,QIND,K,1) * GPDH(K,1) + 
     &              PAR(27) * WCE(PID,2,QIND,K,2) * GPDH(K,2)
        FPWH = FPWH +  
     &              PAR(37) * WCE(PID,3,QIND,K,1) * GPDH(K,1) + 
     &              PAR(47) * WCE(PID,3,QIND,K,2) * GPDH(K,2)
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

        FPWE = WCE(PID,1,QIND,K,1)*PAR(18)*GPDE(K,1) 
     &       + WCE(PID,1,QIND,K,2)*PAR(28)*GPDE(K,2)
        FPWE = FPWE +  
     &          PAR(18) * PAR(117) * WCE(PID,2,QIND,K,1) * GPDE(K,1) + 
     &          PAR(28) * PAR(127) * WCE(PID,2,QIND,K,2) * GPDE(K,2)
        FPWE = FPWE +  
     &          PAR(18) * PAR(137) * WCE(PID,3,QIND,K,1) * GPDE(K,1) + 
     &          PAR(28) * PAR(147) * WCE(PID,3,QIND,K,2) * GPDE(K,2)
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

      IF ( FFTYPE(:7) .EQ. 'SINGLET' ) THEN
        IF (NF .EQ. 3) THEN
           CHARGEFAC =  2.0D0 / 9.0D0
        ELSE IF (NF .EQ. 4) THEN
           CHARGEFAC =  5.0D0 / 18.0D0
        ELSE
           CALL ERROR ('GeParD', ' INIT',
     & 'NF is not integer equal to 3 or 4!                          ',
     & 5, 3)
        END IF
      ELSE
*         -- NONSINGLET --
        IF (NF .EQ. 3) THEN
           CHARGEFAC =  1.0D0 / 9.0D0
        ELSE IF (NF .EQ. 4) THEN
           CHARGEFAC =  1.0D0 / 6.0D0
        ELSE
           CALL ERROR ('GeParD', ' INIT',
     & 'NF is not integer equal to 3 or 4!                          ',
     & 5, 3)
        END IF
      END IF

*         -- when calculating just GPDs --

      IF (PID .LE. -1) CHARGEFAC = 1.0d0

      CFF(P) = CHARGEFAC * DBLCMPLX(RESREALH, RESIMAGH)
      CFFE(P) = CHARGEFAC * DBLCMPLX(RESREALE, RESIMAGE) 

      RETURN
      END
C     ***

