C     ****h* gepard/evolc.F
C  FILE DESCRIPTION
C    calculates "evolved" C's
C
C    $Id$
C     *******


C     ****f* evolc.F/GETQIND
C  NAME
C     GETQIND  --  return QIND for evolved Wilson coefs.
C  DESCRIPTION
C     It returns the QIND index where WCE array
C     will have the evolution corresponding to present Q2
C     and PID.
C     If Q2 is unseen yet, precalculate WCE first
C     for present PID.
C   
C  SYNOPSIS

      INTEGER FUNCTION GETQIND()

      IMPLICIT NONE
*     INTEGER PIDK

C  INPUT
C         PIDK -- process index
C  OUTPUT
C     GETQIND -- index of actual Q2 value in QS(PID, QIND) array
C  PARENTS
C      CFFF, TFFF
C  CHILDREN
C      EVOLC
C  SOURCE
C

      INTEGER K, SEC
#include "header.F"

      GETQIND = 0
      DO 200 K = 1, NQS(PID)
           IF ( QS(PID, K) .EQ. Q2 ) GETQIND = K
 200  CONTINUE

c     --- If Q2 not found, calculate WCE for this Q2
      IF (GETQIND .EQ. 0) THEN
        NQS(PID) = NQS(PID) + 1
        QS(PID, NQS(PID)) = Q2
        DO 210 SEC = 1, NPWMAX
 210      CALL EVOLC(SEC, NQS(PID))
        GETQIND = NQS(PID)
      END IF

      RETURN
      END
C     *****



C     ****s* evolc.F/EVOLC
C  NAME
C     EVOLC  --  "evolved" Wilson coefficients
C  DESCRIPTION
C      Calculates "evolved" Wilson coefficients
C      which means that Wilson coefficients are
C      multiplied by evolution operator.
C   
C  SYNOPSIS

      SUBROUTINE EVOLC (SEC, QIND)

      IMPLICIT NONE
      INTEGER SEC, QIND

C  INPUTS
C        SEC  -- SO(3) PW index
C        QIND -- index of actual Q2 value in QS array
C  PARENTS
C      GETQIND
C  CHILDREN
C      AS2PF, EVOLNSF, EVOLF
C  SOURCE
C

      INTEGER K
      DOUBLE PRECISION R, ASQ02, ASMUR2, ASMUF2
      DOUBLE COMPLEX CEVNSP, CEVNSM, CEV(2)
      DOUBLE COMPLEX EVOLNSA(0:1,0:2)
      DOUBLE COMPLEX EVOLA(0:2,2,2)
#include "header.F"

* ASMUF2 = alpha_s/(2 pi) at factorization scale
      CALL AS2PF (ASMUF2, Q2/RF2, ASP(P), MU02)
* ASMUR2 = alpha_s/(2 pi) at renormalization scale
      CALL AS2PF (ASMUR2, Q2/RR2, ASP(P), MU02)
* ASQ02 = alpha_s/(2 pi) at input scale
      CALL AS2PF (ASQ02, Q02, ASP(P), MU02)

      R = ASMUF2/ASQ02

      DO 123 K = 1, NPTS


**    1. non-singlet

        CALL EVOLNSF (SEC, K, R, EVOLNSA)

*       CZERO puts astrong^0 term to zero for investigation of NLO effects

        CEVNSP = WC(PID, SEC, K, 0, 3) * EVOLNSA(0,0) * CZERO
        CEVNSM = WC(PID, SEC, K, 0, 4) * EVOLNSA(0,0) * CZERO

        IF (P .GE. 1) THEN

         CEVNSP = CEVNSP + ASMUR2 * WC(PID, SEC, K, 1, 3) * EVOLNSA(0,0)
     &                  + ASMUF2 * WC(PID, SEC, K, 0, 3) * EVOLNSA(0,1) 
         CEVNSM = CEVNSM + ASMUR2 * WC(PID, SEC, K, 1, 4) * EVOLNSA(1,0)
     &                  + ASMUF2 * WC(PID, SEC, K, 0, 4) * EVOLNSA(1,1) 
*       FIXME: Wilson coef. might also have NS(-) part different from
*       NS(+) at NLO !

          IF (P .GE. 2) THEN

          CEVNSP = CEVNSP + ASMUR2**2 * WC(PID,SEC,K,2,3)*EVOLNSA(0,0) +
     &          ASMUR2 * ASMUF2 * WC(PID,SEC,K,1,3) * EVOLNSA(0,1) + 
     &                    ASMUF2**2 *  WC(PID,SEC,K,0,3) * EVOLNSA(0,2)
          CEVNSM = CEVNSM + ASMUR2**2 * WC(PID,SEC,K,2,4)*EVOLNSA(1,0) +
     &          ASMUR2 * ASMUF2 * WC(PID,SEC,K,1,4) * EVOLNSA(1,1) + 
     &                    ASMUF2**2 *  WC(PID,SEC,K,0,4) * EVOLNSA(1,2)
*       FIXME: Implement NNLO NS(-) evolution
*       FIXME: Wilson coef. might also have NS(-) part different from
*       NS(+) at NLO and NNLO !
          END IF

        END IF

        WCE(PID, SEC, QIND, K, 3) = CEVNSP
        WCE(PID, SEC, QIND, K, 4) = CEVNSM

*     2. singlet case

       CALL EVOLF (SEC, K, R, EVOLA)

        CEV(1) = ( WC(PID,SEC,K,0,1) * EVOLA(0,1,1) + 
     &             WC(PID,SEC,K,0,2) * EVOLA(0,2,1) ) * CZERO
        CEV(2) = ( WC(PID,SEC,K,0,1) * EVOLA(0,1,2) + 
     &             WC(PID,SEC,K,0,2) * EVOLA(0,2,2) ) * CZERO

        IF (P .GE. 1) THEN

          CEV(1)=CEV(1) + ASMUR2 * ( WC(PID,SEC,K,1,1) * EVOLA(0,1,1) + 
     &                              WC(PID,SEC,K,1,2) * EVOLA(0,2,1) )
     &                 + ASMUF2 * ( WC(PID,SEC,K,0,1) * EVOLA(1,1,1) + 
     &                              WC(PID,SEC,K,0,2) * EVOLA(1,2,1) )

          CEV(2)=CEV(2) + ASMUR2 * ( WC(PID,SEC,K,1,1) * EVOLA(0,1,2) + 
     &                               WC(PID,SEC,K,1,2) * EVOLA(0,2,2) )
     &                  + ASMUF2 * ( WC(PID,SEC,K,0,1) * EVOLA(1,1,2) + 
     &                               WC(PID,SEC,K,0,2) * EVOLA(1,2,2) )


          IF (P .GE. 2) THEN

            CEV(1)=CEV(1)+ASMUR2**2*( WC(PID,SEC,K,2,1) * EVOLA(0,1,1) +
     &                               WC(PID,SEC,K,2,2) * EVOLA(0,2,1) )
     &        + ASMUR2 * ASMUF2 * ( WC(PID,SEC,K,1,1) * EVOLA(1,1,1) + 
     &                              WC(PID,SEC,K,1,2) * EVOLA(1,2,1) )
     &               + ASMUF2**2 * ( WC(PID,SEC,K,0,1) * EVOLA(2,1,1) + 
     &                               WC(PID,SEC,K,0,2) * EVOLA(2,2,1) )

            CEV(2)=CEV(2)+ASMUR2**2*( WC(PID,SEC,K,2,1) * EVOLA(0,1,2) +
     &                               WC(PID,SEC,K,2,2) * EVOLA(0,2,2) )
     &        + ASMUR2 * ASMUF2 * ( WC(PID,SEC,K,1,1) * EVOLA(1,1,2) + 
     &                              WC(PID,SEC,K,1,2) * EVOLA(1,2,2) )
     &               + ASMUF2**2 * ( WC(PID,SEC,K,0,1) * EVOLA(2,1,2) + 
     &                               WC(PID,SEC,K,0,2) * EVOLA(2,2,2) )

          ENDIF

        ENDIF

        WCE(PID, SEC, QIND, K, 1) = CEV(1)
        WCE(PID, SEC, QIND, K, 2) = CEV(2)


 123  CONTINUE

      RETURN
      END
C     *****
      
