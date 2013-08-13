C     ****h* gepard/evolc.F
C  FILE DESCRIPTION
C    calculates "evolved" C's
C
C    $Id$
C     *******

C     ****s* evolc.F/EVOLC
C  NAME
C     EVOLC  --  "evolved" Wilson coefficients
C  DESCRIPTION
C      Calculates "evolved" Wilson coefficients  for both singlet and
C      nonsinglet DVCS and DIS, which means that Wilson coefficients are
C      multiplied by evolution operator.
C   
C  SYNOPSIS

      SUBROUTINE EVOLC (SEC, QIND )

      IMPLICIT NONE
      INTEGER QIND, SEC

C  INPUTS
C        QIND -- index of actual Q2 value in QS or QSDIS array
C  PARENTS
C      INITC
C  CHILDREN
C      AS2PF, EVOLF
C  SOURCE
C

      INTEGER K
      DOUBLE COMPLEX J
      DOUBLE PRECISION R, ASQ02, ASMUR2, ASMUF2
      DOUBLE COMPLEX CEVNS, CEV(2), CEVDIS(2)
      DOUBLE COMPLEX EVOLNSA(0:2)
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

        CALL EVOLNSF (K, R, EVOLNSA)

*       CZERO puts astrong^0 term to zero for investigation of NLO effects

        CEVNS = BIGCNS(K, 0) * EVOLNSA(0) * CZERO

        IF (P .GE. 1) THEN

          CEVNS = CEVNS + ASMUR2 * BIGCNS(K, 1) * EVOLNSA(0) +
     &                    ASMUF2 * BIGCNS(K, 0) * EVOLNSA(1) 

          IF (P .GE. 2) THEN

            CEVNS = CEVNS + ASMUR2**2 * BIGCNS(K,2) * EVOLNSA(0) +
     &          ASMUR2 * ASMUF2 * BIGCNS(K,1) * EVOLNSA(1) + 
     &                      ASMUF2**2 *  BIGCNS(K,0) * EVOLNSA(2)
          END IF

        END IF

        CGRIDNS(QIND, K) = CEVNS

*     2. singlet case

       CALL EVOLF (SEC, K, R, EVOLA)

*       2.1  --- DVCS ---

        CEV(1) = ( BIGC(SEC,K,0,1) * EVOLA(0,1,1) + 
     &             BIGC(SEC,K,0,2) * EVOLA(0,2,1) ) * CZERO
        CEV(2) = ( BIGC(SEC,K,0,1) * EVOLA(0,1,2) + 
     &             BIGC(SEC,K,0,2) * EVOLA(0,2,2) ) * CZERO

        IF (P .GE. 1) THEN

          CEV(1) = CEV(1) + ASMUR2 * ( BIGC(SEC,K,1,1) * EVOLA(0,1,1) + 
     &                                 BIGC(SEC,K,1,2) * EVOLA(0,2,1) )
     &                    + ASMUF2 * ( BIGC(SEC,K,0,1) * EVOLA(1,1,1) + 
     &                                 BIGC(SEC,K,0,2) * EVOLA(1,2,1) )

          CEV(2) = CEV(2) + ASMUR2 * ( BIGC(SEC,K,1,1) * EVOLA(0,1,2) + 
     &                                 BIGC(SEC,K,1,2) * EVOLA(0,2,2) )
     &                    + ASMUF2 * ( BIGC(SEC,K,0,1) * EVOLA(1,1,2) + 
     &                                 BIGC(SEC,K,0,2) * EVOLA(1,2,2) )


          IF (P .GE. 2) THEN

            CEV(1) = CEV(1)+ASMUR2**2*( BIGC(SEC,K,2,1) * EVOLA(0,1,1) +
     &                                 BIGC(SEC,K,2,2) * EVOLA(0,2,1) )
     &        + ASMUR2 * ASMUF2 * ( BIGC(SEC,K,1,1) * EVOLA(1,1,1) + 
     &                              BIGC(SEC,K,1,2) * EVOLA(1,2,1) )
     &                 + ASMUF2**2 * ( BIGC(SEC,K,0,1) * EVOLA(2,1,1) + 
     &                                 BIGC(SEC,K,0,2) * EVOLA(2,2,1) )

            CEV(2) = CEV(2)+ASMUR2**2*( BIGC(SEC,K,2,1) * EVOLA(0,1,2) +
     &                                 BIGC(SEC,K,2,2) * EVOLA(0,2,2) )
     &        + ASMUR2 * ASMUF2 * ( BIGC(SEC,K,1,1) * EVOLA(1,1,2) + 
     &                              BIGC(SEC,K,1,2) * EVOLA(1,2,2) )
     &                 + ASMUF2**2 * ( BIGC(SEC,K,0,1) * EVOLA(2,1,2) + 
     &                                 BIGC(SEC,K,0,2) * EVOLA(2,2,2) )

          ENDIF

        ENDIF

        IF (PROCESS(:4) .EQ. 'DVMP') THEN
          J = N(SEC, K) - 1
          CGRID(SEC, QIND, K, 1) = CEV(1) / NF
          CGRID(SEC, QIND, K, 2) = CEV(2) * 2.d0/CF/(J+3.d0)
        ELSE
          CGRID(SEC, QIND, K, 1) = CEV(1)
          CGRID(SEC, QIND, K, 2) = CEV(2)
        END IF

*        2.2 -- DIS ---

        CEVDIS(1) = ( BIGCF2(K,0,1) * EVOLA(0,1,1) + 
     &             BIGCF2(K,0,2) * EVOLA(0,2,1) ) * CZERO
        CEVDIS(2) = ( BIGCF2(K,0,1) * EVOLA(0,1,2) + 
     &             BIGCF2(K,0,2) * EVOLA(0,2,2) ) * CZERO

        IF (P .GE. 1) THEN

          CEVDIS(1) = CEVDIS(1) + ASMUR2 * ( BIGCF2(K,1,1) * 
     &                  EVOLA(0,1,1) +  BIGCF2(K,1,2) * EVOLA(0,2,1) )
     &                    + ASMUF2 * ( BIGCF2(K,0,1) * EVOLA(1,1,1) + 
     &                                 BIGCF2(K,0,2) * EVOLA(1,2,1) )

          CEVDIS(2) = CEVDIS(2) + ASMUR2 * ( BIGCF2(K,1,1) * 
     &                 EVOLA(0,1,2) +  BIGCF2(K,1,2) * EVOLA(0,2,2) )
     &                    + ASMUF2 * ( BIGCF2(K,0,1) * EVOLA(1,1,2) + 
     &                                 BIGCF2(K,0,2) * EVOLA(1,2,2) )


          IF (P .GE. 2) THEN

          CEVDIS(1) = CEVDIS(1) + ASMUR2**2 * ( BIGCF2(K,2,1) * 
     &                    EVOLA(0,1,1) + BIGCF2(K,2,2) * EVOLA(0,2,1) )
     &        + ASMUR2 * ASMUF2 * ( BIGCF2(K,1,1) * EVOLA(1,1,1) + 
     &                              BIGCF2(K,1,2) * EVOLA(1,2,1) )
     &                  + ASMUF2**2 * ( BIGCF2(K,0,1) * EVOLA(2,1,1) + 
     &                                    BIGCF2(K,0,2) * EVOLA(2,2,1) )

          CEVDIS(2) = CEVDIS(2) + ASMUR2**2 * ( BIGCF2(K,2,1) * 
     &                    EVOLA(0,1,2) + BIGCF2(K,2,2) * EVOLA(0,2,2) )
     &        + ASMUR2 * ASMUF2 * ( BIGCF2(K,1,1) * EVOLA(1,1,2) + 
     &                              BIGCF2(K,1,2) * EVOLA(1,2,2) )
     &                  + ASMUF2**2 * ( BIGCF2(K,0,1) * EVOLA(2,1,2) + 
     &                                    BIGCF2(K,0,2) * EVOLA(2,2,2) )

          ENDIF

        ENDIF

        CGRIDDIS(QIND, K, 1) = CEVDIS(1)
        CGRIDDIS(QIND, K, 2) = CEVDIS(2)


 123  CONTINUE

      RETURN
      END
C     *****
