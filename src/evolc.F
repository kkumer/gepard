C     ****h* gepard/evolc.f
C  FILE DESCRIPTION
C    calculates "evolved" C's
C
C    $Id: parwav.f 71 2007-04-18 20:39:48Z kkumer $
C     *******

C     ****s* evolc.f/EVOLC
C  NAME
C     EVOLC  --  "evolved" Wilson coefficients
C  DESCRIPTION
C      Calculates "evolved" Wilson coefficients  for both singlet and
C      nonsinglet DVCS and DIS, which means that Wilson coefficients are
C      multiplied by evolution operator.
C   
C  SYNOPSIS

      SUBROUTINE EVOLC ( QIND )

      IMPLICIT NONE
      INTEGER QIND

C  INPUTS
C        QIND -- index of actual Q2 value in QS or QSDIS array
C  PARENTS
C      INITC
C  CHILDREN
C      AS2PF, EVOLF
C  SOURCE
C

      INTEGER K
      DOUBLE PRECISION R, ASQ02, ASMUR2, ASMUF2
      DOUBLE COMPLEX CEVNS, CEV(2)
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

*        Evolution operators and partial wave depending on whether we do
*        non-singlet or singlet case. This is specified by using
*        process name that starts 'NS...'.

      IF ( FFTYPE(:7) .EQ. 'NONSING' ) THEN

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

      ELSE
*     singlet case

       CALL EVOLF ( K, R, EVOLA)

       IF ( PROCESS(:3) .EQ. 'DVC' ) THEN
*      --- DVCS ---

        CEV(1) = ( BIGC(K,0,1) * EVOLA(0,1,1) + 
     &             BIGC(K,0,2) * EVOLA(0,2,1) ) * CZERO
        CEV(2) = ( BIGC(K,0,1) * EVOLA(0,1,2) + 
     &             BIGC(K,0,2) * EVOLA(0,2,2) ) * CZERO

        IF (P .GE. 1) THEN

          CEV(1) = CEV(1) + ASMUR2 * ( BIGC(K,1,1) * EVOLA(0,1,1) + 
     &                                 BIGC(K,1,2) * EVOLA(0,2,1) )
     &                    + ASMUF2 * ( BIGC(K,0,1) * EVOLA(1,1,1) + 
     &                                 BIGC(K,0,2) * EVOLA(1,2,1) )

          CEV(2) = CEV(2) + ASMUR2 * ( BIGC(K,1,1) * EVOLA(0,1,2) + 
     &                                 BIGC(K,1,2) * EVOLA(0,2,2) )
     &                    + ASMUF2 * ( BIGC(K,0,1) * EVOLA(1,1,2) + 
     &                                 BIGC(K,0,2) * EVOLA(1,2,2) )


          IF (P .GE. 2) THEN

            CEV(1) = CEV(1)+ASMUR2**2*( BIGC(K,2,1) * EVOLA(0,1,1) +
     &                                 BIGC(K,2,2) * EVOLA(0,2,1) )
     &        + ASMUR2 * ASMUF2 * ( BIGC(K,1,1) * EVOLA(1,1,1) + 
     &                              BIGC(K,1,2) * EVOLA(1,2,1) )
     &                 + ASMUF2**2 * ( BIGC(K,0,1) * EVOLA(2,1,1) + 
     &                                 BIGC(K,0,2) * EVOLA(2,2,1) )

            CEV(2) = CEV(2)+ASMUR2**2*( BIGC(K,2,1) * EVOLA(0,1,2) +
     &                                 BIGC(K,2,2) * EVOLA(0,2,2) )
     &        + ASMUR2 * ASMUF2 * ( BIGC(K,1,1) * EVOLA(1,1,2) + 
     &                              BIGC(K,1,2) * EVOLA(1,2,2) )
     &                 + ASMUF2**2 * ( BIGC(K,0,1) * EVOLA(2,1,2) + 
     &                                 BIGC(K,0,2) * EVOLA(2,2,2) )

          ENDIF

        ENDIF

        CGRID(QIND, K, 1) = CEV(1)
        CGRID(QIND, K, 2) = CEV(2)

       ELSE
*      --- DIS ---

        CEV(1) = ( BIGCF2(K,0,1) * EVOLA(0,1,1) + 
     &             BIGCF2(K,0,2) * EVOLA(0,2,1) ) * CZERO
        CEV(2) = ( BIGCF2(K,0,1) * EVOLA(0,1,2) + 
     &             BIGCF2(K,0,2) * EVOLA(0,2,2) ) * CZERO

        IF (P .GE. 1) THEN

          CEV(1) = CEV(1) + ASMUR2 * ( BIGCF2(K,1,1) * EVOLA(0,1,1) + 
     &                                 BIGCF2(K,1,2) * EVOLA(0,2,1) )
     &                    + ASMUF2 * ( BIGCF2(K,0,1) * EVOLA(1,1,1) + 
     &                                 BIGCF2(K,0,2) * EVOLA(1,2,1) )

          CEV(2) = CEV(2) + ASMUR2 * ( BIGCF2(K,1,1) * EVOLA(0,1,2) + 
     &                                 BIGCF2(K,1,2) * EVOLA(0,2,2) )
     &                    + ASMUF2 * ( BIGCF2(K,0,1) * EVOLA(1,1,2) + 
     &                                 BIGCF2(K,0,2) * EVOLA(1,2,2) )


          IF (P .GE. 2) THEN

          CEV(1) = CEV(1) + ASMUR2**2 * ( BIGCF2(K,2,1) * EVOLA(0,1,1) +
     &                                    BIGCF2(K,2,2) * EVOLA(0,2,1) )
     &        + ASMUR2 * ASMUF2 * ( BIGCF2(K,1,1) * EVOLA(1,1,1) + 
     &                              BIGCF2(K,1,2) * EVOLA(1,2,1) )
     &                  + ASMUF2**2 * ( BIGCF2(K,0,1) * EVOLA(2,1,1) + 
     &                                    BIGCF2(K,0,2) * EVOLA(2,2,1) )

          CEV(2) = CEV(2) + ASMUR2**2 * ( BIGCF2(K,2,1) * EVOLA(0,1,2) +
     &                                    BIGCF2(K,2,2) * EVOLA(0,2,2) )
     &        + ASMUR2 * ASMUF2 * ( BIGCF2(K,1,1) * EVOLA(1,1,2) + 
     &                              BIGCF2(K,1,2) * EVOLA(1,2,2) )
     &                  + ASMUF2**2 * ( BIGCF2(K,0,1) * EVOLA(2,1,2) + 
     &                                    BIGCF2(K,0,2) * EVOLA(2,2,2) )

          ENDIF

        ENDIF

        CGRIDDIS(QIND, K, 1) = CEV(1)
        CGRIDDIS(QIND, K, 2) = CEV(2)

       ENDIF


      END IF

 123  CONTINUE

      RETURN
      END
C     *****
