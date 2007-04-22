C     ****h* gepard/evolc.f
C  FILE DESCRIPTION
C    calculates "evolved" C's
C
C    $Id: parwav.f 71 2007-04-18 20:39:48Z kkumer $
C     *******

C     ****s* parwav.f/PARWAVF
C  NAME
C     PARWAVF  --  J-th conformal partial wave for DVCS and DIS
C  DESCRIPTION
C    calculates J-th conformal partial wave FPW i.e. 
C    |latex $\vec{C}_j \cdot \vec{H}_j$
C    |html  &nbsp; &nbsp;&nbsp;&nbsp;&nbsp; <b><i>C</i></b><sub>j</sub> <b>. <i>H</i></b><sub>j</sub>
C    (together with evolution), and in DVCS case together with
C    prefactor Gamma(5/2+J) / Gamma(3+J)
C  SYNOPSIS
C     SUBROUTINE PARWAVF (K, FPW, PROCESS)
C
C     INTEGER K
C     DOUBLE COMPLEX FPW
C     CHARACTER PROCESS*4
C  INPUTS
C           K -- Mellin-Barnes integration point index
C  OUTPUT
C         FPW -- partial wave
C  PARENTS
C      CFFF, F2F
C  CHILDREN
C      AS2PF, EVOLF, HJ
C  SOURCE
C

      SUBROUTINE EVOLC (SEC, QIND)

      IMPLICIT NONE
      INTEGER K, SEC, QIND
      INTEGER ORD, L, K1
      DOUBLE PRECISION R, ASQ02, ASMUR2, ASMUF2
      DOUBLE COMPLEX CEVNS, CNDNS, CEV(2)
      DOUBLE COMPLEX EVOLNSA(0:2)
      DOUBLE COMPLEX J, EVOLA(0:2,2,2), FCM(2)
      INCLUDE 'header.f'

      PAR(1) = 4.0d0
      PAR(2) = 0.0488d0
      PAR(3) = 2.5d0

* ASMUF2 = alpha_s/(2 pi) at factorization scale
      CALL AS2PF (ASMUF2, Q2/RF2, PAR(2), PAR(3))
* ASMUR2 = alpha_s/(2 pi) at renormalization scale
      CALL AS2PF (ASMUR2, Q2/RR2, PAR(2), PAR(3))
* ASQ02 = alpha_s/(2 pi) at input scale
      CALL AS2PF (ASQ02, PAR(1), PAR(2), PAR(3))

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
     &                  + ASMUF2 * BIGCNS(K, 0) * EVOLNSA(1) 

          IF (P .GE. 2) THEN

            CEVNS = CEVNS + ASMUR2**2 * BIGCNS(K,2) * EVOLNSA(0) +
     &          ASMUR2 * ASMUF2 * BIGCNS(K,1) * EVOLNSA(1) + 
     &                      ASMUF2**2 *  BIGCNS(K,0) * EVOLNSA(2)
          END IF

        END IF

        CGRIDNS(QIND, K) = CEVNS

      ELSE
*     singlet case

       CALL EVOLF (SEC, K, R, EVOLA)

       IF ( PROCESS(:3) .EQ. 'DVC' ) THEN
*      --- DVCS ---

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

        CGRID(SEC, QIND, K, 1) = CEV(1)
        CGRID(SEC, QIND, K, 2) = CEV(2)

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
