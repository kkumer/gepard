C     ****h* gepard/parwav.f
C  FILE DESCRIPTION
C    calculates J-th conformal partial wave for DVCS
C
C    $Id$
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

      SUBROUTINE PARWAVF (SEC, K, FPW)

      IMPLICIT NONE
      INTEGER SEC, K
      DOUBLE COMPLEX FPW
      INTEGER ORD, L, K1
      DOUBLE PRECISION R, ASQ02, ASMUR2, ASMUF2
      DOUBLE COMPLEX CEVNS, CNDNS, CEV(2)
      DOUBLE COMPLEX EVOLNSA(0:2)
      DOUBLE COMPLEX J, EVOLA(0:2,2,2), FCM(2)
      INCLUDE 'header.f'

*   J is needed only for GPD, and GPD is on original
*    MB contour even for shifted (2nd) partial wave
      J = N(1,K) - 1

* ASMUF2 = alpha_s/(2 pi) at factorization scale
      CALL AS2PF (ASMUF2, Q2/RF2, PAR(2), PAR(3))
* ASMUR2 = alpha_s/(2 pi) at renormalization scale
      CALL AS2PF (ASMUR2, Q2/RR2, PAR(2), PAR(3))
* ASQ02 = alpha_s/(2 pi) at input scale
      CALL AS2PF (ASQ02, PAR(1), PAR(2), PAR(3))

      R = ASMUF2/ASQ02

      CALL HJ(J, FCM)


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

        FPW = CEVNS * FCM(1)

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

       ENDIF
       IF (SEC .EQ. 2) THEN
*        second partial wave is needed:
         FPW = PAR(50) * CEV(1) * FCM(1) + PAR(51) * CEV(2) * FCM(2)
       ELSE
         FPW = CEV(1) * FCM(1) + CEV(2) * FCM(2)
       ENDIF

      END IF

      RETURN
      END
C     *****


C     ****s* parwav.f/HJ
C  NAME
C     HJ  --  conformal moment of input-scale singlet GPD  H_{J}
C  DESCRIPTION
C     returns H_{J} for various ansaetze
C  SYNOPSIS
C     SUBROUTINE HJ(J, FCM)
C
C     DOUBLE COMPLEX J, FCM(2) 
C  INPUTS
C           J -- conformal moment
C
C    (Following parameters are inputs only for ANSATZ='FIT')
C
C     ALPHA0SEA, ALPHA0G -- intercepts of Regge trajectories for
C                           "sea" quark and gluon GPDs/PDFs
C               NSEA, NG -- normalization factors
C             MSEA2, MG2 -- positions of DEL2 triple poles
C  OUTPUT
C         FCM -- input scale singlet GPD H_{J}
C  PARENTS
C     PARWAVF
C  CHILDREN
C     CLNGAMMA, CBETA, POCHHAMMER
C  SOURCE
C

      SUBROUTINE HJ(J, FCM)

      IMPLICIT NONE
      DOUBLE COMPLEX J, FCM(2) 
      DOUBLE COMPLEX HSEA, HU, HD
      DOUBLE PRECISION NORMS
      DOUBLE COMPLEX CLNGAMMA, CBETA, POCHSEA, POCHG, POCHHAMMER
      DOUBLE COMPLEX NUM, DENN
      INTEGER LHBETA(4), LHBETAF(4)
      DOUBLE PRECISION LHA(4), LHLAM(4)
      INCLUDE 'header.f'
      DATA LHBETA / 3, 4, 6, 5 /
      DATA LHBETAF / 6, 24, 720, 120 /
      DATA LHA / 5.1072D0, 3.06432D0, 0.1939875D0, 1.7D0 /
      DATA LHLAM / -0.8D0, -0.8D0, 0.1D0, 0.1D0 /

*      'Toy' singlet ansatz
      IF (ANSATZ .EQ. 'TOY') THEN
            FCM(1) = 454760.7514415856 * EXP(CLNGAMMA(0.5d0 + J)) /
     &            EXP(CLNGAMMA(10.6d0 + J))
            FCM(2) = 17.837861981813603 * EXP(CLNGAMMA(-0.1d0 + J)) /
     &            EXP(CLNGAMMA(4.7d0 + J))
*      'Toy' non-singlet ansatz
      ELSE IF (ANSATZ .EQ. 'NSTOY') THEN
            FCM(1) = POCHHAMMER(COMPLEX(0.5d0, 0.0d0), 4) /
     &                  POCHHAMMER(COMPLEX(0.5d0, 0.0d0) + J, 4)
*          'Dummy' gluonic, not used:
            FCM(2) = (0.0d0, 0.0d0)
      ELSE IF ((ANSATZ .EQ. 'HARD') .OR. (ANSATZ .EQ. 'SOFT')) THEN
            PAR(12) = 1.1d0
            PAR(13) = 0.25d0
            PAR(23) = 0.25d0
            POCHSEA = (8.0d0, 0.0d0)
            POCHG = (6.0d0, 0.0d0)
            IF (ANSATZ .EQ. 'HARD') THEN
                  PAR(21) = 0.4d0
                  PAR(22) = PAR(12) + 0.1d0
            ELSE IF (ANSATZ .EQ. 'SOFT') THEN
                  PAR(21) = 0.3d0
                  PAR(22) = PAR(12)
            END IF
            PAR(11) = (2.0d0/3.0d0) - PAR(21)
            FCM(1) = PAR(11) * CBETA(COMPLEX(1.0d0 - PAR(12) - 
     &            PAR(13)*DEL2, 0.0d0) + J, POCHSEA) / CBETA(
     &            COMPLEX(2.0d0 - PAR(12), 0.0d0), POCHSEA) /
     &            (1.0d0 - DEL2)**3
            FCM(2) = PAR(21) * CBETA(COMPLEX(1.0d0-PAR(22)-PAR(23)*DEL2,
     &        0.0d0) + J, POCHG) / CBETA(COMPLEX(2.0d0 - PAR(22),
     &        0.0d0), POCHG) /
     &            (1.0d0 - DEL2)**3
      ELSE IF (ANSATZ .EQ. 'FITOLD') THEN
            FCM(1) = PAR(11) / (1 - DEL2/PAR(14)**2)**3 / POCHHAMMER(
     &            COMPLEX(1.0d0 - PAR(12) - 
     &            PAR(13)*DEL2, 0.0d0) + J, 8) * POCHHAMMER(
     &            COMPLEX(2.0d0 - PAR(12), 0.0d0), 8)
            FCM(2) = PAR(21) / (1 - DEL2/PAR(24)**2)**3 / POCHHAMMER(
     &            COMPLEX(1.0d0 - PAR(22) - PAR(23)*DEL2,
     &        0.0d0) + J, 6) * POCHHAMMER(COMPLEX(2.0d0 - PAR(22),
     &        0.0d0), 6)
      ELSE IF (ANSATZ .EQ. 'FIT') THEN
            HU = PAR(31) * POCHHAMMER(COMPLEX(1.0d0 - PAR(32), 0.0d0) 
     &              , 4) / POCHHAMMER(COMPLEX(1.0d0 - PAR(32), 0.0d0)
     &          + J , 4) * (COMPLEX(1.0d0 - PAR(32), 0.0d0) + J) /
     &           (COMPLEX(1.0d0 - PAR(32) - PAR(33)*DEL2, 0.0d0) + J) /
     &           (1.0d0 - DEL2/(PAR(34)+PAR(35)*J))**PAR(36)
            HD = PAR(41) * POCHHAMMER(COMPLEX(1.0d0 - PAR(42), 0.0d0) 
     &              , 4) / POCHHAMMER(COMPLEX(1.0d0 - PAR(42), 0.0d0)
     &          + J , 4) * (COMPLEX(1.0d0 - PAR(42), 0.0d0) + J) /
     &           (COMPLEX(1.0d0 - PAR(42) - PAR(43)*DEL2, 0.0d0) + J) /
     &           (1.0d0 - DEL2/(PAR(44)+PAR(45)*J))**PAR(46)
* Two options, uncomment only one!
*    1. Take sea normalization as free parameter
            NORMS = PAR(11)
*    2. Constrain sea normalization by momentum sum-rule (it has to
*           be declared as fixed in MINUIT.CMD then).
!            NORMS = 1.0d0 - PAR(21) - 2.0d0*(1.0d0-PAR(32))/(5.0d0 -
!     &               PAR(32)) - (1.0d0-PAR(42))/(5.0d0-PAR(42))
            HSEA = NORMS * POCHHAMMER(COMPLEX(2.0d0 - PAR(12), 0.0d0) 
     &              , 8) / POCHHAMMER(COMPLEX(1.0d0 - PAR(12), 0.0d0)
     &          + J , 8) * (COMPLEX(1.0d0 - PAR(12), 0.0d0) + J) /
     &           (COMPLEX(1.0d0 - PAR(12) - PAR(13)*DEL2, 0.0d0) + J) /
     &           (1.0d0 - DEL2/(PAR(14)+PAR(15)*J))**PAR(16)
          FCM(1) = HSEA + HU + HD
          FCM(2) = PAR(21) * POCHHAMMER(COMPLEX(2.0d0 - PAR(22), 0.0d0) 
     &              , 6) / POCHHAMMER(COMPLEX(1.0d0 - PAR(22), 0.0d0)
     &          + J , 6) * (COMPLEX(1.0d0 - PAR(22), 0.0d0) + J) /
     &           (COMPLEX(1.0d0 - PAR(22) - PAR(23)*DEL2, 0.0d0) + J) /
     &           (1.0d0 - DEL2/(PAR(24)+PAR(25)*J))**PAR(26)
      ELSE IF (ANSATZ .EQ. 'NSFIT') THEN
*    Ugly kludge: NS valence = down
            HD = PAR(41) * POCHHAMMER(COMPLEX(1.0d0 - PAR(42), 0.0d0) 
     &              , 4) / POCHHAMMER(COMPLEX(1.0d0 - PAR(42), 0.0d0)
     &          + J , 4) * (COMPLEX(1.0d0 - PAR(42), 0.0d0) + J) /
     &           (COMPLEX(1.0d0 - PAR(42) - PAR(43)*DEL2, 0.0d0) + J) /
     &           (1.0d0 - DEL2/(PAR(44)+PAR(45)*J))**PAR(46)
            NORMS = PAR(11)
            HSEA = NORMS * POCHHAMMER(COMPLEX(2.0d0 - PAR(12), 0.0d0) 
     &              , 8) / POCHHAMMER(COMPLEX(1.0d0 - PAR(12), 0.0d0)
     &          + J , 8) * (COMPLEX(1.0d0 - PAR(12), 0.0d0) + J) /
     &           (COMPLEX(1.0d0 - PAR(12) - PAR(13)*DEL2, 0.0d0) + J) /
     &           (1.0d0 - DEL2/(PAR(14)+PAR(15)*J))**PAR(16)
          FCM(1) = HD - 0.5d0 / (2.0d0 + 0.5d0) * HSEA
*    'Dummy' gluonic, not used:
          FCM(2) = (0.0d0, 0.0d0)
      ELSE IF (ANSATZ .EQ. 'HOUCHE') THEN
        FCM(1) =
     &      LHA(1)*LHBETAF(1) / POCHHAMMER(J-COMPLEX(LHLAM(1), 0.0D0),
     &        LHBETA(1)+1)
     &    + LHA(2)*LHBETAF(2) / POCHHAMMER(J-COMPLEX(LHLAM(2), 0.0D0),
     &        LHBETA(2)+1)
     &    + LHA(3)*LHBETAF(3) / POCHHAMMER(J-COMPLEX(LHLAM(3), 0.0D0),
     &        LHBETA(3)+1) * 4.8D0
     &    + LHA(3)*LHBETAF(3) / POCHHAMMER(J+1-COMPLEX(LHLAM(3), 0.0D0),
     &        LHBETA(3)+1) * (-2.4D0)
        FCM(2) = 
     &      LHA(4)*LHBETAF(4) / POCHHAMMER(J-COMPLEX(LHLAM(4), 0.0D0),
     &        LHBETA(4)+1) 
      END IF

      RETURN
      END
C     ****
