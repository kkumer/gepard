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
C     PROCESS -- 'DVCS' or 'DIS'
C  OUTPUT
C         FPW -- partial wave
C  PARENTS
C      CFFF, F2F
C  CHILDREN
C      AS2PF, EVOLF, HJ
C  SOURCE
C

      SUBROUTINE PARWAVF (K, FPW, PROCESS)

      IMPLICIT NONE
      INTEGER K
      DOUBLE COMPLEX FPW
      CHARACTER PROCESS*4
      INTEGER ORD, L, K1
      DOUBLE PRECISION R, ASQ2, ASQ02
      DOUBLE COMPLEX CDVCS(0:2,2), CNDNS
      DOUBLE COMPLEX EVOLNSA(0:2)
      DOUBLE COMPLEX J, BIGCL(0:2,2), EVOLA(0:2,2,2), FCM(2)
      INCLUDE 'header.f'

      J = N(K) - 1
      CALL AS2PF (ASQ2, Q2, PAR(2), PAR(3))
      CALL AS2PF (ASQ02, PAR(1), PAR(2), PAR(3))
      R = ASQ2/ASQ02
      CALL HJ(J, FCM)

*  "Big C" Wilson coefficients depending on process type:
*  (They also depend on scheme, but this is taken care of by INIT.)

      DO 10 ORD = 0, P
      DO 10 L = 1, 2
        CDVCS(ORD, L) = (0.0d0, 0.0d0)
        IF ( PROCESS .EQ. 'DVCS' ) THEN
          BIGCL(ORD, L) = BIGC(K, ORD, L)
        ELSE IF ( PROCESS(:3) .EQ. 'DIS' ) THEN
          BIGCL(ORD, L) = BIGCF2(K, ORD, L)
        ELSE 
          CALL ERROR ('GeParD', 'PARWAVF',
     &    'Process ' // PROCESS // ' is unknown! (Need DVCS or DIS)   ',
     &    6, 3)
        END IF
 10   CONTINUE


*  Evolution operators and partial wave depending on whether we do
*  non-singlet or singlet case. This is specified by naming non-singlet
*  ansaetze starting with 'NS...'.

      IF ( ANSATZ(:2) .EQ. 'NS' ) THEN
*   non-singlet case

*     FIXME: In the non-singlet case NNLO Wilson coeffs. are not implemented
*     Also, ND evolution is not treated correctly (ND NLO evol. operator assumes
*     it will multiply C^(0)!)
      IF ( P .GE. 2 ) THEN
          CALL ERROR ('GeParD', 'PARWAVF',
     &    'P>1, which is not implemented in the non-singlet case! ',
     &    7, 2)
      END IF

      CALL EVOLNSF (K, R, EVOLNSA)

      DO 15 ORD = 0, P
      DO 15 K1 = 0, ORD
        CDVCS(ORD, 1) = CDVCS(ORD, 1) + BIGCL(ORD-K1, 1) * EVOLNSA(K1)
 15   CONTINUE

*   Put astrong^0 term to zero if investigating NLO effects
      CDVCS(0, 1) = CZERO * CDVCS(0, 1)

      FPW = (0.0d0, 0.0d0)
      DO 18 ORD=0, P
 18   FPW = FPW + ASQ2**ORD * (CDVCS(ORD,1)*FCM(1))

      ELSE
*   singlet case

      CALL EVOLF (K, R, EVOLA)

      DO 20 ORD = 0, P
      DO 20 K1 = 0, ORD
        CDVCS(ORD, 1) = CDVCS(ORD, 1) + BIGCL(ORD-K1, 1) * EVOLA(K1,1,1)
     &      + BIGCL(ORD-K1,2) * EVOLA(K1,2,1)
        CDVCS(ORD, 2) = CDVCS(ORD, 2) + BIGCL(ORD-K1, 1) * EVOLA(K1,1,2)
     &      + BIGCL(ORD-K1,2) * EVOLA(K1,2,2)
 20   CONTINUE

      FPW = (0.0d0, 0.0d0)
      DO 30 ORD=0, P
 30   FPW = FPW + ASQ2**ORD * (CDVCS(ORD,1)*FCM(1)+CDVCS(ORD,2)*FCM(2))

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
     &            COMPLEX(2.0d0 - PAR(12), 0.0d0), POCHSEA)
            FCM(2) = PAR(21) * CBETA(COMPLEX(1.0d0-PAR(22)-PAR(23)*DEL2,
     &        0.0d0) + J, POCHG) / CBETA(COMPLEX(2.0d0 - PAR(22),
     &        0.0d0), POCHG)
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
            NORMS = PAR(11)
            HSEA = NORMS * POCHHAMMER(COMPLEX(2.0d0 - PAR(12), 0.0d0) 
     &              , 8) / POCHHAMMER(COMPLEX(1.0d0 - PAR(12), 0.0d0)
     &          + J , 8) * (COMPLEX(1.0d0 - PAR(12), 0.0d0) + J) /
     &           (COMPLEX(1.0d0 - PAR(12) - PAR(13)*DEL2, 0.0d0) + J) /
     &           (1.0d0 - DEL2/(PAR(14)+PAR(15)*J))**PAR(16)
*    Ugly kludge: NS valence = down
          FCM(1) = HD - (0.5d0 / 2.0d0 + 0.5d0) * HSEA
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
