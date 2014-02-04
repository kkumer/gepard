C     ****h* gepard/hj.F
C  FILE DESCRIPTION
C    conformal moments of GPD ansaetze
C
C    $Id: parwav.F 71 2007-04-18 20:39:48Z kkumer $
C     *******


C     ****s* hj.F/HJ
C  NAME
C     HJ  --  conformal moment of input-scale singlet GPD  H_{J}
C  DESCRIPTION
C     returns H_{J} for various ansaetze, depending on the variable
C     ANSATZ
C  SYNOPSIS

      SUBROUTINE HJ(J, FCM)

      IMPLICIT NONE
      DOUBLE COMPLEX J, FCM(4) 

C  INPUTS
C           J -- conformal moment
C  OUTPUT
C         FCM -- input scale GPD H_{J}
C         FCM(1) -- singlet quark
C         FCM(2) -- gluon
C         FCM(3) -- non-siglet quark NS(+)
C         FCM(4) -- non-siglet quark NS(-)
C  PARENTS
C     INITGPD, GETMBGPD
C  CHILDREN
C     RATGAMMA, CLNGAMMA, CBETA, POCHHAMMER
C  SOURCE
C

      DOUBLE COMPLEX HSEA, HU, HD
      DOUBLE PRECISION NORMS
      DOUBLE COMPLEX CLNGAMMA, CBETA, POCHSEA, POCHG, POCHHAMMER
      DOUBLE COMPLEX JMA(2), CORR(2)
      DOUBLE COMPLEX QJ, BETADIP, BETAEXP
      INTEGER LHBETA(4), LHBETAF(4)
      DOUBLE PRECISION LHA(4), LHLAM(4), ALPT(2)
#include "header.F"
      DATA LHBETA / 3, 4, 6, 5 /
      DATA LHBETAF / 6, 24, 720, 120 /
      DATA LHA / 5.1072D0, 3.06432D0, 0.1939875D0, 1.7D0 /
      DATA LHLAM / -0.8D0, -0.8D0, 0.1D0, 0.1D0 /

*     Evolves with NS(-) so not used for most ansaetze
*     only exception is NSMHOU
      FCM(4) = ZERO
*      'Toy' singlet ansatz
      IF (ANSATZ .EQ. 'TOY') THEN
            FCM(1) = 454760.7514415856 * EXP(CLNGAMMA(0.5d0 + J)) /
     &            EXP(CLNGAMMA(10.6d0 + J))
            FCM(2) = 17.837861981813603 * EXP(CLNGAMMA(-0.1d0 + J)) /
     &            EXP(CLNGAMMA(4.7d0 + J))
            FCM(3) = ZERO
*      'Toy' non-singlet ansatz
      ELSE IF (ANSATZ .EQ. 'NSTOY') THEN
            FCM(1) = ZERO
            FCM(2) = ZERO
            FCM(3) = POCHHAMMER(DBLCMPLX(0.5d0, 0.0d0), 4) /
     &                  POCHHAMMER(DBLCMPLX(0.5d0, 0.0d0) + J, 4)
*  ## ansaetze used in PLB paper hep-ph/0605237
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
            FCM(1) = PAR(11) * CBETA(
     &      DBLCMPLX(1.0d0 - PAR(12) - PAR(13)*DEL2, 0.0d0) 
     &                                  + J, POCHSEA) / CBETA(
     &            DBLCMPLX(2.0d0 - PAR(12), 0.0d0), POCHSEA) /
     &            (1.0d0 - DEL2)**3
            FCM(2) = PAR(21) * CBETA(    
     &          DBLCMPLX(1.0d0-PAR(22)-PAR(23)*DEL2, 0.0d0) 
     &                                  + J, POCHG) / CBETA(
     &            DBLCMPLX(2.0d0 - PAR(22), 0.0d0), POCHG) /
     &            (1.0d0 - DEL2)**3
            FCM(3) = ZERO
*           Killing E contribution
            PAR(18) = 0.0d0
            PAR(28) = 0.0d0
*  ## ansatz used for testing (programs 'test' and 'accuracy')
      ELSE IF (ANSATZ .EQ. 'TEST') THEN
            FCM(1) = PAR(11) / (1 - DEL2/PAR(14)**2)**3 / POCHHAMMER(
     &            DBLCMPLX(1.0d0 - PAR(12) - PAR(13)*DEL2, 0.0d0)
     &             + J, 8) * POCHHAMMER(
     &            DBLCMPLX(2.0d0 - PAR(12), 0.0d0), 8)
            FCM(2) = PAR(21) / (1 - DEL2/PAR(24)**2)**3 / POCHHAMMER(
     &            DBLCMPLX(1.0d0 - PAR(22) - PAR(23)*DEL2, 0.0d0)
     &             + J, 6) * POCHHAMMER(
     &            DBLCMPLX(2.0d0 - PAR(22), 0.0d0), 6)
            FCM(3) = ZERO
*           Killing E contribution
            PAR(18) = 0.0d0
            PAR(28) = 0.0d0
      ELSE IF ((ANSATZ .EQ. 'FIT') .OR. (ANSATZ .EQ. 'EPH')
     &   .OR. (ANSATZ .EQ. 'EFL')) THEN
* Use first sum-rule constraint
          PAR(21) = 0.6d0 - PAR(11)
          FCM(1) = QJ(J, 9, PAR(11), PAR(12), PAR(13), PAR(19)) *
     &                BETADIP(J, PAR(14), PAR(15), PAR(16))
          FCM(2) = QJ(J, 7, PAR(21), PAR(22), PAR(23), PAR(29)) *
     &                BETADIP(J, PAR(24), PAR(25), PAR(26))
          FCM(3) = ZERO
      ELSE IF ((ANSATZ .EQ. 'FITEXP') .OR. (ANSATZ .EQ. 'EPHEXP')
     &   .OR. (ANSATZ .EQ. 'EFLEXP')) THEN
* Use first sum-rule constraint
          PAR(21) = 0.6d0 - PAR(11)
          FCM(1) = QJ(J, 9, PAR(11), PAR(12), PAR(13), PAR(19)) *
     &                BETAEXP(PAR(14))
          FCM(2) = QJ(J, 7, PAR(21), PAR(22), PAR(23), PAR(29)) *
     &                BETAEXP(PAR(24))
          FCM(3) = ZERO
      ELSE IF (ANSATZ .EQ. 'FITNG') THEN
* Don't use first sum-rule constraint (and not that pochs are different!)
          ALPT(1) = PAR(12) + PAR(13)*DEL2
          ALPT(2) = PAR(22) + PAR(23)*DEL2
          JMA(1) = J - ALPT(1)
          JMA(2) = J - ALPT(2)
          CORR(1) = EXP(   CLNGAMMA(DBLCMPLX(7.0d0 + PAR(12), 0.0d0)) -
     &           CLNGAMMA(DBLCMPLX(7.0d0 + ALPT(1), 0.0d0))   )
          CORR(2) = EXP(   CLNGAMMA(DBLCMPLX(7.0d0 + PAR(22), 0.0d0)) -
     &           CLNGAMMA(DBLCMPLX(7.0d0 + ALPT(2), 0.0d0))   )
        IF ( PID .EQ. 0 ) THEN 
*         ----  DIS  ----
          FCM(1) = PAR(11) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(12),0.0d0) 
     &              , 8) / POCHHAMMER(1.0d0 - PAR(12)
     &          + J , 8) * (1.0d0 - PAR(12) + J) / (1.0d0 + JMA(1))
          FCM(2) = PAR(21) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(22),0.0d0)
     &              , 6) / POCHHAMMER(1.0d0 - PAR(22)
     &          + J , 6) * (1.0d0 - PAR(22) + J) / (1.0d0 + JMA(2))
          FCM(3) = ZERO
        ELSE  
          FCM(1) = PAR(11) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(12),0.0d0)
     &              , 8) / POCHHAMMER(1.0d0 - PAR(12)
     &          + J , 8) * (1.0d0 - PAR(12) + J) / (1.0d0 + JMA(1)) /
     &           (1.0d0 - DEL2/(PAR(14)+PAR(15)*J))**PAR(16) *
     &        ( 1.0d0 - PAR(19)*PI*(XI/2.0d0)**(JMA(1)+1.0d0)*SQRT(PI)*
     &         2.0d0**(-2.0d0*ALPT(1))*EXP( CLNGAMMA(1.0d0+J) +
     &          CLNGAMMA(1.0d0+J+ALPT(1)) - CLNGAMMA(1.5d0+J) - 
     &          CLNGAMMA(1.0d0+JMA(1)) - 2.0d0*CLNGAMMA(
     &           DBLCMPLX(ALPT(1), 0.0d0)) ) / SIN(PI*ALPT(1)) *
     &          (1.0d0 - (SIN(PI*ALPT(1)) -
     &           SIN(J*PI))/SIN(PI*JMA(1)) ) * CORR(1)  )
          FCM(2) = PAR(21) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(22),0.0d0)
     &              , 6) / POCHHAMMER(1.0d0 - PAR(22)
     &          + J , 6) * (1.0d0 - PAR(22) + J) / (1.0d0 + JMA(2)) /
     &           (1.0d0 - DEL2/(PAR(24)+PAR(25)*J))**PAR(26) *
     &        ( 1.0d0 - PAR(29)*PI*(XI/2.0d0)**(JMA(2)+1.0d0)*SQRT(PI)*
     &         2.0d0**(-2.0d0*ALPT(2))*EXP( CLNGAMMA(1.0d0+J) +
     &          CLNGAMMA(1.0d0+J+ALPT(2)) - CLNGAMMA(1.5d0+J) - 
     &          CLNGAMMA(1.0d0+JMA(2)) - 2.0d0*CLNGAMMA(
     &           DBLCMPLX(ALPT(2), 0.0d0)) ) / SIN(PI*ALPT(2)) *
     &          (1.0d0 - (SIN(PI*ALPT(2)) -
     &           SIN(J*PI))/SIN(PI*JMA(2)) ) * CORR(2) )
          FCM(3) = ZERO
        ENDIF
*  ## ansatz used for "big paper"  hep-ph/0703179 
      ELSE IF (ANSATZ .EQ. 'FITBP') THEN
            HU = PAR(31) * POCHHAMMER(DBLCMPLX(1.0d0 - PAR(32), 0.0d0) 
     &              , 4) / POCHHAMMER(DBLCMPLX(1.0d0 - PAR(32), 0.0d0)
     &          + J , 4) * (DBLCMPLX(1.0d0 - PAR(32), 0.0d0) + J) /
     &           (DBLCMPLX(1.0d0 - PAR(32) - PAR(33)*DEL2, 0.0d0) + J) /
     &           (1.0d0 - DEL2/(PAR(34)+PAR(35)*J))**PAR(36)
            HD = PAR(41) * POCHHAMMER(DBLCMPLX(1.0d0 - PAR(42), 0.0d0) 
     &              , 4) / POCHHAMMER(DBLCMPLX(1.0d0 - PAR(42), 0.0d0)
     &          + J , 4) * (DBLCMPLX(1.0d0 - PAR(42), 0.0d0) + J) /
     &           (DBLCMPLX(1.0d0 - PAR(42) - PAR(43)*DEL2, 0.0d0) + J) /
     &           (1.0d0 - DEL2/(PAR(44)+PAR(45)*J))**PAR(46)
* Two options, uncomment only one!
*    1. Take sea normalization as free parameter
            NORMS = PAR(11)
*    2. Constrain sea normalization by momentum sum-rule (it has to
*           be declared as fixed in MINUIT.CMD then).
!            NORMS = 1.0d0 - PAR(21) - 2.0d0*(1.0d0-PAR(32))/(5.0d0 -
!     &               PAR(32)) - (1.0d0-PAR(42))/(5.0d0-PAR(42))
            HSEA = NORMS * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(12), 0.0d0) 
     &              , 8) / POCHHAMMER(DBLCMPLX(1.0d0 - PAR(12), 0.0d0)
     &          + J , 8) * (DBLCMPLX(1.0d0 - PAR(12), 0.0d0) + J) /
     &           (DBLCMPLX(1.0d0 - PAR(12) - PAR(13)*DEL2, 0.0d0) + J) /
     &           (1.0d0 - DEL2/(PAR(14)+PAR(15)*J))**PAR(16)
         FCM(1) = HSEA + HU + HD
         FCM(2) = PAR(21) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(22), 0.0d0) 
     &              , 6) / POCHHAMMER(DBLCMPLX(1.0d0 - PAR(22), 0.0d0)
     &          + J , 6) * (DBLCMPLX(1.0d0 - PAR(22), 0.0d0) + J) /
     &           (DBLCMPLX(1.0d0 - PAR(22) - PAR(23)*DEL2, 0.0d0) + J) /
     &           (1.0d0 - DEL2/(PAR(24)+PAR(25)*J))**PAR(26)
         FCM(3) = ZERO
*           Killing E contribution
            PAR(18) = 0.0d0
            PAR(28) = 0.0d0
      ELSE IF (ANSATZ .EQ. 'NSFIT') THEN
          FCM(1) = ZERO
          FCM(2) = ZERO
*     NS valence = down
            HD = PAR(41) * POCHHAMMER(DBLCMPLX(1.0d0 - PAR(42), 0.0d0) 
     &              , 4) / POCHHAMMER(DBLCMPLX(1.0d0 - PAR(42), 0.0d0)
     &          + J , 4) * (DBLCMPLX(1.0d0 - PAR(42), 0.0d0) + J) /
     &           (DBLCMPLX(1.0d0 - PAR(42) - PAR(43)*DEL2, 0.0d0) + J) /
     &           (1.0d0 - DEL2/(PAR(44)+PAR(45)*J))**PAR(46)
            NORMS = PAR(11)
            HSEA = NORMS * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(12), 0.0d0) 
     &              , 8) / POCHHAMMER(DBLCMPLX(1.0d0 - PAR(12), 0.0d0)
     &          + J , 8) * (DBLCMPLX(1.0d0 - PAR(12), 0.0d0) + J) /
     &           (DBLCMPLX(1.0d0 - PAR(12) - PAR(13)*DEL2, 0.0d0) + J) /
     &           (1.0d0 - DEL2/(PAR(14)+PAR(15)*J))**PAR(16)
          FCM(3) = HD - 0.5d0 / (2.0d0 + 0.5d0) * HSEA
*  ## ansatz for comparison with Les Houches PDF evolution benchmark
      ELSE IF (ANSATZ .EQ. 'HOUCHE') THEN
        FCM(1) =
     &      LHA(1)*LHBETAF(1) / POCHHAMMER(J-DBLCMPLX(LHLAM(1), 0.0D0),
     &        LHBETA(1)+1)
     &  + LHA(2)*LHBETAF(2) / POCHHAMMER(J-DBLCMPLX(LHLAM(2), 0.0D0),
     &        LHBETA(2)+1)
     &  + LHA(3)*LHBETAF(3) / POCHHAMMER(J-DBLCMPLX(LHLAM(3), 0.0D0),
     &        LHBETA(3)+1) * 4.8D0
     &  + LHA(3)*LHBETAF(3) / POCHHAMMER(J+1-DBLCMPLX(LHLAM(3), 0.0D0),
     &        LHBETA(3)+1) * (-2.4D0)
        FCM(2) = 
     &      LHA(4)*LHBETAF(4) / POCHHAMMER(J-DBLCMPLX(LHLAM(4), 0.0D0),
     &        LHBETA(4)+1) 
*           Killing E contribution
            PAR(18) = 0.0d0
            PAR(28) = 0.0d0
        FCM(3) = ZERO
      ELSE IF (ANSATZ .EQ. 'NSMHOU') THEN
*      (u_valence - d_valence), which is evolved by NS(-)
        FCM(1) = ZERO
        FCM(2) = ZERO
        FCM(3) = ZERO
        FCM(4) =
     &      LHA(1)*LHBETAF(1) / POCHHAMMER(J-DBLCMPLX(LHLAM(1), 0.0D0),
     &        LHBETA(1)+1)
     &  - LHA(2)*LHBETAF(2) / POCHHAMMER(J-DBLCMPLX(LHLAM(2), 0.0D0),
     &        LHBETA(2)+1)
      ELSE IF (ANSATZ .EQ. 'NSPHOU') THEN
*      (u_valence + 2ubar) - (d_valence + 2dbar), which is evolved by NS(+)
        FCM(1) =  ZERO
        FCM(2) =  ZERO
        FCM(3) =
     &      LHA(1)*LHBETAF(1) / POCHHAMMER(J-DBLCMPLX(LHLAM(1), 0.0D0),
     &        LHBETA(1)+1)
     &  - LHA(2)*LHBETAF(2) / POCHHAMMER(J-DBLCMPLX(LHLAM(2), 0.0D0),
     &        LHBETA(2)+1)
     &  - LHA(3)*LHBETAF(3) / POCHHAMMER(J+1-DBLCMPLX(LHLAM(3), 0.0D0),
     &        LHBETA(3)+1) * 2.D0
*           Killing E contribution
            PAR(18) = 0.0d0
            PAR(28) = 0.0d0
      END IF

      RETURN
      END
C     ****

C     ****f* hj.F/QJ
C  NAME
C     QJ  --  GPD ansatz from NPB841 (without residual t-dependence)
C  SYNOPSIS

      DOUBLE COMPLEX FUNCTION QJ(J, POCH, NORM, AL0, ALP, SKEW)

      IMPLICIT NONE
      INTEGER POCH
      DOUBLE COMPLEX J
      DOUBLE PRECISION NORM, AL0, ALP, SKEW

C  INPUTS
C          J -- conformal moment
C       POCH -- Pochammer index  
C     (rest) -- model parameters
C  RETURN VALUE
C         QJ -- Q_J  (conformal moment of GPD with Reggeized t-dep,
C                     but WITHOUT residual \beta(t) from Eq. (19) )
C  PARENTS
C     HJ
C  CHILDREN
C     CLNGAMMA, POCHHAMMER
C  SOURCE
C

      DOUBLE COMPLEX CLNGAMMA, POCHHAMMER
      DOUBLE COMPLEX JMA, CORR
      DOUBLE PRECISION ALPT
#include "header.F"

*     NG = 0.6d0 - N
      ALPT = AL0 + ALP*DEL2
      JMA = J - ALPT
      IF ( PID .EQ. 0 ) THEN 
*       ---  DIS ----
        QJ = NORM * POCHHAMMER(DBLCMPLX(2.0d0 - AL0,0.0d0) 
     &            , POCH) / POCHHAMMER(1.0d0 - AL0
     &        + J , POCH) * (1.0d0 - AL0 + J) / (1.0d0 + JMA)
      ELSE
*       correction factor from Eq. (57), same for quarks and G
        CORR = EXP(   CLNGAMMA(DBLCMPLX(7.0d0 + AL0, 0.0d0)) -
     &         CLNGAMMA(DBLCMPLX(7.0d0 + ALPT, 0.0d0))   )
        IF ((ABS(SKEW) .LE. 1e-6) .OR. (XI .LE. 1e-12)) THEN
*          don't need resummed subleading waves
          QJ = NORM * POCHHAMMER(DBLCMPLX(2.0d0 - AL0, 0.0d0)
     &            , POCH) / POCHHAMMER(1.0d0 - AL0
     &        + J , POCH) * (1.0d0 - AL0 + J) / (1.0d0 + JMA)
        ELSE
          QJ = NORM * POCHHAMMER(DBLCMPLX(2.0d0 - AL0, 0.0d0)
     &            , POCH) / POCHHAMMER(1.0d0 - AL0
     &        + J , POCH) * (1.0d0 - AL0 + J) / (1.0d0 + JMA) *
     &      ( 1.0d0 - SKEW*PI*(XI/2.0d0)**(JMA+1.0d0)*SQRT(PI)*
     &       2.0d0**(-2.0d0*ALPT)*EXP( CLNGAMMA(1.0d0+J) +
     &        CLNGAMMA(1.0d0+J+ALPT) - CLNGAMMA(1.5d0+J) - 
     &        CLNGAMMA(1.0d0+JMA) - 2.0d0*CLNGAMMA(
     &         DBLCMPLX(ALPT, 0.0d0)) ) / SIN(PI*ALPT) *
     &        (1.0d0 - (SIN(PI*ALPT) -
     &         SIN(J*PI))/SIN(PI*JMA) ) * CORR  )
        END IF
      END IF
      RETURN
      END
C     ***

C     ****f* hj.F/BETADIP
C  NAME
C     BETADIP  -- residual t-dependence function \beta from Eq.(19)
C                 (dipole version)
C  SYNOPSIS

      DOUBLE COMPLEX FUNCTION BETADIP(J, M02, DELM2, PP)

      IMPLICIT NONE
      DOUBLE COMPLEX J
      DOUBLE PRECISION M02, DELM2, PP

C  INPUTS
C          J -- conformal moment
C     (rest) -- model parameters
C  RETURN VALUE
C     BETADIP -- residual t-dep   1/(1-t/M2)^p ;  M2=M02+J*DELM2
C  PARENTS
C     HJ
C  CHILDREN
C     CLNGAMMA, POCHHAMMER
C  SOURCE
C

#include "header.F"

      BETADIP =  1.0d0 / (1.0d0 - DEL2/(M02+DELM2*J))**PP
      RETURN
      END
C     ***

C     ****f* hj.F/BETAEXP
C  NAME
C     BETAEXP  -- residual t-dependence function \beta from Eq.(19)
C                 (exponential version)
C  SYNOPSIS

      DOUBLE COMPLEX FUNCTION BETAEXP(M02)

      IMPLICIT NONE
      DOUBLE PRECISION M02

C  INPUTS
C          J -- conformal moment
C     (rest) -- model parameters
C  RETURN VALUE
C     BETAEXP -- residual t-dep   exp(t/2/M2) ; 
C  PARENTS
C     HJ
C  CHILDREN
C     CLNGAMMA, POCHHAMMER
C  SOURCE
C

#include "header.F"

      BETAEXP =  EXP(DEL2 / (2.d0 * M02))
      RETURN
      END
C     ***
