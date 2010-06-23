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
      DOUBLE COMPLEX J, FCM(2) 

C  INPUTS
C           J -- conformal moment
C  OUTPUT
C         FCM -- input scale singlet GPD H_{J}
C  PARENTS
C     INITGPD, GETMBGPD
C  CHILDREN
C     RATGAMMA, CLNGAMMA, CBETA, POCHHAMMER
C  SOURCE
C

      DOUBLE COMPLEX HSEA, HU, HD
      DOUBLE PRECISION NORMS
      DOUBLE COMPLEX CLNGAMMA, CBETA, POCHSEA, POCHG, POCHHAMMER
      DOUBLE COMPLEX RATGAMMA, JMA(2)
      DOUBLE COMPLEX MULTIGAMMA, CORR(2)
      INTEGER LHBETA(4), LHBETAF(4)
      DOUBLE PRECISION LHA(4), LHLAM(4), ALPT(2)
#include "header.F"
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
            FCM(1) = POCHHAMMER(DBLCMPLX(0.5d0, 0.0d0), 4) /
     &                  POCHHAMMER(DBLCMPLX(0.5d0, 0.0d0) + J, 4)
*          'Dummy' gluonic, not used:
            FCM(2) = (0.0d0, 0.0d0)
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
*           Killing E contribution
            PAR(18) = 0.0d0
            PAR(28) = 0.0d0
      ELSE IF (ANSATZ .EQ. 'FIT') THEN
* Use first sum-rule constraint
          PAR(21) = 0.6d0 - PAR(11)
          ALPT(1) = PAR(12) + PAR(13)*DEL2
          ALPT(2) = PAR(22) + PAR(23)*DEL2
          JMA(1) = J - ALPT(1)
          JMA(2) = J - ALPT(2)
          CORR(1) = EXP(   CLNGAMMA(DBLCMPLX(7.0d0 + PAR(12), 0.0d0)) -
     &           CLNGAMMA(DBLCMPLX(7.0d0 + ALPT(1), 0.0d0))   )
          CORR(2) = EXP(   CLNGAMMA(DBLCMPLX(7.0d0 + PAR(22), 0.0d0)) -
     &           CLNGAMMA(DBLCMPLX(7.0d0 + ALPT(2), 0.0d0))   )
        IF ( PROCESS .EQ. 'DIS' ) THEN 
          FCM(1) = PAR(11) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(12),0.0d0) 
     &              , 9) / POCHHAMMER(1.0d0 - PAR(12)
     &          + J , 9) * (1.0d0 - PAR(12) + J) / (1.0d0 + JMA(1))
          FCM(2) = PAR(21) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(22),0.0d0)
     &              , 7) / POCHHAMMER(1.0d0 - PAR(22)
     &          + J , 7) * (1.0d0 - PAR(22) + J) / (1.0d0 + JMA(2))
        ELSE  
          FCM(1) = PAR(11) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(12),0.0d0)
     &              , 9) / POCHHAMMER(1.0d0 - PAR(12)
     &          + J , 9) * (1.0d0 - PAR(12) + J) / (1.0d0 + JMA(1)) /
     &           (1.0d0 - DEL2/(PAR(14)+PAR(15)*J))**PAR(16) *
     &        ( 1.0d0 - PAR(19)*PI*(XI/2.0d0)**(JMA(1)+1.0d0)*SQRT(PI)*
     &         2.0d0**(-2.0d0*ALPT(1))*EXP( CLNGAMMA(1.0d0+J) +
     &          CLNGAMMA(1.0d0+J+ALPT(1)) - CLNGAMMA(1.5d0+J) - 
     &          CLNGAMMA(1.0d0+JMA(1)) - 2.0d0*CLNGAMMA(
     &           DBLCMPLX(ALPT(1), 0.0d0)) ) / SIN(PI*ALPT(1)) *
     &          (1.0d0 - (SIN(PI*ALPT(1)) -
     &           SIN(J*PI))/SIN(PI*JMA(1)) ) * CORR(1)  )
          FCM(2) = PAR(21) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(22),0.0d0)
     &              , 7) / POCHHAMMER(1.0d0 - PAR(22)
     &          + J , 7) * (1.0d0 - PAR(22) + J) / (1.0d0 + JMA(2)) /
     &           (1.0d0 - DEL2/(PAR(24)+PAR(25)*J))**PAR(26) *
     &        ( 1.0d0 - PAR(29)*PI*(XI/2.0d0)**(JMA(2)+1.0d0)*SQRT(PI)*
     &         2.0d0**(-2.0d0*ALPT(2))*EXP( CLNGAMMA(1.0d0+J) +
     &          CLNGAMMA(1.0d0+J+ALPT(2)) - CLNGAMMA(1.5d0+J) - 
     &          CLNGAMMA(1.0d0+JMA(2)) - 2.0d0*CLNGAMMA(
     &           DBLCMPLX(ALPT(2), 0.0d0)) ) / SIN(PI*ALPT(2)) *
     &          (1.0d0 - (SIN(PI*ALPT(2)) -
     &           SIN(J*PI))/SIN(PI*JMA(2)) ) * CORR(2) )
        ENDIF
      ELSE IF (ANSATZ .EQ. 'FIT86') THEN
* Use first sum-rule constraint
          PAR(21) = 0.6d0 - PAR(11)
          ALPT(1) = PAR(12) + PAR(13)*DEL2
          ALPT(2) = PAR(22) + PAR(23)*DEL2
          JMA(1) = J - ALPT(1)
          JMA(2) = J - ALPT(2)
          CORR(1) = EXP(   CLNGAMMA(DBLCMPLX(7.0d0 + PAR(12), 0.0d0)) -
     &           CLNGAMMA(DBLCMPLX(7.0d0 + ALPT(1), 0.0d0))   )
          CORR(2) = EXP(   CLNGAMMA(DBLCMPLX(7.0d0 + PAR(22), 0.0d0)) -
     &           CLNGAMMA(DBLCMPLX(7.0d0 + ALPT(2), 0.0d0))   )
        IF ( PROCESS .EQ. 'DIS' ) THEN 
          FCM(1) = PAR(11) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(12),0.0d0) 
     &              , 8) / POCHHAMMER(1.0d0 - PAR(12)
     &          + J , 8) * (1.0d0 - PAR(12) + J) / (1.0d0 + JMA(1))
          FCM(2) = PAR(21) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(22),0.0d0)
     &              , 6) / POCHHAMMER(1.0d0 - PAR(22)
     &          + J , 6) * (1.0d0 - PAR(22) + J) / (1.0d0 + JMA(2))
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
        ENDIF
      ELSE IF (ANSATZ .EQ. 'FITMUL') THEN
* Use first sum-rule constraint
          PAR(21) = 0.6d0 - PAR(11)
          ALPT(1) = PAR(12) + PAR(13)*DEL2
          ALPT(2) = PAR(22) + PAR(23)*DEL2
          JMA(1) = J - ALPT(1)
          JMA(2) = J - ALPT(2)
        IF ( PROCESS .EQ. 'DIS' ) THEN 
          FCM(1) = PAR(11) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(12),0.0d0) 
     &              , 9) / POCHHAMMER(1.0d0 - PAR(12)
     &          + J , 9) * (1.0d0 - PAR(12) + J) / (1.0d0 + JMA(1))
          FCM(2) = PAR(21) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(22),0.0d0)
     &              , 7) / POCHHAMMER(1.0d0 - PAR(22)
     &          + J , 7) * (1.0d0 - PAR(22) + J) / (1.0d0 + JMA(2))
        ELSE  
          FCM(1) = PAR(11) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(12),0.0d0)
     &              , 9) / POCHHAMMER(1.0d0 - PAR(12)
     &          + J , 9) * (1.0d0 - PAR(12) + J) / (1.0d0 + JMA(1)) /
     &           (1.0d0 - DEL2/(PAR(14)+PAR(15)*J))**PAR(16) *
     &        ( 1.0d0 - PAR(19)*PI*(XI/2.0d0)**(JMA(1)+1.0d0)*SQRT(PI)*
     &         2.0d0**(-2.0d0*ALPT(1))*
     &    MULTIGAMMA(
     &       1.0d0+J, 1.0d0+J+ALPT(1), DBLCMPLX(7.0d0 + PAR(12), 0.0d0),
     &        1.5d0+J, 1.0d0+JMA(1), DBLCMPLX(7.0d0 + ALPT(1), 0.0d0),
     &        DBLCMPLX(ALPT(1), 0.0d0) ) / SIN(PI*ALPT(1)) *
     &          (1.0d0 - (SIN(PI*ALPT(1)) -
     &           SIN(J*PI))/SIN(PI*JMA(1)) ) )
          FCM(2) = PAR(21) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(22),0.0d0)
     &              , 7) / POCHHAMMER(1.0d0 - PAR(22)
     &          + J , 7) * (1.0d0 - PAR(22) + J) / (1.0d0 + JMA(2)) /
     &           (1.0d0 - DEL2/(PAR(24)+PAR(25)*J))**PAR(26) *
     &        ( 1.0d0 - PAR(29)*PI*(XI/2.0d0)**(JMA(2)+1.0d0)*SQRT(PI)*
     &         2.0d0**(-2.0d0*ALPT(2))*
     &    MULTIGAMMA(
     &      1.0d0+J, 1.0d0+J+ALPT(2), DBLCMPLX(7.0d0 + PAR(22), 0.0d0),
     &        1.5d0+J, 1.0d0+JMA(2), DBLCMPLX(7.0d0 + ALPT(2), 0.0d0),
     &        DBLCMPLX(ALPT(2), 0.0d0) ) / SIN(PI*ALPT(2)) *
     &          (1.0d0 - (SIN(PI*ALPT(2)) -
     &           SIN(J*PI))/SIN(PI*JMA(2)) ) )
        ENDIF
      ELSE IF (ANSATZ .EQ. 'FITRAT') THEN
* Use first sum-rule constraint
          PAR(21) = 0.6d0 - PAR(11)
          ALPT(1) = PAR(12) + PAR(13)*DEL2
          ALPT(2) = PAR(22) + PAR(23)*DEL2
          JMA(1) = J - ALPT(1)
          JMA(2) = J - ALPT(2)
          CORR(1) = RATGAMMA(DBLCMPLX(7.0d0 + PAR(12), 0.0d0),
     &           DBLCMPLX(7.0d0 + ALPT(1), 0.0d0))
          CORR(2) = RATGAMMA(DBLCMPLX(7.0d0 + PAR(22), 0.0d0),
     &           DBLCMPLX(7.0d0 + ALPT(2), 0.0d0))
        IF ( PROCESS .EQ. 'DIS' ) THEN 
          FCM(1) = PAR(11) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(12),0.0d0) 
     &              , 9) / POCHHAMMER(1.0d0 - PAR(12)
     &          + J , 9) * (1.0d0 - PAR(12) + J) / (1.0d0 + JMA(1))
          FCM(2) = PAR(21) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(22),0.0d0)
     &              , 7) / POCHHAMMER(1.0d0 - PAR(22)
     &          + J , 7) * (1.0d0 - PAR(22) + J) / (1.0d0 + JMA(2))
        ELSE  
          FCM(1) = PAR(11) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(12),0.0d0)
     &              , 9) / POCHHAMMER(1.0d0 - PAR(12)
     &          + J , 9) * (1.0d0 - PAR(12) + J) / (1.0d0 + JMA(1)) /
     &           (1.0d0 - DEL2/(PAR(14)+PAR(15)*J))**PAR(16) *
     &        ( 1.0d0 - PAR(19)*PI*(XI/2.0d0)**(JMA(1)+1.0d0)*SQRT(PI)*
     &         2.0d0**(-2.0d0*ALPT(1))*RATGAMMA(1.0d0+J, 1.5d0+J)*
     &         RATGAMMA(1.0d0+J+ALPT(1), 1.0d0+JMA(1))*
     &          EXP(- 2.0d0*CLNGAMMA(
     &           DBLCMPLX(ALPT(1), 0.0d0)) ) / SIN(PI*ALPT(1)) *
     &          (1.0d0 - (SIN(PI*ALPT(1)) -
     &           SIN(J*PI))/SIN(PI*JMA(1)) ) * CORR(1)  )
          FCM(2) = PAR(21) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(22),0.0d0)
     &              , 7) / POCHHAMMER(1.0d0 - PAR(22)
     &          + J , 7) * (1.0d0 - PAR(22) + J) / (1.0d0 + JMA(2)) /
     &           (1.0d0 - DEL2/(PAR(24)+PAR(25)*J))**PAR(26) *
     &        ( 1.0d0 - PAR(29)*PI*(XI/2.0d0)**(JMA(2)+1.0d0)*SQRT(PI)*
     &         2.0d0**(-2.0d0*ALPT(2))*RATGAMMA(1.0d0+J, 1.5d0+J)*
     &         RATGAMMA(1.0d0+J+ALPT(2), 1.0d0+JMA(2))*
     &          EXP(- 2.0d0*CLNGAMMA(
     &           DBLCMPLX(ALPT(2), 0.0d0)) ) / SIN(PI*ALPT(2)) *
     &          (1.0d0 - (SIN(PI*ALPT(2)) -
     &           SIN(J*PI))/SIN(PI*JMA(2)) ) * CORR(2)  )
        ENDIF
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
        IF ( PROCESS .EQ. 'DIS' ) THEN 
          FCM(1) = PAR(11) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(12),0.0d0) 
     &              , 8) / POCHHAMMER(1.0d0 - PAR(12)
     &          + J , 8) * (1.0d0 - PAR(12) + J) / (1.0d0 + JMA(1))
          FCM(2) = PAR(21) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(22),0.0d0)
     &              , 6) / POCHHAMMER(1.0d0 - PAR(22)
     &          + J , 6) * (1.0d0 - PAR(22) + J) / (1.0d0 + JMA(2))
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
        ENDIF
      ELSE IF (ANSATZ .EQ. 'FITE86') THEN
* Use first sum-rule constraint
          PAR(21) = 0.6d0 - PAR(11)
          ALPT(1) = PAR(12) + PAR(13)*DEL2
          ALPT(2) = PAR(22) + PAR(23)*DEL2
          JMA(1) = J - ALPT(1)
          JMA(2) = J - ALPT(2)
          CORR(1) = EXP(   CLNGAMMA(DBLCMPLX(7.0d0 + PAR(12), 0.0d0)) -
     &           CLNGAMMA(DBLCMPLX(7.0d0 + ALPT(1), 0.0d0))   )
          CORR(2) = EXP(   CLNGAMMA(DBLCMPLX(7.0d0 + PAR(22), 0.0d0)) -
     &           CLNGAMMA(DBLCMPLX(7.0d0 + ALPT(2), 0.0d0))   )
        IF ( PROCESS .EQ. 'DIS' ) THEN 
          FCM(1) = PAR(11) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(12),0.0d0) 
     &              , 8) / POCHHAMMER(1.0d0 - PAR(12)
     &          + J , 8) * (1.0d0 - PAR(12) + J) / (1.0d0 + JMA(1))
          FCM(2) = PAR(21) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(22),0.0d0)
     &              , 6) / POCHHAMMER(1.0d0 - PAR(22)
     &          + J , 6) * (1.0d0 - PAR(22) + J) / (1.0d0 + JMA(2))
        ELSE  
          FCM(1) = PAR(11) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(12),0.0d0)
     &              , 8) / POCHHAMMER(1.0d0 - PAR(12)
     &          + J , 8) * (1.0d0 - PAR(12) + J) / (1.0d0 + JMA(1)) *
     &           EXP(DEL2 / (2.d0 * PAR(14))) *
     &        ( 1.0d0 - PAR(19)*PI*(XI/2.0d0)**(JMA(1)+1.0d0)*SQRT(PI)*
     &         2.0d0**(-2.0d0*ALPT(1))*EXP( CLNGAMMA(1.0d0+J) +
     &          CLNGAMMA(1.0d0+J+ALPT(1)) - CLNGAMMA(1.5d0+J) - 
     &          CLNGAMMA(1.0d0+JMA(1)) - 2.0d0*CLNGAMMA(
     &           DBLCMPLX(ALPT(1), 0.0d0)) ) / SIN(PI*ALPT(1)) *
     &          (1.0d0 - (SIN(PI*ALPT(1)) -
     &           SIN(J*PI))/SIN(PI*JMA(1)) ) * CORR(1)  )
          FCM(2) = PAR(21) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(22),0.0d0)
     &              , 6) / POCHHAMMER(1.0d0 - PAR(22)
     &          + J , 6) * (1.0d0 - PAR(22) + J) / (1.0d0 + JMA(2)) *
     &           EXP(DEL2 / (2.d0 * PAR(24))) *
     &        ( 1.0d0 - PAR(29)*PI*(XI/2.0d0)**(JMA(2)+1.0d0)*SQRT(PI)*
     &         2.0d0**(-2.0d0*ALPT(2))*EXP( CLNGAMMA(1.0d0+J) +
     &          CLNGAMMA(1.0d0+J+ALPT(2)) - CLNGAMMA(1.5d0+J) - 
     &          CLNGAMMA(1.0d0+JMA(2)) - 2.0d0*CLNGAMMA(
     &           DBLCMPLX(ALPT(2), 0.0d0)) ) / SIN(PI*ALPT(2)) *
     &          (1.0d0 - (SIN(PI*ALPT(2)) -
     &           SIN(J*PI))/SIN(PI*JMA(2)) ) * CORR(2) )
        ENDIF
      ELSE IF (ANSATZ .EQ. 'FITEXP') THEN
* Use first sum-rule constraint
          PAR(21) = 0.6d0 - PAR(11)
          ALPT(1) = PAR(12) + PAR(13)*DEL2
          ALPT(2) = PAR(22) + PAR(23)*DEL2
          JMA(1) = J - ALPT(1)
          JMA(2) = J - ALPT(2)
          CORR(1) = EXP(   CLNGAMMA(DBLCMPLX(7.0d0 + PAR(12), 0.0d0)) -
     &           CLNGAMMA(DBLCMPLX(7.0d0 + ALPT(1), 0.0d0))   )
          CORR(2) = EXP(   CLNGAMMA(DBLCMPLX(7.0d0 + PAR(22), 0.0d0)) -
     &           CLNGAMMA(DBLCMPLX(7.0d0 + ALPT(2), 0.0d0))   )
        IF ( PROCESS .EQ. 'DIS' ) THEN 
          FCM(1) = PAR(11) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(12),0.0d0) 
     &              , 9) / POCHHAMMER(1.0d0 - PAR(12)
     &          + J , 9) * (1.0d0 - PAR(12) + J) / (1.0d0 + JMA(1))
          FCM(2) = PAR(21) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(22),0.0d0)
     &              , 7) / POCHHAMMER(1.0d0 - PAR(22)
     &          + J , 7) * (1.0d0 - PAR(22) + J) / (1.0d0 + JMA(2))
        ELSE  
          FCM(1) = PAR(11) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(12),0.0d0)
     &              , 9) / POCHHAMMER(1.0d0 - PAR(12)
     &          + J , 9) * (1.0d0 - PAR(12) + J) / (1.0d0 + JMA(1)) *
     &           EXP(DEL2 / (2.d0 * PAR(14))) *
     &        ( 1.0d0 - PAR(19)*PI*(XI/2.0d0)**(JMA(1)+1.0d0)*SQRT(PI)*
     &         2.0d0**(-2.0d0*ALPT(1))*EXP( CLNGAMMA(1.0d0+J) +
     &          CLNGAMMA(1.0d0+J+ALPT(1)) - CLNGAMMA(1.5d0+J) - 
     &          CLNGAMMA(1.0d0+JMA(1)) - 2.0d0*CLNGAMMA(
     &           DBLCMPLX(ALPT(1), 0.0d0)) ) / SIN(PI*ALPT(1)) *
     &          (1.0d0 - (SIN(PI*ALPT(1)) -
     &           SIN(J*PI))/SIN(PI*JMA(1)) ) * CORR(1)  )
          FCM(2) = PAR(21) * POCHHAMMER(DBLCMPLX(2.0d0 - PAR(22),0.0d0)
     &              , 7) / POCHHAMMER(1.0d0 - PAR(22)
     &          + J , 7) * (1.0d0 - PAR(22) + J) / (1.0d0 + JMA(2)) *
     &           EXP(DEL2 / (2.d0 * PAR(24))) *
     &        ( 1.0d0 - PAR(29)*PI*(XI/2.0d0)**(JMA(2)+1.0d0)*SQRT(PI)*
     &         2.0d0**(-2.0d0*ALPT(2))*EXP( CLNGAMMA(1.0d0+J) +
     &          CLNGAMMA(1.0d0+J+ALPT(2)) - CLNGAMMA(1.5d0+J) - 
     &          CLNGAMMA(1.0d0+JMA(2)) - 2.0d0*CLNGAMMA(
     &           DBLCMPLX(ALPT(2), 0.0d0)) ) / SIN(PI*ALPT(2)) *
     &          (1.0d0 - (SIN(PI*ALPT(2)) -
     &           SIN(J*PI))/SIN(PI*JMA(2)) ) * CORR(2) )
        ENDIF
      ELSE IF (ANSATZ .EQ. 'SPLICE') THEN
              CALL SPLICE(J, FCM)
*           Killing E contribution
            PAR(18) = 0.0d0
            PAR(28) = 0.0d0
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
*           Killing E contribution
            PAR(18) = 0.0d0
            PAR(28) = 0.0d0
      ELSE IF (ANSATZ .EQ. 'NSFIT') THEN
*    Ugly kludge: NS valence = down
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
          FCM(1) = HD - 0.5d0 / (2.0d0 + 0.5d0) * HSEA
*    'Dummy' gluonic, not used:
          FCM(2) = (0.0d0, 0.0d0)
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
      END IF

      RETURN
      END
C     ****
