C     ****h* gepard/dvem.F
C  FILE DESCRIPTION
C    calculation of NLO Wilson coefficients for DVEM 
C    according to DVEM-notes.pdf
C
C    $Id$
C     *******


C     ****s* dvem.F/CDVEMF
C  NAME
C     CDVEMF  --  "big C" NLO Wilson coefficient for DVEM
C  DESCRIPTION
C    calculates (unnormalized) Wilson coefficients for DVEM
C    i.e. they are normalized w.r.t. LO C=1.
C    Presently parities are fixed (to vector production case)
C    OUTPUT MCU(K) K=1,2,3 for quark, pure singlet quark and gluon contrib
C  SYNOPSIS

      SUBROUTINE CDVEMF (SGNTR, J, K, MCU)

!f2py integer intent(in) :: SGNTR
!f2py complex(8) intent(in) :: J, K
!f2py complex(8) intent(out) :: MCU

      IMPLICIT NONE
      INTEGER SGNTR
      DOUBLE COMPLEX MCU(3), MCQ, J, K, JAUX, CLNGAMMA
*     DOUBLE COMPLEX NQNORM, NQBNORM, NGNORM
      INTEGER PJ, PK, PJAUX, TWO
      EXTERNAL CLNGAMMA

C  SOURCE
C

      DOUBLE COMPLEX MCMUFQ1F, MCMURQ1F, MCBETQ1F, MCFQ1F, MCGAQ1F
      DOUBLE COMPLEX MCBMUFQ1F, MCBAQ1F, MCBBQ1F
      DOUBLE COMPLEX MCGBIAQ1F, MCGBIBQ1F, MCGBIIQ1F
      DOUBLE COMPLEX MCG1F
      DOUBLE COMPLEX NCONF, NLNCONF, FINT
      DOUBLE PRECISION LRGPDF2, LRDAF2, LRR2
      DOUBLE COMPLEX HS1, HS2, HS3, HSM1, HSM2, HSM3, HS1M2
      EXTERNAL HS1, HS2, HS3, HSM1, HSM2, HSM3, HS1M2
      INTEGER PTY, LI
      EXTERNAL PTY, NCONF, NLNCONF, FINT
      PARAMETER (TWO=2)
#include "header.F"

C     J = N(1, K) - 1

      PJ = (1+SGNTR)/2
      PK = 2
      CALL BETAF

      LRGPDF2 = LOG(RGPDF2)
      LRDAF2 = LOG(RDAF2)
      LRR2 = LOG(RR2)

*      NQNORM=(0.5641895835477563*2.**(3.+J)*CF*(3.+2.*K)*
*     &  (-1.+PTY(PJ+PK))*PTY(PJ+PK))/((1.+K)*(2.+K)) *
*     &  EXP(CLNGAMMA(J+2.5d0) - CLNGAMMA(J+3.0d0))
*
*      NQBNORM=(0.5641895835477563*2.**(5.+J)*(3.+2.*K))/((1.+K)*(2.+K))*
*     &  EXP(CLNGAMMA(J+2.5d0) - CLNGAMMA(J+3.0d0))
*
*      NGNORM=(0.5641895835477563*2.**(5.+J)*(3.+2.*K))/((
*     &  3.+J)*(1.+K)*(2.+K))*
*     &  EXP(CLNGAMMA(J+2.5d0) - CLNGAMMA(J+3.0d0))

*  Spliced from Mathematica:

*  ... quark part

        MCMUFQ1F=2.d0*CF*LRGPDF2*(1.5d0+1/((1.d0+J)*(2.d0+J))-2.d0*HS1(1
     &  .d0+J))+2.d0*CF*LRDAF2*(1.5d0+1/((1.d0+K)*(2.d0+K))-2.d0*HS1(1.d
     &  0+K))

        MCMURQ1F=LRR2*BETA0(NF)

        MCBETQ1F=BETA0(NF)*(-1.6666666666666667+1/((1.+J)*(2.+J))+1/((1.
     &  +K)*(2.+K))-2.*HS1(1.+J)-2.*HS1(1.+K))

        MCFQ1F=CF*(-9.333333333333334+4./((1.+J)*(2.+J))+4./((1.+K)*(2.+
     &  K))+(1/(2.+J)-1/(1.+K))**2+(1/(1.+J)+1/(1.+K))**2+(1/(1.+J)-1/(2
     &  .+K))**2+(1/(2.+J)+1/(2.+K))**2+2.*(HS1(1.+J)+HS1(1.+K))*(-5.-2.
     &  /((1.+J)*(2.+J))-2./((1.+K)*(2.+K))+2.*HS1(1.+J)+2.*HS1(1.+K)))

        MCGAQ1F=2.*CG*(-3.3333333333333335-2.*HSM2(1.+J,1.+PJ)-2.*HSM2(1
     &  .+K,1.+PK)+(2.*HS1(1.+K)*PTY(PJ))/((1.+J)*(2.+J))-((1.+1/((1.+J)
     &  *(2.+J)))*(1.+PTY(PJ)))/((1.+J)*(2.+J))+(2.*HS1(1.+J)*PTY(PK))/(
     &  (1.+K)*(2.+K))-((1.+1/((1.+K)*(2.+K)))*(1.+PTY(PK)))/((1.+K)*(2.
     &  +K))-(PTY(PJ)+PTY(PK)+2.*PTY(PJ+PK))/((1.+J)*(2.+J)*(1.+K)*(2.+K
     &  )))

*     PRINT *, 'MCMUFQ1F = ', MCMUFQ1F
*     PRINT *, 'MCMURQ1F = ', MCMURQ1F
*     PRINT *, 'MCBETQ1F = ', MCBETQ1F
*     PRINT *, 'MCFQ1F = ', MCFQ1F
*     PRINT *, 'MCGAQ1F = ', MCGAQ1F

      MCQ = MCMUFQ1F + MCMURQ1F + MCBETQ1F + MCFQ1F + MCGAQ1F

      DO 10 LI=1,2
      IF (LI .EQ. 2) THEN
*      J <--> K
          PJAUX = PJ
          JAUX =J
          PJ = PK
          J = K
          PK = PJAUX
          K = JAUX
      END IF

        MCGBIAQ1F=CG*(-4.289868133696452-4.*HSM2(1.+K,1.+PK)+2.*(1.+K)*(
     &  2.+K)*(ZETA3+3.2898681336964524*HS1(1.+K)+4.*HS1M2(1.+K,1.+PK)-2
     &  .*HS3(1.+K)-2.*HSM3(1.+K,1.+PK))-2.*((1.+K)**(-2)+(2.+K)**(-2))*
     &  (1.+PTY(PK))+(4.+PTY(PJ)+PTY(PK))/((1.+K)*(2.+K))-2.*(1/((1.+J)*
     &  (2.+J))-2.*HS1(1.+J))*(1.-PTY(PK)/((1.+K)*(2.+K))))

      MCQ = MCQ + MCGBIAQ1F 

        MCGBIBQ1F=2.*CG*PTY(PK)*((1.-1/((1.+J)*(2.+J))+2.*HS1(1.+J))*(1.
     &  -0.5*(1.+K)*(2.+K)*(1.6449340668482262-HS2(0.5*K))+0.5*(1.+K)*(2
     &  .+K)*(1.6449340668482262-HS2(0.5*(1.+K))))-2.*(1.+K)*(2.+K)*(-ZE
     &  TA3+3.2898681336964524*HS1(1.+K)+4.*HS1M2(1.+K,1.+PK)-2.*HSM3(1.
     &  +K,1.+PK))*PTY(PK))

      MCQ = MCQ + MCGBIBQ1F 

        MCGBIIQ1F=CG*(-2.*(1.+J)*(2.+J)*(1.+K)*(2.+K)*FINT(0.5,J,0.5,K,3
     &  .)*NCONF(0.5,K)*NLNCONF(0.5,J)+(1.+J)*(2.+J)*FINT(0.5,J,1.5,K,1.
     &  )*NCONF(1.5,K)*NLNCONF(0.5,J)-(1.+K)*(2.+K)*FINT(1.5,J,0.5,K,1.)
     &  *NCONF(0.5,K)*NLNCONF(1.5,J)-2.*FINT(1.5,J,1.5,K,2.)*NCONF(1.5,K
     &  )*NLNCONF(1.5,J))*PTY(PK)

      MCQ = MCQ + MCGBIIQ1F

*     PRINT *, 'LI, MCGBIAQ1F = ', LI,  MCGBIAQ1F
*     PRINT *, 'LI, MCGB1BQ1F = ', LI,  MCGBIBQ1F
*     PRINT *, 'LI, MCGBIIQ1F = ', LI,  MCGBIIQ1F

 10   CONTINUE
 
*     PRINT *, ' NQNORM = ',  NQNORM
*     PRINT *, ' unnorm Q = ',  MCQ
*     change J and K back to their values
        PJAUX = PJ
        JAUX =J
        PJ = PK
        J = K
        PK = PJAUX
        K = JAUX

C  ... quark part starting at NLO

        MCBMUFQ1F=(2.*CF*(1.+2./((1.+J)*(2.+J)))*LRGPDF2)/(J*(3.+J))

        MCBAQ1F=2.*CF*(-(1/(J*(1.+J)*(2.+J)*(3.+J)))+((1.+2./((1.+J)*(2.
     &  +J)))*(1.+1/((1.+J)*(2.+J))+3./(J*(3.+J))-2.*HS1(1.+J)-HS1(1.+K)
     &  ))/(J*(3.+J)))

        MCBBQ1F=CF*(-1.6449340668482262-1/((1.+J)*(2.+J)*(1.+K)*(2.+K))-
     &  (8.*(3.+3.*J+J*J))/(J**2*(1.+J)*(2.+J)*(3.+J)**2)+2.*(1/((1.+J)*
     &  (2.+J))-2./(J*(3.+J)))*HS1(1.+K)-2.*HSM2(1.+J,2.))+CF*(-(FINT(1.
     &  5,J,1.5,K,4.)*NCONF(1.5,K)*NLNCONF(1.5,J))-0.5*(1.+K)*(2.+K)*(FI
     &  NT(0.5,K,1.5,J,1.)*NCONF(1.5,J)*NLNCONF(0.5,K)+FINT(1.5,J,0.5,K,
     &  1.)*NCONF(0.5,K)*NLNCONF(1.5,J))+FINT(1.5,K,1.5,J,4.)*NCONF(1.5,
     &  J)*NLNCONF(1.5,K))

*     PRINT *, ' MCBMUFQ1F = ',  MCBMUFQ1F
*     PRINT *, ' MCBAQ1F = ',  MCBAQ1F
*     PRINT *, ' MCBBQ1F = ',  MCBBQ1F

*     PRINT *, ' NQBNORM = ',  NQBNORM
*     PRINT *, ' unnorm NLO QUARK = ',  (MCBMUFQ1F+MCBAQ1F+MCBBQ1F)
*     PRINT *, ' NLO QUARK = ',  NQBNORM * (MCBMUFQ1F+MCBAQ1F+MCBBQ1F)

C   ... gluon part

        MCG1F=LRR2*BETA0(NF)+LRGPDF2*(2.*CF*(0.5+1/((1.+J)*(2.+J)))-BETA
     &  0(NF)+2.*CA*(-2./((1.+J)*(2.+J))+6./(J*(3.+J))-2.*HS1(1.+J)))+2.
     &  *CF*LRDAF2*(1.5+1/((1.+K)*(2.+K))-2.*HS1(1.+K))+CF*(-5.+2./((1.+
     &  J)**2*(2.+J)**2)+7./((1.+J)*(2.+J))+2./((1.+K)**2*(2.+K)**2)+5./
     &  ((1.+K)*(2.+K))+2./((1.+J)*(2.+J)*(1.+K)*(2.+K))-2.*HS1(1.+J)*(3
     &  .+2./((1.+J)*(2.+J))+1/((1.+K)*(2.+K))-2.*HS1(1.+K))-2.*(3.+2./(
     &  (1.+J)*(2.+J))+2./((1.+K)*(2.+K)))*HS1(1.+K)+4.*HS1(1.+K)**2)+2.
     &  *CA*(1.+4./((1.+J)*(2.+J))+36./(J**2*(1.+J)*(2.+J)*(3.+J)**2)+4.
     &  /(J*(3.+J))+4./(J*(1.+J)**2*(2.+J)**2*(3.+J))+30./(J*(1.+J)*(2.+
     &  J)*(3.+J))+2.*HS1(1.+J)**2+2.*HS1(1.+J)*HS1(1.+K)-2.*(1/((1.+J)*
     &  (2.+J))+1/(J*(3.+J))+4./(J*(1.+J)*(2.+J)*(3.+J)))*(2.*HS1(1.+J)+
     &  HS1(1.+K)))+CA*(5./(J*(3.+J))-(6.*(3.+2.*J))/(J**2*(3.+J)**2)-(J
     &  *(5.+3.*J))/((1.+J)**2*(2.+J)**2)+HS2(0.5*(2.+J))-HS2(0.5*(3.+J)
     &  )+0.25*(1.+K)*(2.+K)*(HS2(0.5*K)-HS2(0.5*(1.+K))))+CA*(-2./J**2-
     &  2./(1.+J)**2-2./(2.+J)**2-2./(3.+J)**2+0.5*(-1.-3./((1.+J)*(2.+J
     &  ))+1.6449340668482262*(1.+J)*(2.+J)+7./(J*(3.+J)))+(6.+3.*J+J*J)
     &  /(J*(1.+J)*(2.+J)*(3.+J)*(1.+K)*(2.+K))-(8.*(3.+3.*J+J*J)*HS1(1.
     &  +K))/(J*(1.+J)*(2.+J)*(3.+J))+2.*HS1(1.+J)*(1.+HS1(1.+K))+J*(3.+
     &  J)*HSM2(1.+J,2.))+CF*(-2.644934066848226-2./((1.+K)**2*(2.+K)**2
     &  )-(1.-1/((1.+J)*(2.+J))-3./(J*(3.+J))+2.*HS1(1.+J))/((1.+K)*(2.+
     &  K))-2.*HSM2(1.+K,1.))+CF*(0.25*J*(3.+J)*(1.+K)*(2.+K)*(FINT(0.5,
     &  K,1.5,J,1.)*NCONF(1.5,J)*NLNCONF(0.5,K)+FINT(1.5,J,0.5,K,1.)*NCO
     &  NF(0.5,K)*NLNCONF(1.5,J))-FINT(1.5,K,2.5,J,5.)*NCONF(2.5,J)*NLNC
     &  ONF(1.5,K)-FINT(2.5,J,1.5,K,5.)*NCONF(1.5,K)*NLNCONF(2.5,J))-0.5
     &  *CA*(FINT(1.5,K,2.5,J,2.)*NCONF(2.5,J)*NLNCONF(1.5,K)+FINT(2.5,J
     &  ,1.5,K,2.)*NCONF(1.5,K)*NLNCONF(2.5,J)+(1.+K)*(2.+K)*(FINT(0.5,K
     &  ,2.5,J,5.)*NCONF(2.5,J)*NLNCONF(0.5,K)-FINT(2.5,J,0.5,K,5.)*NCON
     &  F(0.5,K)*NLNCONF(2.5,J)))

*     PRINT *, ' NGNORM = ',  NGNORM
*     PRINT *, ' unnorm G = ',  MCG1F

*   dividing all by 2 to get to (alphas/2pi) expansion
      MCU(1) = MCQ / 2.d0
*   dividing pure singlet by 2/CF for agreement with DM
      MCU(2) = (MCBMUFQ1F+MCBAQ1F+MCBBQ1F)/ 2.d0 * (2.d0/CF)
      MCU(3) = MCG1F / 2.d0

*     MC(1) = NQNORM * MCU(1)
*     MC(2) = NQBNORM * MCU(2)
*     MC(3) = NGNORM * MCU(3)

      RETURN
      END
C     *******



      INTEGER FUNCTION PTY(K)
       
      IMPLICIT NONE
      INTEGER K

      PTY = (1 - 2*MOD(K, 2))
      RETURN
      END

      DOUBLE COMPLEX FUNCTION NLNCONF(MU, J)

      IMPLICIT NONE
      REAL MU
      DOUBLE COMPLEX J, CLNGAMMA
      EXTERNAL CLNGAMMA

*     WRITE (*,*) 'MU,NINT(2*MU) = ', MU, NINT(2*MU)
      IF (NINT(2*MU) .EQ. 1) THEN
        NLNCONF = 1.0d0/(2.0d0 + J)
      ELSE IF (NINT(2*MU) .EQ. 3) THEN
        NLNCONF = (J + 2.0d0)/(J + 1.0d0)
      ELSE IF (NINT(2*MU) .EQ. 5) THEN
        NLNCONF = (J + 2.0d0) * (J + 3.0d0) / J
      ELSE
        NLNCONF = 1.0d0
      END IF
      
*     WRITE (*,*) 'MU, J, NLNCONF = ', MU, J, NLNCONF
      RETURN
      END

      DOUBLE COMPLEX FUNCTION NCONF(NU, K)
      REAL NU
      DOUBLE COMPLEX K, CLNGAMMA
      EXTERNAL CLNGAMMA

      NCONF =  CLNGAMMA(K+2.0d0) - CLNGAMMA(2.0d0*K + 4.0d0) +
     &         CLNGAMMA(K + 1.5d0 + NU) 
      NCONF = EXP(NCONF)
*     WRITE (*,*) 'NU, K, NCONF = ', NU, K, NCONF
      RETURN
      END
