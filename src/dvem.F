C     ****h* gepard/dvem.F
C  FILE DESCRIPTION
C    calculation of  Wilson coefficients for DVEM in MSbar scheme
C    according to DVEM-notes.pdf
C
C    $Id$
C     *******


C     ****s* dvem.F/MCQ
C  NAME
C     MCQ       --  "big C" Wilson coefficient for DVEM - quark contribution
C  DESCRIPTION
C    calculates Wilson coefficients for  ...
C  SYNOPSIS

      DOUBLE COMPLEX FUNCTION MCQ (J, K)

      IMPLICIT NONE
      DOUBLE COMPLEX J, K, JAUX
      INTEGER PJ, PK, PJAUX, TWO

C  SOURCE
C

      DOUBLE COMPLEX MCMUFQ1F, MCMURQ1F, MCBETQ1F, MCFQ1F, MCGAQ1F
      DOUBLE COMPLEX MCGBIAQ1F, MCGBIBQ1F, MCGBIIQ1F, MCSINGQ1F
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

      PJ = 1
      PK = 2
      RGPDF2 = 2.0d0
      RDAF2 = 2.1d0
      RR2 = 2.2d0
      NF = 3
      CALL BETAF

      LRGPDF2 = LOG(RGPDF2)
      LRDAF2 = LOG(RDAF2)
      LRR2 = LOG(RR2)

*  Spliced from Mathematica:

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
          PJAUX = PJ
          JAUX =J
          PJ = PK
          J = K
          PK = PJAUX
          K = JAUX
      END IF

        MCGBIAQ1F=CG*(-4.289868133696452-4.*HSM2(1.+K,1.+PK)+2.*
     &  (1.+K)*(2.+K)*(ZETA3+3.2898681336964524*HS1(1.+K)+4.*HS1
     &  M2(1.+K,1.+PK)-2.*HS3(1.+K)-2.*HSM3(1.+K,1.+PK))-2.*((1.
     &  +K)**(-2)+(2.+K)**(-2))*(1.+PTY(PK))+(4.+PTY(PJ)+PTY(PK)
     &  )/((1.+K)*(2.+K))-2.*(1/((1.+J)*(2.+J))-2.*HS1(1.+J))*(1
     &  .-PTY(PK)/((1.+K)*(2.+K))))

      MCQ = MCQ + MCGBIAQ1F 

        MCGBIBQ1F=2.*CG*((1.-1/((1.+J)*(2.+J))+2.*HS1(1.+J))*(1.
     &  -0.5*(1.+K)*(2.+K)*(1.6449340668482262-HS2(0.5*K))+0.5*(
     &  1.+K)*(2.+K)*(1.6449340668482262-HS2(0.5*(1.+K))))-2.*(1
     &  .+K)*(2.+K)*(-0.25*ZETA3+3.2898681336964524*HS1(1.+K)+4.
     &  *HS1M2(1.+K,PK)-HSM3(1.+K,PK))*PTY(K))*PTY(PK)

      MCQ = MCQ + MCGBIBQ1F 

        MCGBIIQ1F=CG*(-2.*(1.+J)*(2.+J)*(1.+K)*(2.+K)*FINT(0.5,J
     &  ,0.5,K,3.)*NCONF(0.5,K)*NLNCONF(0.5,J)+(1.+J)*(2.+J)*FIN
     &  T(0.5,J,1.5,K,1.)*NCONF(1.5,K)*NLNCONF(0.5,J)-(1.+K)*(2.
     &  +K)*FINT(1.5,J,0.5,K,1.)*NCONF(0.5,K)*NLNCONF(1.5,J)-2.*
     &  FINT(1.5,J,1.5,K,2.)*NCONF(1.5,K)*NLNCONF(1.5,J))*PTY(PK
     &  )

      MCQ = MCQ + MCGBIIQ1F

*     PRINT *, 'LI, MCGBIAQ1F = ', LI,  MCGBIAQ1F
*     PRINT *, 'LI, MCGB1BQ1F = ', LI,  MCGBIBQ1F
*     PRINT *, 'LI, MCGBIIQ1F = ', LI,  MCGBIIQ1F

 10   CONTINUE

        MCSINGQ1F=(2.*CF*(1.+2./((1.+J)*(2.+J)))*LRGPDF2)/(J*(3.
     &  +J))+2.*CF*(-(1/(J*(1.+J)*(2.+J)*(3.+J)))+((1.+2./((1.+J
     &  )*(2.+J)))*(1.+1/((1.+J)*(2.+J))+3./(J*(3.+J))-2.*HS1(1.
     &  +J)-HS1(1.+K)))/(J*(3.+J)))+CF*(-1.6449340668482262-1/((
     &  1.+J)*(2.+J)*(1.+K)*(2.+K))-(8.*(3.+3.*J+J*J))/(J**2*(1.
     &  +J)*(2.+J)*(3.+J)**2)+2.*(1/((1.+J)*(2.+J))-2./(J*(3.+J)
     &  ))*HS1(1.+K)-2.*HSM2(1.+J,2.))+CF*(-(FINT(1.5,J,1.5,K,4.
     &  )*NCONF(1.5,K)*NLNCONF(1.5,J))-0.5*(1.+K)*(2.+K)*(FINT(0
     &  .5,K,1.5,J,1.)*NCONF(1.5,J)*NLNCONF(0.5,K)+FINT(1.5,J,0.
     &  5,K,1.)*NCONF(0.5,K)*NLNCONF(1.5,J))+FINT(1.5,K,1.5,J,4.
     &  )*NCONF(1.5,J)*NLNCONF(1.5,K))

      PRINT *, ' MCSINGQ1F = ',  MCSINGQ1F

        MCG1F=LRR2*BETA0(NF)+LRGPDF2*(2.*CF*(0.5+1/((1.+J)*(2.+J
     &  )))-BETA0(NF)+2.*CA*(-2./((1.+J)*(2.+J))+6./(J*(3.+J))-2
     &  .*HS1(1.+J)))+2.*CF*LRDAF2*(1.5+1/((1.+K)*(2.+K))-2.*HS1
     &  (1.+K))+CF*(-5.+2./((1.+J)**2*(2.+J)**2)+7./((1.+J)*(2.+
     &  J))+2./((1.+K)**2*(2.+K)**2)+5./((1.+K)*(2.+K))+2./((1.+
     &  J)*(2.+J)*(1.+K)*(2.+K))-2.*HS1(1.+J)*(3.+2./((1.+J)*(2.
     &  +J))+1/((1.+K)*(2.+K))-2.*HS1(1.+K))-2.*(3.+2./((1.+J)*(
     &  2.+J))+2./((1.+K)*(2.+K)))*HS1(1.+K)+4.*HS1(1.+K)**2)+2.
     &  *CA*(1.+4./((1.+J)*(2.+J))+36./(J**2*(1.+J)*(2.+J)*(3.+J
     &  )**2)+4./(J*(3.+J))+4./(J*(1.+J)**2*(2.+J)**2*(3.+J))+30
     &  ./(J*(1.+J)*(2.+J)*(3.+J))+2.*HS1(1.+J)**2+2.*HS1(1.+J)*
     &  HS1(1.+K)-2.*(1/((1.+J)*(2.+J))+1/(J*(3.+J))+4./(J*(1.+J
     &  )*(2.+J)*(3.+J)))*(2.*HS1(1.+J)+HS1(1.+K)))+CA*(5./(J*(3
     &  .+J))-(6.*(3.+2.*J))/(J**2*(3.+J)**2)-(J*(5.+3.*J))/((1.
     &  +J)**2*(2.+J)**2)+HS2(0.5*(2.+J))-HS2(0.5*(3.+J))+0.25*(
     &  1.+K)*(2.+K)*(HS2(0.5*K)-HS2(0.5*(1.+K))))+CA*(-2./J**2-
     &  2./(1.+J)**2-2./(2.+J)**2-2./(3.+J)**2+0.5*(-1.-3./((1.+
     &  J)*(2.+J))+1.6449340668482262*(1.+J)*(2.+J)+7./(J*(3.+J)
     &  ))+(6.+3.*J+J*J)/(J*(1.+J)*(2.+J)*(3.+J)*(1.+K)*(2.+K))-
     &  (8.*(3.+3.*J+J*J)*HS1(1.+K))/(J*(1.+J)*(2.+J)*(3.+J))+2.
     &  *HS1(1.+J)*(1.+HS1(1.+K))+J*(3.+J)*HSM2(1.+J,2.))+CF*(-2
     &  .644934066848226-2./((1.+K)**2*(2.+K)**2)-(1.-1/((1.+J)*
     &  (2.+J))-3./(J*(3.+J))+2.*HS1(1.+J))/((1.+K)*(2.+K))-2.*H
     &  SM2(1.+K,1.))+CF*(0.25*J*(3.+J)*(1.+K)*(2.+K)*(FINT(0.5,
     &  K,1.5,J,1.)*NCONF(1.5,J)*NLNCONF(0.5,K)+FINT(1.5,J,0.5,K
     &  ,1.)*NCONF(0.5,K)*NLNCONF(1.5,J))-FINT(1.5,K,2.5,J,5.)*N
     &  CONF(2.5,J)*NLNCONF(1.5,K)-FINT(2.5,J,1.5,K,5.)*NCONF(1.
     &  5,K)*NLNCONF(2.5,J))-0.5*CA*(FINT(1.5,K,2.5,J,2.)*NCONF(
     &  2.5,J)*NLNCONF(1.5,K)+FINT(2.5,J,1.5,K,2.)*NCONF(1.5,K)*
     &  NLNCONF(2.5,J)+(1.+K)*(2.+K)*(FINT(0.5,K,2.5,J,5.)*NCONF
     &  (2.5,J)*NLNCONF(0.5,K)-FINT(2.5,J,0.5,K,5.)*NCONF(0.5,K)
     &  *NLNCONF(2.5,J)))

      PRINT *, ' MCG1F = ',  MCG1F

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
