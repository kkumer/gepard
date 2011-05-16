C     ****h* gepard/dvem.F
C  FILE DESCRIPTION
C    calculation of  Wilson coefficients for DVEM in MSbar scheme
C    according to DVEM-notes.pdf
C
C    $Id$
C     *******


C     ****s* dvem.F/MC
C  NAME
C     MC       --  "big C" Wilson coefficient for DVEM
C  DESCRIPTION
C    calculates Wilson coefficients for  ...
C  SYNOPSIS

      DOUBLE COMPLEX FUNCTION MC (J, K)

      IMPLICIT NONE
      DOUBLE COMPLEX J, K
      INTEGER PJ, PK, TWO

C  SOURCE
C

      DOUBLE COMPLEX MCMUFQ1F, MCMURQ1F, MCBETQ1F, MCFQ1F, MCGAQ1F
      DOUBLE COMPLEX MCGBIAQ1F, MCGBIBQ1F, MCGBIIQ1F
      DOUBLE COMPLEX NCONF, NLNCONF, FINT
      DOUBLE PRECISION LRGPDF2, LRDAF2, LRR2
      DOUBLE COMPLEX HS1, HS2, HS3, HSM1, HSM2, HSM3, HS1M2
      EXTERNAL HS1, HS2, HS3, HSM1, HSM2, HSM3, HS1M2
      INTEGER PTY
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

        MCGBIAQ1F=CG*(-4.289868133696452-4.*HSM2(1.+K,1.+PK)+2.*(1.+K)*(
     &  2.+K)*(ZETA3+3.2898681336964524*HS1(1.+K)+4.*HS1M2(1.+K,1.+PK)-2
     &  .*HS3(1.+K)-2.*HSM3(1.+K,1.+PK))-2.*((1.+K)**(-2)+(2.+K)**(-2))*
     &  (1.+PTY(PK))+(4.+PTY(PJ)+PTY(PK))/((1.+K)*(2.+K))-2.*(1/((1.+J)*
     &  (2.+J))-2.*HS1(1.+J))*(1.-PTY(PK)/((1.+K)*(2.+K))))

        MCGBIBQ1F=2.*CG*((1.-1/((1.+J)*(2.+J))+2.*HS1(1.+J))*(1.-0.5*(1.
     &  +K)*(2.+K)*(1.6449340668482262-HS2(0.5*K))+0.5*(1.+K)*(2.+K)*(1.
     &  6449340668482262-HS2(0.5*(1.+K))))-2.*(1.+K)*(2.+K)*(-0.25*ZETA3
     &  +3.2898681336964524*HS1(1.+K)+4.*HS1M2(1.+K,PK)-HSM3(1.+K,PK))*P
     &  TY(K))*PTY(PK)

        MCGBIIQ1F=CG*(-2.*(1.+J)*(2.+J)*(1.+K)*(2.+K)*FINT(0.5,J,0.5,K,3
     &  .)*NCONF(0.5,K)*NLNCONF(0.5,J)+(1.+J)*(2.+J)*FINT(0.5,J,1.5,K,1.
     &  )*NCONF(1.5,K)*NLNCONF(0.5,J)-(1.+K)*(2.+K)*FINT(1.5,J,0.5,K,1.)
     &  *NCONF(0.5,K)*NLNCONF(1.5,J)-2.*FINT(1.5,J,1.5,K,2.)*NCONF(1.5,K
     &  )*NLNCONF(1.5,J))*PTY(PK)

      PRINT *, 'MCMUFQ1F = ', MCMUFQ1F
      PRINT *, 'MCMURQ1F = ', MCMURQ1F
      PRINT *, 'MCBETQ1F = ', MCBETQ1F
      PRINT *, 'MCFQ1F = ', MCFQ1F
      PRINT *, 'MCGAQ1F = ', MCGAQ1F
      PRINT *, 'MCGBIAQ1F = ', MCGBIAQ1F
      PRINT *, 'MCGB1BQ1F = ', MCGBIBQ1F
      PRINT *, 'MCGBIIQ1F = ', MCGBIIQ1F

      MC = MCMUFQ1F + MCMURQ1F + MCBETQ1F + MCFQ1F + MCGAQ1F
     &  + MCGBIAQ1F + MCGBIBQ1F + MCGBIIQ1F

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

      WRITE (*,*) 'MU,NINT(MU) = ', MU, NINT(MU)
      IF (NINT(2*MU) .EQ. 1) THEN
        NLNCONF = 1.0d0/(2.0d0 + J)
      ELSE IF (NINT(2*MU) .EQ. 3) THEN
        NLNCONF = (J + 2.0d0)/(J + 1.0d0)
      ELSE IF (NINT(2*MU) .EQ. 5) THEN
        NLNCONF = (J + 2.0d0) * (J + 3.0d0) / J
      ELSE
        NLNCONF = 1.0d0
      END IF
      
      WRITE (*,*) 'NLNCONF = ', NLNCONF
      RETURN
      END

      DOUBLE COMPLEX FUNCTION NCONF(NU, K)
      REAL NU
      DOUBLE COMPLEX K, CLNGAMMA
      EXTERNAL CLNGAMMA

      NCONF =  CLNGAMMA(K+2.0d0) - CLNGAMMA(2.0d0*K + 4.0d0) +
     &         CLNGAMMA(K + 1.5d0 + NU) 
      WRITE (*,*) 'NCONF = ', NCONF
      NCONF = EXP(NCONF)
      WRITE (*,*) 'NCONF = ', NCONF
      RETURN
      END
