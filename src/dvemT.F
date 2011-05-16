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

<* MCMUFQ1Fexpr *>

<* MCMURQ1Fexpr *>

<* MCBETQ1Fexpr *>

<* MCFQ1Fexpr *>

<* MCGAQ1Fexpr *>

<* MCGBIAQ1Fexpr *>

<* MCGBIBQ1Fexpr *>

<* MCGBIIQ1Fexpr *>

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
