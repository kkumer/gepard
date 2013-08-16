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
C    OUTPUT MCU(K) K=1,2,3 for quark, pure singlet quark and gluon contrib
C  SYNOPSIS

      SUBROUTINE CDVEMF (J, K, MCU)

!f2py complex(8) intent(in) :: J, K
!f2py complex(8) intent(out) :: MCU

      IMPLICIT NONE
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

      PJ = 1
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

<* MCMUFQ1Fexpr *>

<* MCMURQ1Fexpr *>

<* MCBETQ1Fexpr *>

<* MCFQ1Fexpr *>

<* MCGAQ1Fexpr *>

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

<* MCGBIAQ1Fexpr *>

      MCQ = MCQ + MCGBIAQ1F 

<* MCGBIBQ1Fexpr *>

      MCQ = MCQ + MCGBIBQ1F 

<* MCGBIIQ1Fexpr *>

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

<* MCBMUFQ1Fexpr *>

<* MCBAQ1Fexpr *>

<* MCBBQ1Fexpr *>

*     PRINT *, ' MCBMUFQ1F = ',  MCBMUFQ1F
*     PRINT *, ' MCBAQ1F = ',  MCBAQ1F
*     PRINT *, ' MCBBQ1F = ',  MCBBQ1F

*     PRINT *, ' NQBNORM = ',  NQBNORM
*     PRINT *, ' unnorm NLO QUARK = ',  (MCBMUFQ1F+MCBAQ1F+MCBBQ1F)
*     PRINT *, ' NLO QUARK = ',  NQBNORM * (MCBMUFQ1F+MCBAQ1F+MCBBQ1F)

C   ... gluon part

<* MCG1Fexpr *>

*     PRINT *, ' NGNORM = ',  NGNORM
*     PRINT *, ' unnorm G = ',  MCG1F

      MCU(1) = MCQ
      MCU(2) = (MCBMUFQ1F+MCBAQ1F+MCBBQ1F)
      MCU(3) = MCG1F

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
