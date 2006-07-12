C     ****h* gepard/beta.f
C  FILE DESCRIPTION
C    coefficients of beta function of QCD
C
C    $Id$
C  NOTES
C    This is modification of file betafct.f
C    from package `QCD-Pegasus'  by  Andreas Vogt.
C     *******


C     ****s* beta.f/BETAF
C  NAME
C     BETAF  --  coefficients of beta function of QCD
C  DESCRIPTION
C     The subroutine BETAF for the coefficients  BETA0...BETA3  of the 
C     beta function of QCD up to order alpha_s^5 (N^3LO), normalized by 
C
C          d a_s / d ln mu_r^2  =   BETA0 a_s^2 + BETA1 a_s^3 + ... 
C
C      with  a_s = alpha_s/(4*pi). 
C
C     NB: This normalization is different then Vogt's: BETAx -> -BETAx !!
C
C     The MSbar coefficients are written to the common-block  BETA  for 
C     NF = 3...6  (parameters NFMIN, NFMAX) quark flavours.
C
C      Beyond NLO the QCD colour factors are hard-wired in this routine,
C      and the numerical coefficients are truncated to six digits.
C  SYNOPSIS
C     SUBROUTINE BETAF
C  PARENTS
C     INIT
C  SOURCE
C

      SUBROUTINE BETAF

      IMPLICIT NONE
      INTEGER NFMIN, NFMAX, NF
      DOUBLE PRECISION CF, CA, TR
      DOUBLE PRECISION BETA0, BETA1, BETA2, BETA3
      DOUBLE PRECISION B00, B01, B10, B11
      PARAMETER (NFMIN = 3, NFMAX = 6)
      PARAMETER (CA = 3.D0, CF = 4./3.D0, TR = 0.5 D0)

*   Output common-block

       COMMON / BETABLK / BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX),
     &                    BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)

*   Colour factors
* ..The full LO and NLO coefficients 

       B00 =  11./3.D0 * CA
       B01 =  -4./3.D0 * TR
       B10 =  34./3.D0 * CA**2
       B11 = -20./3.D0 * CA*TR - 4.* CF*TR

* ..Flavour-number loop and output to the array

       DO 1 NF = NFMIN, NFMAX

       BETA0(NF) = - B00 - B01 * NF
       BETA1(NF) = - B10 - B11 * NF

       BETA2(NF) = - 1428.50 + 279.611 * NF - 6.01852 * NF**2
       BETA3(NF) = - 29243.0 + 6946.30 * NF - 405.089 * NF**2 
     &             - 1.49931 * NF**3

  1    CONTINUE
       RETURN
       END
C     *******
