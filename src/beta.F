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
C     NFL = 3...6  (parameters NFMIN, NFMAX) quark flavours.
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
      INTEGER NFL
      DOUBLE PRECISION B00, B01, B10, B11
#include "header.F"


*   Colour factors
* ..The full LO and NLO coefficients 

       B00 =  11./3.D0 * CA
       B01 =  -4./3.D0 * TF
       B10 =  34./3.D0 * CA**2
       B11 = -20./3.D0 * CA*TF - 4.* CF*TF

* ..Flavour-number loop and output to the array

       DO 1 NFL = NFMIN, NFMAX

       BETA0(NFL) = - B00 - B01 * NFL
       BETA1(NFL) = - B10 - B11 * NFL

       BETA2(NFL) = - 1428.50 + 279.611 * NFL - 6.01852 * NFL**2
       BETA3(NFL) = - 29243.0 + 6946.30 * NFL - 405.089 * NFL**2 
     &             - 1.49931 * NFL**3

  1    CONTINUE
       RETURN
       END
C     *******
