C     ****h* gepard/dvcs.F
C  FILE DESCRIPTION
C    calculation of partial and total DVCS cross sections
C
C    $Id$
C     *******

C     ****f* dvcs.F/PARSIGMA
C  NAME
C    PARSIGMA  --  calculation of partial DVCS cross section
C                  d sigma / d t 
C  DESCRIPTION
C    Calls CFFF for calculation of Compton form factor,
C    and multiplies by required kinematical factors. It is intended
C    to be called by chi-square minimization program for
C    fitting the experimental data.
C  SYNOPSIS

      DOUBLE PRECISION FUNCTION PARSIGMA ()

C  PARENTS
C      SIGMA, PROCDATA
C  CHILDREN
C      CFFF
C  SOURCE
C


      IMPLICIT NONE
#include "header.F"

*     Scales and kinematics

      W2 = (Q2 - XI * Q2) / (2.0d0 * XI)

*     Calculating CFF
      
      DEL2 = - MTS(MTIND)
      CALL CFFF 

*     260.5633976788416 = 4 Pi alpha^2 * (GeV^-2 -> nbarn)
*     3.52142 = 4 MP^2

      PARSIGMA = 260.5633976788416d0 * W2 * (
     &  ABS(CFF(P))**2 - DEL2/3.52142d0*ABS(CFFE(P))**2 ) / (
     &  (W2 + Q2) * (2.0d0 * W2 + Q2)**2 )


      RETURN
      END
C     ****


C     ****f* dvcs.F/SIGMA
C  NAME
C       SIGMA  --  calculation of total DVCS cross section
C  DESCRIPTION
C    Calls PARSIGMA for calculation of partial cross section
C    and integrates over t numerically, using Gauss method
C  SYNOPSIS

      DOUBLE PRECISION FUNCTION SIGMA ()

C  PARENTS
C     PROCDATA
C  CHILDREN
C     PARSIGMA   -- integrand, partial cross section  d sigma / d |t|
C  SOURCE
C

      IMPLICIT NONE
      DOUBLE PRECISION PARSIGMA
#include "header.F"

*     Integrating
      SIGMA = 0.0d0
      DO 10 MTIND = 1, NMTS
          SIGMA = SIGMA + MTWG(MTIND) * PARSIGMA()
 10   CONTINUE

      RETURN
      END
C     ***
