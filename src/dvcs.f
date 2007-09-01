C     ****h* gepard/dvcs.f
C  FILE DESCRIPTION
C    calculation of partial and total DVCS cross sections
C
C    $Id$
C     *******

C     ****f* dvcs.f/PARSIGMA
C  NAME
C    PARSIGMA  --  calculation of partial DVCS cross section
C                  d sigma / d t as a function of t, GPD
C                  parameters and Q0
C  DESCRIPTION
C    Calls CFFF for calculation of Compton form factor,
C    and multiplies by required kinematical factors. It is intened
C    to be called by chi-square minimization program for
C    fitting the experimental data.
C  SYNOPSIS
C     DOUBLE PRECISION FUNCTION PARSIGMA (MT)
C
C     DOUBLE PRECISION MT
C  INPUTS 
C                  MT  --  (- t) = (- DEL2), minus
C                          DVCS asymmetry parameter t=(P2-P1)^2
C              FITPAR  --  Array with fit parameters, see FIT.CMD
C                          for description (via common block)
C  OUTPUT
C            PARSIGMA  --  partial cross section d sigma / d t
C  IDENTIFIERS
C
C          W2 -- (gamma*-proton) invariant mass squared
C          XI -- DVCS scaling parameter
C        DEL2 -- DVCS asymmetry parameter (P2-P1)^2
C          Q2 -- photon virtuality squared
C         Q02 -- initial scale squared
C  PARENTS
C      SIGMA, PROCDATA
C  CHILDREN
C      CFFF
C  BUGS
C       dimension of FITPAR is hardcoded 
C  SOURCE
C

      DOUBLE PRECISION FUNCTION PARSIGMA ()

      IMPLICIT NONE
      DOUBLE PRECISION W2
      INCLUDE 'header.f'

*     Scales and kinematics

      W2 = (Q2 - XI * Q2) / (2.0d0 * XI)

*     Calculating CFF

*     260.5633976788416 = 4 Pi alpha^2 * (GeV^-2 -> nbarn)
      
      DEL2 = - MTS(MTIND)
      CALL CFFF 

*     Multiplying with charge factor

      IF (NF .EQ. 3) THEN
         CFF(P) = CFF(P) * 2.0D0 / 9.0D0
      ELSE IF (NF .EQ. 4) THEN
         CFF(P) = CFF(P) * 5.0D0 / 18.0D0
      ELSE
         CALL ERROR ('GeParD', '  CFF',
     &   'NF is not integer equal to 3 or 4!                          ',
     &   5, 3)
      END IF

      PARSIGMA = 260.5633976788416d0 * W2 * ABS(CFF(P))**2 / (
     &  (W2 + Q2) * (2.0d0 * W2 + Q2)**2 )


      RETURN
      END
C     ****


C     ****f* dvcs.f/SIGMA
C  NAME
C       SIGMA  --  calculation of total DVCS cross section
C                  sigma as a function of GPD
C                  parameters and Q0
C  SYNOPSIS
C     DOUBLE PRECISION FUNCTION SIGMA ()
C  DESCRIPTION
C    Calls PARSIGMA for calculation of partial cross section
C    and integrates over t numerically, using Gauss method
C
C  INPUTS 
C              FITPAR  --  Array with fit parameters, see FIT.CMD
C                          for description
C  OUTPUT
C               SIGMA  --  total cross section sigma
C
C  PARENTS
C     PROCDATA
C  CHILDREN
C     PARSIGMA   -- integrand, partial cross section  d sigma / d |t|
C  SOURCE
C

      DOUBLE PRECISION FUNCTION SIGMA ()

      IMPLICIT NONE
      DOUBLE PRECISION PARSIGMA
      INCLUDE 'header.f'

*     Integrating
      SIGMA = 0.0d0
      DO 10 MTIND = 1, NMTS
          SIGMA = SIGMA + MTWG(MTIND) * PARSIGMA()
 10   CONTINUE

      RETURN
      END
C     ***
