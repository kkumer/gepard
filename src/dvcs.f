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
C
C  IDENTIFIERS
C
C                   W2 -- (gamma*-proton) invariant mass squared
C
C  CHILDREN
C      CFFF
C
C  BUGS
C       dimension of FITPAR is hardcoded 
C
C  SOURCE
C

      DOUBLE PRECISION FUNCTION PARSIGMA (MT)

      IMPLICIT NONE
      DOUBLE PRECISION MT
      INTEGER P
      DOUBLE PRECISION NGIN, FITPAR(10)
      DOUBLE PRECISION XI, DEL2, Q2, Q02, W2
      DOUBLE PRECISION NG, NSEA, MG, MSEA, ALPHA0G, ALPHA0SEA
      DOUBLE COMPLEX CFF(0:2)
      CHARACTER SCHEME*5, ANSATZ*6

*     Input common-blocks 

      COMMON / FITPARAMS  /  FITPAR

*     Output common-blocks 

      COMMON / KINEMATICS /  XI, DEL2, Q2, Q02
      COMMON / APPROX     /  P
      COMMON / LABELS     /  SCHEME, ANSATZ
      COMMON / GPD        /  NG, NSEA, MG, MSEA, ALPHA0G, ALPHA0SEA
      COMMON / CFF        /  CFF


*     We do fitting at NLO in CSBAR scheme. It is here that this
*     is specified and written in common blocks.

      ANSATZ = 'FIT'
      P = 1

*     Scales and kinematics

      Q02 = FITPAR(7)

      DEL2 = - MT
      W2 = (Q2 - XI * Q2) / (2.0d0 * XI)

*     GPD parameters

      NG = FITPAR(1)
      NSEA = FITPAR(2)
      MG = FITPAR(3)
      MSEA = FITPAR(4)
      ALPHA0G = FITPAR(5)
      ALPHA0SEA = FITPAR(6)

*     Calculating CFF

*     260.5633976788416 = 4 Pi alpha^2 * (GeV^-2 -> nbarn)
      
      CALL CFFF 
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
C    and integrates over t numerically, using simple Bode 
C    formula. Seems to be accurate to 5 percent.
C
C  INPUTS 
C              FITPAR  --  Array with fit parameters, see FIT.CMD
C                          for description
C  OUTPUT
C               SIGMA  --  total cross section sigma
C
C  CHILDREN
C     PARSIGMA   -- integrand, partial cross section  d sigma / d |t|
C  SOURCE
C

      DOUBLE PRECISION FUNCTION SIGMA ()

      IMPLICIT NONE
      DOUBLE PRECISION FITPAR(10)
*   Integration parameters:
      INTEGER IER, IWORK, LAST, LENW, LIMIT, NEVAL
      DOUBLE PRECISION A, TCUT, LAM, ABSERR, EPSABS, EPSREL, WORK
      DOUBLE PRECISION RES
      PARAMETER (LIMIT = 400, LENW = 4 * LIMIT)
      DIMENSION IWORK(LIMIT), WORK(LENW)
      DOUBLE PRECISION PARSIGMA

*   Input common-blocks 

      COMMON / FITPARAMS  /  FITPAR

*   Integration using Bode formula

      SIGMA = (14.0d0 * PARSIGMA(0.0d0) + 64.0d0 * PARSIGMA(0.25d0) +
     &  24.0d0 * PARSIGMA(0.5d0) + 64.0d0 * PARSIGMA(75.0d0) +
     &  14.0d0 * PARSIGMA(1.0d0)) / 180d0

      RETURN
      END
C     ***
