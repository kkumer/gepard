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

      DOUBLE PRECISION FUNCTION PARSIGMA (MT)

      IMPLICIT NONE
      DOUBLE PRECISION MT
      INTEGER SPEED, ACC, P, NF
      CHARACTER SCHEME*5, ANSATZ*6
      DOUBLE PRECISION NGIN, FITPAR(10)
      DOUBLE PRECISION XI, DEL2, Q2, Q02, W2
      DOUBLE PRECISION NG, NSEA, MG, MSEA, ALPHA0G, ALPHA0SEA
      DOUBLE COMPLEX CFF(0:2)

*     Input common-blocks 

      COMMON / PARINT /  SPEED, ACC, P, NF
      COMMON / PARCHR /  SCHEME, ANSATZ
      
      COMMON / FITPARAMS  /  FITPAR

*     Output common-blocks 

      COMMON / KINEMATICS /  XI, DEL2, Q2, Q02
      COMMON / GPD        /  NG, NSEA, MG, MSEA, ALPHA0G, ALPHA0SEA
      COMMON / CFF        /  CFF


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
      INTEGER SPEED, ACC, P, NF
      INTEGER K3
      DOUBLE PRECISION FITPAR(10)
      DOUBLE PRECISION ABSCISSAS(8), WEIGHTS(8), ABS2, WGH2
      DOUBLE PRECISION PARSIGMA, ABS4(2), WGH4(2)

*   Input common-blocks 

      COMMON / PARINT /  SPEED, ACC, P, NF

      COMMON / FITPARAMS  /  FITPAR

*   Abscissas and weights for 8 point Gauss quadrature 
*   according to Abramowitz and Stegun Eq. (25.4.30)

       DATA ABSCISSAS
     &  /-0.96028 98564 97536 D0, -0.79666 64774 13627 D0, 
     &   -0.52553 24099 16329 D0, -0.18343 46424 95650 D0, 
     &    0.18343 46424 95650 D0,  0.52553 24099 16329 D0,
     &    0.79666 64774 13627 D0,  0.96028 98564 97536 D0/

       DATA WEIGHTS
     &  / 0.10122 85362 90376 D0,  0.22238 10344 53374 D0, 
     &    0.31370 66458 77887 D0,  0.36268 37833 78362 D0, 
     &    0.36268 37833 78362 D0,  0.31370 66458 77887 D0,
     &    0.22238 10344 53374 D0,  0.10122 85362 90376 D0/

*   Abscissas and weights for 4 point Gauss quadrature

       ABS4(1) = 0.86113 63115 94052 D0
       ABS4(2) = 0.33998 10435 84856 D0
       WGH4(1) = 0.34785 48451 37453 D0
       WGH4(2) = 0.65214 51548 62546 D0

*   Abscissas and weights for 2 point Gauss quadrature

      ABS2 = 0.57735 02691 89625 D0
      WGH2 = 1.00 D0

      IF ( SPEED .EQ. 1) THEN
        SIGMA = 0.D0
        DO 20 K3 = 1, 8
          SIGMA = SIGMA + WEIGHTS(K3) * 
     &            PARSIGMA( 0.5D0 * ABSCISSAS(K3) + 0.5D0 )
 20     CONTINUE
      ELSEIF ( SPEED .EQ. 2) THEN
        SIGMA = ( WGH4(1) * ( PARSIGMA(- 0.5D0 * ABS4(1) + 0.5D0) +
     &              PARSIGMA(0.5D0 * ABS4(1) + 0.5D0) ) +
     &            WGH4(2) * ( PARSIGMA(- 0.5D0 * ABS4(2) + 0.5D0) +
     &              PARSIGMA(0.5D0 * ABS4(2) + 0.5D0) ) )
      ELSE
        SIGMA = WGH2 * ( PARSIGMA(- 0.5D0 * ABS2 + 0.5D0) +
     &            PARSIGMA(0.5D0 * ABS2 + 0.5D0) )
      END IF

      SIGMA = 0.5D0 * SIGMA

      RETURN
      END
C     ***
