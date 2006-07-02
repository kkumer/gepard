C     ****h* gepard/plotsigma.f
C  FILE DESCRIPTION
C    Simple fast plotting with all parameters fixed
C
C    $Id$
C     *******

C     ****p* plotsigma.f/PLOTSIGMA
C  NAME
C    PLOTSIGMA   --  Plots sigma(Q2)
C            
C  DESCRIPTION
C             Specifies fixed values for parameters
C
C  OUTPUT
C          File SIGMA.PLT with plot points
C
C  CHILDREN
C      SIGMA
C
C  SOURCE
C


      PROGRAM PLOTSIGMA

      IMPLICIT NONE
      INTEGER N, NMAX, NPAR, IFLAG, I
      PARAMETER (NMAX=40, NPAR=7)
      DOUBLE PRECISION FITPAR(10)
      DOUBLE PRECISION X(NMAX), FITFN
      DOUBLE PRECISION W2, XI, DEL2, Q2, Q02
      DOUBLE PRECISION NG, NSEA, MG2, MSEA2, ALPHA0G, ALPHA0SEA
      DOUBLE PRECISION SIGMA, ISIGMA
      DATA (X(I), I=1,6) /3.0, 5.25, 8.75, 15.5, 25.0, 55.0/

      COMMON / KINEMATICS /  XI, DEL2, Q2, Q02
      COMMON / FITPARAMS  /  FITPAR

C     For Test/plotsigma.test
!      NG = 0.1
!      NSEA = 0.07
!      MG2 = 1.1
!      MSEA2 = 1.1
!      ALPHA0G = 1.0
!      ALPHA0SEA = 1.0
!      Q02 = 1.5

      NG = 0.10168
      NSEA = 0.052312
      MG2 = 1.0
      MSEA2 = 1.0
      ALPHA0G = 1.1
      ALPHA0SEA = 0.8
      Q02 = 1.0


      FITPAR(1) = NG 
      FITPAR(2) = NSEA 
      FITPAR(3) = MG2 
      FITPAR(4) = MSEA2 
      FITPAR(5) = ALPHA0G 
      FITPAR(6) = ALPHA0SEA 
      FITPAR(7) = Q02 


*      Ploting sigma(Q2) corresponding to:
*     DATASET 3  [H1,  Eur.Phys.J.C44:1-11,2005, hep-ex/0505061]
*       X = Q2,   Y = SIG   (TCUT = -1 GeV)
*       W = 82 GeV

      W2 = 82.0d0**2

      N = 1
      DO 10 I= 1, N
      Q2 = X(I)
      XI = Q2 / ( 2.0d0 * W2 + Q2)
      FITFN = SIGMA()
      WRITE (*, *) X(I), FITFN
 10   CONTINUE

      STOP
      END
C     ****
