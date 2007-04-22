C     ****h* gepard/initgpd.f
C  FILE DESCRIPTION
C    Initialization of GPDs
C
C    $Id: fit.F 63 2006-12-21 14:05:07Z kkumer $
C     *******

C     ****s* initgpd.f/INIGPD
C  NAME
C     INITGPD  --  initialize GPD values
C            
C  DESCRIPTION
C  SYNOPSIS
C     SUBROUTINE INITC
C
C  INPUTS 
C                NPAR  --  Number of variable parameters
C                   A  --  Values of all parameters
C               IFLAG  --  If .EQ. 3 then print results out, see minuit doc
C  OUTPUT
C                   F  --  Value of chi-square
C  IDENTIFIERS
C     FITPAR  --  Equal to A
C  CHISQPART  --  chi-square produced by single data set
C
C  CHILDREN
C      PROCDATA, PGPLOT routines
C
C  PARENTS
C      MINUIT, FIT (via MINUIT)
C
C  SOURCE
C


      SUBROUTINE INITGPD

      IMPLICIT NONE
      INTEGER K, EFFACC, NBNDMAX
      PARAMETER (NBNDMAX = 1)
      DOUBLE PRECISION BND(NBNDMAX+1), MTSEXP(4)
      INCLUDE 'header.f'
      DOUBLE COMPLEX J, FCM(2)

*     [0] MT=0 for forward (DIS) case
*     [1] MT's occuring in experimental data
      NMTSEXP = 4
      DATA MTSEXP / 0.1d0, 0.3d0, 0.5d0, 0.8d0 /

*     [2] MT's occuring as abscissae of Gaussian integration of PARSIGMA(MT)

*     Integration is done on the interval:
      DATA BND /0.0D0, 1.0D0/ 

*     Integrand PARSIGMA(MT) is nicely-behaved function and we
*     integrate it by equidistant Gauss-Legendre integration.
*     Effective numerical accuracy used for this integration is:
      EFFACC = ACC - SPEED + 1
*     which means that resulting number of integration points is:
      NMTS = 2**EFFACC

*     Calculating abscissas and weights
      CALL INTEGRAF(EFFACC, 1, BND, NBNDMAX, MTS, MTWG)

C     Initialize grids with values of GPD's for all MT's
      DO 20  K = 1, NPTS
        J = N(1,K) - 1

*       [1] forward:
        DEL2 = 0.0d0
        CALL HJ(J, FCM)
        HGRID(0, K, 1) = FCM(1)
        HGRID(0, K, 2) = FCM(2)

*       [2] Gaussian:
        DO 10 MTIND = 1, NMTS
          DEL2 = - MTS(MTIND)
          CALL HJ(J, FCM)
          HGRID(MTIND, K, 1) = FCM(1)
          HGRID(MTIND, K, 2) = FCM(2)
   10   CONTINUE

*       [3] PARSIGMA experimental data:
        DO 20 MTIND = NMTS+1, NMTS+NMTSEXP
          DEL2 = - MTSEXP(MTIND-NMTS)
          CALL HJ(J, FCM)
          HGRID(MTIND, K, 1) = FCM(1)
          HGRID(MTIND, K, 2) = FCM(2)
 20   CONTINUE

      RETURN
      END
C     ****
