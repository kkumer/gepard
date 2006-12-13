C     ****h* gepard/integra.f
C  FILE DESCRIPTION
C    Gaussian integration between array of points
C
C    $Id: dvcs.f 11 2006-07-12 07:50:28Z kuk05260 $
C     *******

C     ****s* integra.f/INTEGRA
C  NAME
C        INTEGRA  -- Gaussian integration between array of points
C  DESCRIPTION
C     Returns array of abscissas and weights for Gaussian integration
C     over some region which is subdivided by an array of points so
C     that points are more dense in the subregion where function changes
C     more violently. Subintegration intervals are defined as:
C     DIVISION(1) .. DIVISION(1 + SPEED) .. DIVISION(1 + 2*SPEED) ..
C        ... DIVISION(MAXDIV + 1), 
C     and on each interval 2**ACC Gauss-Legendre points are used.
C  SYNOPSIS
C     SUBROUTINE INTEGRA
C   INPUT
C      ACC  -- accuracy. each subintegration is done on 2**ACC
C              Gauss-Legendre points
C    SPEED  -- speed. Only each SPEEDth point from array DIVISION is
C              used for specification of subintegration intervals
C DIVISION  -- array with abscissae defining subintegration intervals
C NINTGMAX  -- how many subintervals for SPEED=1 (should be divisible by 4)
C              DIVISION should be of length NINTGMAX+1
C  OUTPUT
C        X  -- values of abscissae
C       WG  -- weights
C  CHILDREN
C      GAUSS
C  PARENTS
C      INIT
C  SOURCE
C

      SUBROUTINE INTEGRAF(ACC, SPEED, DIVISION, NINTGMAX, X, WG)
* 
      IMPLICIT NONE
      INTEGER ACC, SPEED, NINTGMAX
      INTEGER NGAUSSMAX, ACCMAX, NINTGMAXMAX, K, K2, K3
      DOUBLE PRECISION NGAUSS, SUMM, DIFF
      PARAMETER ( NGAUSSMAX = 64, ACCMAX = 6, NINTGMAXMAX = 12 )
      DOUBLE PRECISION DIVISION( NINTGMAXMAX+1 )
      DOUBLE PRECISION X( NINTGMAXMAX * NGAUSSMAX ) 
      DOUBLE PRECISION WG( NINTGMAXMAX * NGAUSSMAX )
      DOUBLE PRECISION ABSC( ACCMAX, NGAUSSMAX ) 
      DOUBLE PRECISION WGHT( ACCMAX, NGAUSSMAX )
      DOUBLE PRECISION DOWN( NINTGMAXMAX + 1 ), UP( NINTGMAXMAX )


      CALL GAUSS (ABSC, WGHT)

*   Abscissas X(K) and weights WG(K) for the whole contour. Note that
*   the factor (upper limit - lower limit)/2 is put into the weights

      K = 0
      DO 20 K2 = 1, NINTGMAX, SPEED
        SUMM = DIVISION(K2 + SPEED) + DIVISION(K2) 
        DIFF = DIVISION(K2 + SPEED) - DIVISION(K2) 
      DO 20 K3 = 1, 2**ACC
        K = K + 1
        X(K) = (DIFF * ABSC(ACC, K3) + SUMM) * 0.5D0
        WG(K) = 0.5D0 * DIFF * WGHT(ACC, K3) 
 20   CONTINUE


      RETURN
      END 
C     ***
