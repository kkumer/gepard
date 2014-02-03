C     ****h* gepard/getmbgpd.F
C  FILE DESCRIPTION
C    Initialization of GPDs
C
C    $Id: fit.F 63 2006-12-21 14:05:07Z kkumer $
C     *******

C     ****s* getmbgpd.F/GETMBGPD
C  NAME
C     GETMBGPD  --  calculate GPD values on the all MB contour points
C            
C  SYNOPSIS

      SUBROUTINE GETMBGPD

C  PARENTS
C    FCN, CFFF, XSPACE
C  CHILDREN
C    HJ, EJ
C  SOURCE
C
      IMPLICIT NONE
      INTEGER K
#include "header.F"
      DOUBLE COMPLEX J, FCM(4)


c$OMP parallel do PRIVATE(J, FCM)
      DO 30  K = 1, NPTS
        J = N(1, K) - 1
        CALL HJ(J, FCM)
        GPDH(K, 1) = FCM(1)
        GPDH(K, 2) = FCM(2)
        GPDH(K, 3) = FCM(3)
        GPDH(K, 4) = FCM(4)
        IF ( (ANSATZ .EQ. 'EFL') .OR. (ANSATZ .EQ. 'EFLEXP') ) THEN
*           we model E independently, so call EJ to overwrite FCM
            CALL EJ(J, FCM)
        END IF
*           otherwise, we use E \propto H 
        GPDE(K, 1) = FCM(1)
        GPDE(K, 2) = FCM(2)
        GPDE(K, 3) = FCM(3)
        GPDE(K, 4) = FCM(4)
 30   CONTINUE
c$OMP end parallel do

      RETURN
      END
C     ****
