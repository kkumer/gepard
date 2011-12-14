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
C    HJ
C  SOURCE
C
      IMPLICIT NONE
      INTEGER K
#include "header.F"
      DOUBLE COMPLEX J, FCM(2)


c$OMP parallel do PRIVATE(J, FCM)
      DO 30  K = 1, NPTS
        J = N(1, K) - 1
        CALL HJ(J, FCM)
        MBGPD(K, 1) = FCM(1)
        MBGPD(K, 2) = FCM(2)
        IF (PAR(111) .GT. 0.0d0) THEN
*           we model E independently
            CALL EJ(J, FCM)
            MBGPD(K, 3) = FCM(1)
            MBGPD(K, 4) = FCM(2)
        ELSE
*           use E \propto H
            MBGPD(K, 3) = FCM(1)
            MBGPD(K, 4) = FCM(2)
        END IF
 30   CONTINUE
c$OMP end parallel do

      RETURN
      END
C     ****
