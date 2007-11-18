C     ****h* gepard/getmbgpd.f
C  FILE DESCRIPTION
C    Initialization of GPDs
C
C    $Id: fit.F 63 2006-12-21 14:05:07Z kkumer $
C     *******

C     ****s* initgpd.f/GETMBGPD
C  NAME
C     GETMBGPD  --  calculate GPD values on the all MB contour points
C            
C  SYNOPSIS

      SUBROUTINE GETMBGPD

C  SOURCE
C
      IMPLICIT NONE
      INTEGER K
      INCLUDE 'header.f'
      DOUBLE COMPLEX J, FCM(2)


C     Initialize grid with values of GPD's for forward case
      DO 30  K = 1, NPTS
        J = N(K) - 1
        CALL HJ(J, FCM)
        MBGPD(K, 1) = FCM(1)
        MBGPD(K, 2) = FCM(2)
 30   CONTINUE

      RETURN
      END
C     ****
