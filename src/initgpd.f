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
C  SOURCE
C


      SUBROUTINE INITGPD ( MT )

      IMPLICIT NONE
      DOUBLE PRECISION MT
      INTEGER K
      INCLUDE 'header.f'
      DOUBLE COMPLEX J, FCM(2)


C     Initialize grid with values of GPD's for all MT's
      DO 30  K = 1, NPTS
        J = N(1,K) - 1

        DO 10 MTIND = 0, NMTS+NMTSEXP
          DEL2 = - MTS(MTIND)
          CALL HJ(J, FCM)
          HGRID(MTIND, K, 1) = FCM(1)
          HGRID(MTIND, K, 2) = FCM(2)
 10     CONTINUE

*       Value from subroutine argument is put at the end of array
        DEL2 = - MT
        CALL HJ(J, FCM)
        MTIND = MTIND + 1
        HGRID(MTIND, K, 1) = FCM(1)
        HGRID(MTIND, K, 2) = FCM(2)

 30   CONTINUE

      RETURN
      END
C     ****
