C     ****h* gepard/aux.f
C  FILE DESCRIPTION
C    Program for various temporary calculations and checks
C
C    $Id: test.f 27 2006-07-21 21:36:05Z kuk05260 $
C     *******

C     ****p* aux.f/AUX
C  NAME
C    AUX  --  Auxilliary program
C            
C  SOURCE
C

      PROGRAM AUX

      IMPLICIT NONE
      INTEGER K
      DOUBLE PRECISION W2, AUX
      DOUBLE PRECISION XIA(11), REF(11)
      DOUBLE PRECISION PARSIGMA, SIGMA
      INCLUDE 'header.f'

      DATA XIA /1.D-7, 1.D-6, 1.D-5, 1.D-4, 1.D-3,
     &  1.D-2, 1.D-1, 3.D-1, 5.D-1, 7.D-1, 9.D-1 /


      OPEN (UNIT = 11, FILE = 'aux.dat', STATUS = 'UNKNOWN')

      CALL READPAR

      PAR(21) = 0.4d0
      PAR(11) = 2.0d0/3.0d0 - PAR(21)
      PAR(12) = 1.1d0
      PAR(13) = 0.25d0
      PAR(14) = 1.1d0
      PAR(22) = 1.2d0
      PAR(23) = 0.25d0
      PAR(24) = 1.2d0

      PAR(1) = 1.0d0
      PAR(2) = 0.05d0
      PAR(3) = 2.5d0


      Q2 = 10.0D0

      DO 30 C = -0.9d0, 0.9d0, 0.4d0
      WRITE (11, *) '# C = ', C
      CALL INIT
      DO 20 K = 1, 11
        XI = XIA(K)
        CALL CFFF
        WRITE (11, 901) XI, ABS(CFF(P))
 20   CONTINUE
      WRITE (11, *) 
 30   CONTINUE

 901  FORMAT (1X, E7.1, 8X, E20.9)
      STOP
      END
C     ****
