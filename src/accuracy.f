C     ****h* gepard/accuracy.f
C  FILE DESCRIPTION
C    Program for measuring accuracy for various SPEED
C    settings
C
C    $Id: test.f 27 2006-07-21 21:36:05Z kuk05260 $
C     *******

C     ****p* accuracy.f/ACCURACY
C  NAME
C    ACCURACY  --  Testing gepard routines
C            
C  DESCRIPTION
C    Gives relative errors for range of XI, and
C    for various SPEED settings.
C
C  CHILDREN
C    READPAR, INIT, SIGMA
C  SOURCE
C

      PROGRAM ACCURACY

      IMPLICIT NONE
      INTEGER K
      DOUBLE PRECISION W2, AUX
      DOUBLE PRECISION XIA(21), REF(21)
      DOUBLE PRECISION PARSIGMA, SIGMA
      INCLUDE 'header.f'


      PROCESS = 'DVCS'
      FFTYPE = 'SINGLET'

      DATA XIA / 1.D-7,  5.D-7,  1.D-6,  5.D-6,  1.D-5,  5.D-5,
     ,           1.D-4,  5.D-4,  1.D-3,  5.D-3,  1.D-2,  5.D-2,
     ,           1.D-1,  2.D-1,  3.D-1,  4.D-1,  5.D-1,  6.D-1,
     ,           7.D-1,  8.D-1,  9.D-1 / 

*   Referent data obtained with ACC = 6, SPEED = 1
*   C = 0.5, PHI = 1.9

      DATA REF
     &/ 0.16930189D+03, 0.11699769D+03, 0.99838143D+02, 0.69202573D+02,
     &  0.59155984D+02, 0.41230115D+02, 0.35353641D+02, 0.24844213D+02,
     &  0.21362640D+02, 0.14871283D+02, 0.12479480D+02, 0.69631732D+01,
     &  0.44675379D+01, 0.21702898D+01, 0.11316877D+01, 0.60441352D+00,
     &  0.32397091D+00, 0.17133662D+00, 0.87230831D-01, 0.40557338D-01,
     &  0.14541146D-01 /


      OPEN (UNIT = 11, FILE = 'acc.dat', STATUS = 'UNKNOWN')

      PAR(21) = 0.4d0
      PAR(11) = 2.0d0/3.0d0 - PAR(21)
      PAR(12) = 1.1d0
      PAR(13) = 0.25d0
      PAR(14) = 1.1d0
      PAR(22) = 1.2d0
      PAR(23) = 0.25d0
      PAR(24) = 1.2d0

      Q02 = 1.0d0

      CALL READPAR

      Q2 = 4.0D0
      NQS = 1
      QS(1) = Q2

       

      DO 30 SPEED = 1, 4
      WRITE (11, *) '# SPEED = ', SPEED 
      DO 20 K = 1, 21
        XI = XIA(K)
        CALL INIT
        CALL INITGPD(1)
        CALL EVOLC(1)
        AUX = SIGMA ()
        WRITE (11, 901) XI, ABS((AUX-REF(K))/REF(K))
 20   CONTINUE
      WRITE (11, *) 
 30   CONTINUE

 901  FORMAT (1X, E7.1, 8X, E20.9)
 902  FORMAT (1X, E7.1, 8X, E20.9)
      STOP
      END
C     ****
