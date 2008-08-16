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
      DOUBLE PRECISION AUX
      DOUBLE PRECISION XIA(21), REF(21)
      DOUBLE PRECISION PARSIGMA, SIGMA
#include "header.F"


      PROCESS = 'DVCS'
      FFTYPE = 'SINGLET'
      ANSATZ = 'FITNG'

      DATA XIA / 1.D-7,  5.D-7,  1.D-6,  5.D-6,  1.D-5,  5.D-5,
     ,           1.D-4,  5.D-4,  1.D-3,  5.D-3,  1.D-2,  5.D-2,
     ,           1.D-1,  2.D-1,  3.D-1,  4.D-1,  5.D-1,  6.D-1,
     ,           7.D-1,  8.D-1,  9.D-1 / 

*   Referent data obtained with ACC = 6, SPEED = 1
*   C = 0.5, PHI = 1.9 and integration up to Im(J)=250

      DATA REF
     &/ 0.1854666381799963E+03, 0.1044769683142155E+03,
     &  0.8134685795104900E+02, 0.4519037829947333E+02,
     &  0.3498995720825778E+02, 0.1922766065702439E+02,
     &  0.1483832568476377E+02, 0.8128902820774689E+01,
     &  0.6277182222944339E+01, 0.3428757155565274E+01,
     &  0.2607753291117352E+01, 0.1177576737039321E+01,
     &  0.6923991912747306E+00, 0.3083028778065646E+00,
     &  0.1530298461442887E+00, 0.7917463326306280E-01,
     &  0.4157713979249601E-01, 0.2172365817604435E-01,
     &  0.1099830785906399E-01, 0.5110478432107896E-02,
     &  0.1837500616132106E-02/

      OPEN (UNIT = 11, FILE = 'acc.dat', STATUS = 'UNKNOWN')

      PAR(21) = 0.4d0
      PAR(11) = 2.0d0/3.0d0 - PAR(21)
      PAR(12) = 1.1d0
      PAR(13) = 0.25d0
      PAR(14) = 1.1d0
      PAR(15) = 0.1d0
      PAR(16) = 3.0d0
      PAR(18) = 1.1d0
      PAR(19) = 0.9d0

      PAR(22) = 1.2d0
      PAR(23) = 0.25d0
      PAR(24) = 1.2d0
      PAR(14) = 1.2d0
      PAR(15) = 0.1d0
      PAR(16) = 1.0d0
      PAR(18) = 0.9d0
      PAR(19) = 1.1d0

      Q02 = 4.0d0

      CALL READPAR

      Q2 = 12.0D0
      NQS = 1
      QS(1) = Q2

       

      DO 30 SPEED = 1, 4
      WRITE (11, *) '# SPEED = ', SPEED 
      DO 20 K = 1, 21
        XI = XIA(K)
        CALL INIT
        CALL EVOLC(1)
        AUX = SIGMA ()
        WRITE (11, 901) XI, ABS((AUX-REF(K))/REF(K))
*       WRITE (11, 903) AUX
 20   CONTINUE
      WRITE (11, *) 
 30   CONTINUE

 901  FORMAT (1X, E7.1, 8X, E20.14)
 902  FORMAT (1X, E7.1, 8X, E20.14)
 903  FORMAT (1X, E22.16)
      STOP
      END
C     ****
