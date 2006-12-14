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
      DOUBLE PRECISION XIA(11), REF(11)
      DOUBLE PRECISION PARSIGMA, SIGMA
      INCLUDE 'header.f'


      PROCESS = 'DVCS'
      FFTYPE = 'SINGLET'

      DATA XIA /1.D-7, 1.D-6, 1.D-5, 1.D-4, 1.D-3,
     &  1.D-2, 1.D-1, 3.D-1, 5.D-1, 7.D-1, 9.D-1 /

*   Referent data obtained with ACC = 6, and integration
*   extended to Y = 80 (FIXME: one should also use more Gaussian
*   points for t-integration)

!      DATA REF
!     &/ 0.342836286E+04, 0.202172267E+04, 0.119790900E+04,
!     &  0.715911446E+03, 0.432593578E+03, 0.252709543E+03,
!     &  0.904676567E+02, 0.229166796E+02, 0.656041179E+01,
!     &  0.176642456E+01, 0.294458260E+00 /

      DATA REF
     &/ 0.169301759E+03, 0.998381005E+02, 0.591559717E+02,
     &  0.353536377E+02, 0.213626391E+02, 0.124794806E+02,
     &  0.446753786E+01, 0.113168772E+01, 0.323970911E+00,
     &  0.872308312E-01, 0.145411468E-01 /


      OPEN (UNIT = 11, FILE = 'acc.dat', STATUS = 'UNKNOWN')

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

      CALL READPAR

      Q2 = 4.0D0
      ACC = 3

      DO 30 SPEED = 1, 4
      WRITE (11, *) '# SPEED = ', SPEED 
      CALL INIT
      DO 20 K = 1, 11
        XI = XIA(K)
        AUX = SIGMA ()
        WRITE (11, 901) XI, ABS((AUX-REF(K))/REF(K))
 20   CONTINUE
      WRITE (11, *) 
 30   CONTINUE

 901  FORMAT (1X, E7.1, 8X, E20.9)
      STOP
      END
C     ****
