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
      INTEGER SPEED, ACC, P, NF
      INTEGER K
      DOUBLE PRECISION NG, NSEA, MG, MSEA, ALPHA0G, ALPHA0SEA
      DOUBLE PRECISION FITPAR(10)
      DOUBLE PRECISION W2, XI, DEL2, Q2, Q02, AUX
      DOUBLE PRECISION F2(0:2), XIA(11), REF(11)
      DOUBLE PRECISION PARSIGMA, SIGMA

      COMMON / PARINT /  SPEED, ACC, P, NF

      COMMON / FITPARAMS  /  FITPAR
      COMMON / KINEMATICS /  XI, DEL2, Q2, Q02
      COMMON / F2P        /  F2

      DATA XIA /1.D-7, 1.D-6, 1.D-5, 1.D-4, 1.D-3,
     &  1.D-2, 1.D-1, 3.D-1, 5.D-1, 7.D-1, 9.D-1 /

*   Referent data obtained with ACC = 6, and integration
*   extended to Y = 80 (FIXME: one should also use more Gaussian
*   points for t-integration)

      DATA REF
     &/ 0.342836286E+04, 0.202172267E+04, 0.119790900E+04,
     &  0.715911446E+03, 0.432593578E+03, 0.252709543E+03,
     &  0.904676567E+02, 0.229166796E+02, 0.656041179E+01,
     &  0.176642456E+01, 0.294458260E+00 /

      OPEN (UNIT = 11, FILE = 'acc.dat', STATUS = 'UNKNOWN')

      NG = 0.4
      NSEA = 2./3. - NG
      MG = 1.2
      MSEA = 1.1
      ALPHA0G = 1.2
      ALPHA0SEA = 1.1
      Q02 = 1.0

      FITPAR(1) = NG 
      FITPAR(2) = NSEA 
      FITPAR(3) = MG 
      FITPAR(4) = MSEA 
      FITPAR(5) = ALPHA0G 
      FITPAR(6) = ALPHA0SEA 
      FITPAR(7) = Q02 

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
