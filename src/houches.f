C     ****h* gepard/houches.f
C  FILE DESCRIPTION
C    Program for testing evolution
C
C    $Id: test.f 27 2006-07-21 21:36:05Z kuk05260 $
C     *******

C     ****p* houches.f/HOUCHES
C  NAME
C    HOUCHES  --  Testing gepard routines
C            
C  DESCRIPTION
C     Calculates evolution of PDFs (defined
C     as in Les Houches accord) by using special
C     'EVOLQ' and 'EVOLG' schemes where C = (1, 0) and
C      (0, 1), respecively, to all orders.
C
C  CHILDREN
C     INIT, PARWAVF, AS2PF
C  SOURCE
C

      PROGRAM HOUCHES

      IMPLICIT NONE
      INTEGER SPEED, ACC, P, NF
      DOUBLE PRECISION AS0, RF2, RR2
      CHARACTER SCHEME*5, ANSATZ*6
      DOUBLE PRECISION XI, DEL2, Q2, Q02
      INTEGER NPTS, NPTSMAX, K,  K1, K2
      DOUBLE COMPLEX EPH, J, FPW
      DOUBLE PRECISION PHI, C, F2IMAG, RES, XG(4)
      PARAMETER (NPTSMAX = 768)
      DOUBLE PRECISION Y(NPTSMAX), WG(NPTSMAX), PI, AS
      DOUBLE COMPLEX N(NPTSMAX)
      DOUBLE PRECISION LHXGINP(11), LHXGLO(11), LHXGNLO(11)
      DOUBLE PRECISION LHXGNNLO(11), XIA(11), REF(11)
      PARAMETER (PI = 3.14159265358979 D0)

      DATA XIA /1.D-7, 1.D-6, 1.D-5, 1.D-4, 1.D-3,
     &  1.D-2, 1.D-1, 3.D-1, 5.D-1, 7.D-1, 9.D-1 /


*  Data produced by QCD-Pegasus

      DATA LHXGLO
     &/ 1.3162D+03, 6.0008D+02, 2.5419D+02, 9.7371D+01, 3.2078D+01,
     &  8.0546D+00, 8.8766D-01, 8.2676D-02, 7.9240D-03, 3.7311D-04,
     &  1.0918D-06 /

      DATA LHXGNLO
     &/ 1.0636D+03, 5.1144D+02, 2.2761D+02, 9.1215D+01, 3.1275D+01,
     &  8.1038D+00, 9.0281D-01, 8.4047D-02, 8.1200D-03, 3.9204D-04,
     &  1.2428D-06 /

      DATA LHXGNNLO
     &/ 9.6982D+02, 4.8505D+02, 2.2190D+02, 9.0517D+01, 3.1321D+01,
     &  8.1365D+00, 9.0687D-01, 8.4334D-02, 8.1290D-03, 3.9034D-04,
     &  1.2191D-06 /

*  Double precision referent data produced by this very program. 
*  Agrees with Pegasus to all decimal places produced by him

      DATA REF
     &/ 0.106364828E+04, 0.511439128E+03, 0.227608411E+03,
     &  0.912151640E+02, 0.312752364E+02, 0.810377586E+01,
     &  0.902814476E+00, 0.840469129E-01, 0.812004733E-02,
     &  0.392036567E-03, 0.124277587E-05 /


*   Input common-blocks 

      COMMON / CONTOUR    /  NPTS
      COMMON / CPHI       /  C, PHI
      COMMON / POINTS     /  Y, WG
      COMMON / NPOINTS    /  N

*   Output common-blocks 

      COMMON / PARINT /  SPEED, ACC, P, NF
      COMMON / PARFLT /  AS0, RF2, RR2
      COMMON / PARCHR /  SCHEME, ANSATZ

      COMMON / KINEMATICS /  XI, DEL2, Q2, Q02

*   Parameters of Les Houches benchmark

      NF = 4
      RF2 = 1.0D0
      RR2 = 1.0D0
      ANSATZ = 'HOUCHE'
      Q02 = 2.0D0

*   Following artifical renorm. scheme has C = (0, 1) at
*   all orders, thus effectively putting F_2 -> Q_s^2 x g(x)

      SCHEME = 'EVOLG'

      OPEN (UNIT = 11, FILE = 'NLOACC.dat', STATUS = 'UNKNOWN') 
      OPEN (UNIT = 12, FILE = 'NNLOACC.dat', STATUS = 'UNKNOWN') 


      ACC = 3

      DO 210 SPEED = 1, 4
      WRITE (11, *) '# SPEED = ', SPEED
      WRITE (12, *) '# SPEED = ', SPEED
      WRITE (*,*) 'SPEED = ', SPEED

*  Gluon part - NLO evolution
      DO 35 K1 = 1, 11
        XI = XIA(K1)
        RES = 0.0D0
        P = 1
        Q2 = 1.0D4
        AS0 = 0.0525303832D0
        CALL INIT
        EPH = EXP ( COMPLEX(0.0d0, PHI) )
        DO 30 K = 1, NPTS
          J = N(K) - 1
          CALL PARWAVF (K, FPW, 'DIS')
          F2IMAG = IMAGPART(EPH * (1.0d0/XI)**(J-C) * FPW )
          RES = RES + WG(K)*F2IMAG
 30     CONTINUE
        XG(3) = (1.0d0/XI)**(C) * RES / 3.1415926535897932
!        WRITE (11, 803) XI, ABS((XG(3)-LHXGNLO(K1))/LHXGNLO(K1))
        WRITE (11, 803) XI, ABS((XG(3)-REF(K1))/REF(K1))
 35   CONTINUE
      CALL AS2PF (AS, Q2, AS0, 2.5D0, NF, P)  
      WRITE (*,904) Q2, AS * 2.0D0 * PI

*  Gluon part - NNLO evolution
      DO 45 K1 = 1, 11
        XI = XIA(K1)
        RES = 0.0D0
        P = 2
        Q2 = 1.0D4
        AS0 = 0.0524395701
        CALL INIT
        EPH = EXP ( COMPLEX(0.0d0, PHI) )
        DO 40 K = 1, NPTS
          J = N(K) - 1
          CALL PARWAVF (K, FPW, 'DIS')
          F2IMAG = IMAGPART(EPH * (1.0d0/XI)**(J-C) * FPW )
          RES = RES + WG(K)*F2IMAG
 40     CONTINUE
        XG(4) = (1.0d0/XI)**(C) * RES / 3.1415926535897932
        WRITE (12, 803) XI, ABS((XG(4)-LHXGNNLO(K1))/LHXGNNLO(K1))
 45   CONTINUE
      CALL AS2PF (AS, Q2, AS0, 2.5D0, NF, P)  
      WRITE (*,904) Q2, AS * 2.0D0 * PI

      WRITE (11, *)
      WRITE (12, *)
 210  CONTINUE

 803  FORMAT (1X, E7.1, 8X, E20.9)
 901  FORMAT (1X, 'x', 15X, 'x g(x)', 15X, 'rel. diff.', /
     &  1X, 61('-'))
 902  FORMAT (1X, E7.1, 8X, F9.4, 15X, E8.1)
 903  FORMAT (1X, E7.1, 8X, E10.5, 15X, E8.1)
 904  FORMAT (1X, 'alpha_strong(', E8.1, ') = ', F8.6)

      STOP
      END
C     ****
