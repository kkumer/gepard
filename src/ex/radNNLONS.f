C     ****h* gepard/radNNLONS.f
C  FILE DESCRIPTION
C    calculation of relative radiative corrections to {\cal H} non-singlet
C    DVCS form factor
C    Produces data files for Fig. 11. of  KMPKS06b paper
C
C    $Id$
C     *******


C     ****p* radNNLONS.f/RADNNLONS
C  NAME
C     RADNNLONS  --  Program producing data for Figure radNNLONS
C  DESCRIPTION
C    calculation of relative radiative NNLONS and NLO corrections to {\cal H} non-singlet
C    DVCS form factor in CSBAR scheme
C    Produces data files for Fig. 11 of  KMPKS06b
C  OUTPUT
C       radNNLONS[0-3].dat  --  4 files (one for each panel of Figure
C                            radNNLONS) with 4 sets of point coordinates
C                            (\xi, \delta K) or (\xi, \delta\varphi),
C                            corresponding to 4 lines on a graph
C  IDENTIFIERS
C          NPOINTS -- number of points on each line
C       LOGXISTART -- log(\xi) where line begins
C         LOGXIEND -- log(\xi) where line ends
C              XIS -- array(NPOINTS), values of XI
C           POINTS -- array(0:3, 4, NPOINTS), holds results for
C                 four panels (Fig 7, in order  a, b, c and d) and
C                 4 K_{\lambda,arg} lines depicted on these panels
C                 Order is:  1.  NNLO, CSBAR, HARD 
C                            2.  NNLO, CSBAR, SOFT
C                            3.   NLO, CSBAR, HARD
C                            4.   NLO, CSBAR, SOFT
C
C                P -- approximation order N^{P}LO P=0,1,2
C
C  CHILDREN
C      READPAR, INIT, CFFF, DCARG
C  SOURCE
C

      PROGRAM RADNNLONS

      IMPLICIT NONE
      INTEGER PT, NPOINTS, LN, NDEL
      DOUBLE PRECISION XISTART, XIEND, XISTEP, MP
      DOUBLE PRECISION DCARG
      PARAMETER ( NPOINTS = 60 )
      DOUBLE PRECISION POINTS(0:3, 4, NPOINTS)
      DOUBLE PRECISION XIS(NPOINTS)
      PARAMETER ( XISTART = 0.01d0, XIEND = 0.5d0,
     &       XISTEP = (XIEND - XISTART) / (NPOINTS - 1) )
*     Proton mass:
      PARAMETER (MP = 0.938272d0 )

      INCLUDE '../header.f'


      CALL READPAR
    
      NF = 4

      ANSATZ = 'NSFIT'
      INCLUDE 'ansatz.f'


*     Files that will hold results

      OPEN (UNIT = 10, FILE = "radNNLONS0.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 11, FILE = "radNNLONS1.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 12, FILE = "radNNLONS2.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 13, FILE = "radNNLONS3.dat", STATUS = "UNKNOWN")

      DO 5 NDEL = 0, 3
  5         WRITE (10 + NDEL, *) '# Output of radNNLONS.f. See prolog of 
     & that program'


      DEL2 = -0.25d0
      SCHEME = 'CSBAR'

*     Looping over two different walues of Q^2

      DO 40 NDEL = 0, 1

      Q2 = 2.5d0
      IF (NDEL .EQ. 1) Q2 = 10.0d0

      XI = XISTART
      DO 20 PT = 1, NPOINTS
      XIS(PT) = XI

*     Looping over four lines for each panel

      DO 10 LN = 1, 4

*     Recognizing parameters for each line

      IF (LN .LE. 2) THEN
            P = 2
      ELSE
            P = 1
      END IF
      IF ((LN .EQ. 1) .OR. (LN .EQ. 3)) THEN
!          'VALENCE + SEA'
          PAR(11) = 4.0d0 / 15.0d0
      ELSE
!          'VALENCE ONLY'
          PAR(11) = 0.0d0
      END IF


*     Calculating two CFFs needed for present line and point ...

      CALL INIT
      CALL CFFF 
      P = P - 1
      CALL CFFF 
      P = P + 1

*     ... and saving them to array

      POINTS(NDEL, LN, PT) = (ABS(CFF(P)) / ABS(CFF(P-1)) -
     &        1.0d0) * 100.0d0
      POINTS(2+NDEL, LN, PT) = DCARG( CFF(P) / CFF(P-1) )

 10   CONTINUE

 20   XI = XI + XISTEP

*     Printing all the results from arrays to files

      DO 40 LN = 1, 4
      DO 30 PT = 1, NPOINTS
         WRITE (UNIT=10+NDEL,FMT=998) XIS(PT), POINTS(NDEL,LN,PT)
         WRITE (UNIT=12+NDEL,FMT=998) XIS(PT), POINTS(2+NDEL,LN,PT)
 30   CONTINUE
*     Empty line is put between the sets
      WRITE (UNIT=10+NDEL, FMT=999)
      WRITE (UNIT=12+NDEL, FMT=999)
 40   CONTINUE

998   FORMAT (F12.7,5X,F12.7)
999   FORMAT (1X)

      STOP
      END
C     ****
