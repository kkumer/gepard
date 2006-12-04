C     ****h* gepard/radcorr.f
C  FILE DESCRIPTION
C    calculation of relative radiative corrections to {\cal H} singlet
C    DVCS form factor
C    Produces data files for Fig. 1. of  KMPKS06 Letter
C
C    $Id$
C     *******


C     ****p* radcorr.f/RADCORR
C  NAME
C     RADCORR  --  Main program
C  DESCRIPTION
C    calculation of relative radiative corrections to {\cal H} singlet
C    DVCS form factor
C    Produces data files for Fig. 1. of  KMPKS06 Letter
C  OUTPUT
C       FIG1[A-D].dat  --  4 files with six sets of point coordinates
C                          (\xi, K_\lambda) or (\xi, K_arg),
C                          corresponding to six lines on a graph
C  IDENTIFIERS
C          NPOINTS -- number of points on each line
C       LOGXISTART -- log(\xi) where line begins
C         LOGXIEND -- log(\xi) where line ends
C              XIS -- array(NPOINTS), values of XI
C           POINTS -- array(4, 6, NPOINTS), holds results for
C                 four panels (Fig 1, in order  a, b, c and d) and
C                 six K_{\lambda,arg} lines depicted on these panels
C                 Order is:  1. NNLO, CSBAR, HARD 
C                            2.  NLO, CSBAR, HARD
C                            3.  NLO, MSBAR, HARD
C                            4. NNLO, CSBAR, SOFT
C                            5.  NLO, CSBAR, SOFT
C                            6.  NLO, MSBAR, SOFT
C
C                P -- approximation order N^{P}LO P=0,1,2
C
C  CHILDREN
C      READPAR, INIT, CFFF, DCARG
C  SOURCE
C

      PROGRAM RADCORR

      IMPLICIT NONE
      INTEGER PT, NPOINTS, LN, NDEL
      DOUBLE PRECISION LOGXI, LOGXISTART, LOGXIEND, LOGXISTEP
      DOUBLE PRECISION DCARG
      PARAMETER ( NPOINTS = 2 )
      DOUBLE PRECISION POINTS(4, 6, NPOINTS)
      DOUBLE PRECISION XIS(NPOINTS)
      PARAMETER ( LOGXISTART = -5.0d0, LOGXIEND = -1.30103d0,
     &       LOGXISTEP = (LOGXIEND - LOGXISTART) / (NPOINTS - 1)  )
      INCLUDE 'header.f'


      CALL READPAR

*     Files that will hold results

      OPEN (UNIT = 11, FILE = "FIG1A.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 12, FILE = "FIG1B.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 13, FILE = "FIG1C.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 14, FILE = "FIG1D.dat", STATUS = "UNKNOWN")

      DO 5 NDEL = 1, 4
  5         WRITE (10 + NDEL, *) '# Output of radcorr.f. See prolog of 
     & that program'

*     Scales 

      Q2 = 2.5d0
      PAR(1) = 2.5d0
      PAR(2) = 0.05d0
      PAR(3) = 2.5d0

*     Looping over two different walues of \Delta^2

      DO 40 NDEL = 0, 1

      DEL2 = 0.0d0
      IF (NDEL .EQ. 1) DEL2 = -0.5d0

      LOGXI = LOGXISTART
      DO 20 PT = 1, NPOINTS
      XI = 10**LOGXI
      XIS(PT) = XI

*     Looping over six lines for each panel

      DO 10 LN = 1, 6

*     Recognizing parameters for each line

      IF ((LN .EQ. 1) .OR. (LN .EQ. 4)) THEN
            P = 2
      ELSE
            P = 1
      END IF
      IF ((LN .EQ. 3) .OR. (LN .EQ. 6)) THEN
            SCHEME = 'MSBAR'
      ELSE
            SCHEME = 'CSBAR'
      END IF
      IF (LN .LE. 3) THEN
            ANSATZ = 'HARD'
      ELSE
            ANSATZ = 'SOFT'
      END IF

      CALL INIT

*     Calculating two CFFs needed for present line and point ...

      CALL CFFF 
      P = P - 1
      CALL CFFF 

*     ... and saving them to array

      POINTS(1+NDEL, LN, PT) = LOG(ABS( CFF(P+1) )) / LOG(ABS( CFF(P) ))
      POINTS(3+NDEL, LN, PT) = DCARG( CFF(P+1) ) / DCARG( CFF(P) )

 10   CONTINUE

 20   LOGXI = LOGXI + LOGXISTEP

*     Printing all the results from arrays to files

      DO 40 LN = 1, 6
      DO 30 PT = 1, NPOINTS
         WRITE (UNIT=11+NDEL,FMT=998) XIS(PT), POINTS(1+NDEL, LN, PT)
         WRITE (UNIT=13+NDEL,FMT=998) XIS(PT), POINTS(3+NDEL, LN, PT)
 30   CONTINUE
*     Empty line is put between the sets
      WRITE (UNIT=11+NDEL, FMT=999)
      WRITE (UNIT=13+NDEL, FMT=999)
 40   CONTINUE

998   FORMAT (F12.7,5X,F12.7)
999   FORMAT (1X)

      STOP
      END
C     ****
