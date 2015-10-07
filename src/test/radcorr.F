C     ****h* gepard/radcorr.F
C  FILE DESCRIPTION
C    calculation of relative radiative corrections to {\cal H} singlet
C    DVCS form factor
C    Produces data files for Fig. 1. of  KMPKS06 Letter
C
C    $Id$
C     *******


C     ****p* radcorr.F/RADCORR
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
C          NNPOINTS -- number of points on each line
C       LOGXISTART -- log(\xi) where line begins
C         LOGXIEND -- log(\xi) where line ends
C              XIS -- array(NNPOINTS), values of XI
C           POINTS -- array(4, 6, NNPOINTS), holds results for
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
      INTEGER PT, NNPOINTS, LN, NDEL
      DOUBLE PRECISION LOGXI, LOGXISTART, LOGXIEND, LOGXISTEP
      DOUBLE PRECISION DCARG
      PARAMETER ( NNPOINTS = 2 )
      DOUBLE PRECISION POINTS(4, 6, NNPOINTS)
      DOUBLE PRECISION XIS(NNPOINTS)
      PARAMETER ( LOGXISTART = -5.0d0, LOGXIEND = -1.30103d0,
     &       LOGXISTEP = (LOGXIEND - LOGXISTART) / (NNPOINTS - 1)  )
#include "header.F"

      PID = 1

      CALL READPAR

*     Files that will hold results

      OPEN (UNIT = 11, FILE = "FIG1A.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 12, FILE = "FIG1B.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 13, FILE = "FIG1C.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 14, FILE = "FIG1D.dat", STATUS = "UNKNOWN")

      DO 5 NDEL = 1, 4
  5         WRITE (10 + NDEL, *) '# Output of radcorr.F. See prolog of 
     & that program'

*     Scales 

      Q02 = 2.5d0
      Q2 = 2.5d0

*     Looping over two different walues of \Delta^2

      DO 40 NDEL = 0, 1

      LOGXI = LOGXISTART
      DO 20 PT = 1, NNPOINTS
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
*     At the time of writing this paper we didn't have non-diagonal
*     MSBAR evolution, thus:
            SCHEME = 'MSBDI'
      ELSE
            SCHEME = 'CSBAR'
      END IF
      IF (LN .LE. 3) THEN
            ANSATZ = 'HARD'
      ELSE
            ANSATZ = 'SOFT'
      END IF

      CALL INIT
      MTIND = 0
      IF (NDEL .EQ. 1) MTIND = NMTS + 3
      DEL2 = -MTS(MTIND)

*     Calculating two CFFs needed for present line and point ...

      CALL CFFF 
      P = P - 1
      CALL INIT
      CALL CFFF 

*     ... and saving them to array

      POINTS(1+NDEL, LN, PT) = LOG(ABS( CFF(P+1) )) / LOG(ABS( CFF(P) ))
      POINTS(3+NDEL, LN, PT) = DCARG( CFF(P+1) ) / DCARG( CFF(P) )

 10   CONTINUE

 20   LOGXI = LOGXI + LOGXISTEP

*     Printing all the results from arrays to files

      DO 40 LN = 1, 6
      DO 30 PT = 1, NNPOINTS
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
