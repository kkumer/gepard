C     ****h* gepard/radNNLO.f
C  FILE DESCRIPTION
C    calculation of relative radiative corrections to {\cal H} singlet
C    DVCS form factor
C    Produces data files for Fig. 7. of  KMPKS06b paper
C
C    $Id$
C     *******


C     ****p* radNNLO.f/RADNNLO
C  NAME
C     RADNNLO  --  Program producing data for Figure radNNLO
C  DESCRIPTION
C    calculation of relative radiative NNLO and NLO corrections to {\cal H} singlet
C    DVCS form factor in CSBAR scheme
C    Produces data files for Fig. 9 of  KMPKS06b
C  SYNOPSIS

      PROGRAM RADNNLO

C  OUTPUT
C       radNNLO[0-3].dat  --  4 files (one for each panel of Figure
C                            radNNLO) with 4 sets of point coordinates
C                            (\xi, K_\lambda) or (\xi, \delta\varphi),
C                            corresponding to 4 lines on a graph
C  IDENTIFIERS
C          LNPOINTS -- number of points on each line
C       LOGXISTART -- log(\xi) where line begins
C         LOGXIEND -- log(\xi) where line ends
C              XIS -- array(LNPOINTS), values of XI
C           POINTS -- array(0:3, 4, LNPOINTS), holds results for
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
      IMPLICIT NONE
      INTEGER PT, LNPOINTS, LN, NDEL
      DOUBLE PRECISION LOGXI, LOGXISTART, LOGXIEND, LOGXISTEP, MP
      DOUBLE PRECISION DCARG
      PARAMETER ( LNPOINTS = 40 )
      DOUBLE PRECISION POINTS(0:3, 4, LNPOINTS)
      DOUBLE PRECISION XIS(LNPOINTS)
      PARAMETER ( LOGXISTART = -5.0d0, LOGXIEND = -0.30103d0,
     &       LOGXISTEP = (LOGXIEND - LOGXISTART) / (LNPOINTS - 1) )
*     Proton mass:
      PARAMETER (MP = 0.938272d0 )

#include "../header.F"


      PROCESS = 'DVCS'
      FFTYPE = 'SINGLET'

      CALL READPAR

#include "ansatz.F"
      ANSATZ = 'FITBP'
*     ANSATZ = 'NSFIT'
      NQS = 1


*     Files that will hold results

      OPEN (UNIT = 10, FILE = "radNNLO0.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 11, FILE = "radNNLO1.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 12, FILE = "radNNLO2.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 13, FILE = "radNNLO3.dat", STATUS = "UNKNOWN")

      DO 5 NDEL = 0, 3
  5         WRITE (10 + NDEL, *) '# Output of radNNLO.f. See prolog of 
     & that program'


      DEL2 = -0.25d0
      SCHEME = 'CSBAR'

*     Looping over two different walues of Q^2

      DO 40 NDEL = 0, 1

      Q2 = 2.5d0
      IF (NDEL .EQ. 1) Q2 = 100.0d0
      NQS = 1
      QS(1) = Q2

      LOGXI = LOGXISTART
      DO 20 PT = 1, LNPOINTS
      XI = 10**LOGXI
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
!          'HARD'
            PAR(21) = 0.4d0
            PAR(22) = PAR(12) + 0.05d0
      ELSE
!           'SOFT'
            PAR(21) = 0.3d0
            PAR(22) = PAR(12) - 0.2d0
      END IF
      PAR(11) = (2.0d0/3.0d0) - PAR(21)


*     Calculating two CFFs needed for present line and point ...

      CALL INIT
      DEL2 = -0.25d0
      CALL EVOLC(1,1)
      CALL CFFF 
      P = P - 1
      CALL EVOLC(1,1)
      CALL CFFF 
      P = P + 1

*     ... and saving them to array

*     IF ((XI .LE. 1.1e-5) .AND. (Q2 .LE. 10.) .AND.
*    & ((LN .EQ. 1) .OR. (LN .EQ. 3))) THEN
*         WRITE(*,*) XI
*         WRITE(*,*) CFF
*         WRITE(*,*) (ABS(CFF(P)) / ABS(CFF(P-1)) -
*    &        1.0d0) * 100.0d0
*     END IF
      POINTS(NDEL, LN, PT) = (ABS(CFF(P)) / ABS(CFF(P-1)) -
     &        1.0d0) * 100.0d0
      POINTS(2+NDEL, LN, PT) = DCARG(CFF(P)) - DCARG(CFF(P-1))

 10   CONTINUE

 20   LOGXI = LOGXI + LOGXISTEP

*     Printing all the results from arrays to files

      DO 40 LN = 1, 4
      DO 30 PT = 1, LNPOINTS
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
