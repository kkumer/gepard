C     ****h* gepard/radNLONS.f
C  FILE DESCRIPTION
C    calculation of relative radiative corrections to {\cal H} non-singlet
C    DVCS form factor
C    Produces data files for Fig. 7. of  KMPKS06b paper
C
C    $Id$
C     *******


C     ****p* radNLONS.f/RADNLONS
C  NAME
C     RADNLONS  --  Program producing data for Figure radNLONS
C  DESCRIPTION
C    calculation of relative radiative NLO corrections to {\cal H} non-singlet
C    DVCS form factor in MSBAR and CSBAR schemes
C    Produces data files for Fig. 7 of  KMPKS06b
C  SYNOPSIS

      PROGRAM RADNLONS

C  OUTPUT
C       radNLONS[0-3].dat  --  4 files (one for each panel of Figure
C                            radNLONS) with 4 sets of point coordinates
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
C                 Order is:  1.  NLO, CSBAR, HARD 
C                            2.  NLO, CSBAR, SOFT
C                            3.  NLO, MSBAR, HARD
C                            4.  NLO, MSBAR, SOFT
C
C                P -- approximation order N^{P}LO P=0,1,2
C
C  CHILDREN
C      READPAR, INIT, CFFF, DCARG
C  SOURCE
C
      IMPLICIT NONE
      INTEGER PT, LNPOINTS, LN, NDEL
      DOUBLE PRECISION XISTART, XIEND, XISTEP, MP
      DOUBLE PRECISION DCARG
      PARAMETER ( LNPOINTS = 40 )
      DOUBLE PRECISION POINTS(0:3, 4, LNPOINTS)
      DOUBLE PRECISION XIS(LNPOINTS)
      PARAMETER ( XISTART = 0.01d0, XIEND = 0.5d0,
     &       XISTEP = (XIEND - XISTART) / (LNPOINTS - 1) )
*     Proton mass:
      PARAMETER (MP = 0.938272d0 )

#include "../header.F"

      PROCESS = 'DVCS'
      FFTYPE = 'NONSINGLET'

      CALL READPAR

#include "ansatz.F"
      ANSATZ = 'NSFIT'
      NQS = 1

*     Files that will hold results

      OPEN (UNIT = 10, FILE = "radNLONS0.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 11, FILE = "radNLONS1.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 12, FILE = "radNLONS2.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 13, FILE = "radNLONS3.dat", STATUS = "UNKNOWN")

      DO 5 NDEL = 0, 3
  5         WRITE (10 + NDEL, *) '# Output of radNLONS.f. See prolog of 
     & that program'

*     Scales 

      Q2 = 2.5d0

*     Looping over two different walues of \Delta^2

      DO 40 NDEL = 0, 1

      DEL2 = 0.0d0
      IF (NDEL .EQ. 1) DEL2 = -1.0d0

      XI = XISTART
      DO 20 PT = 1, LNPOINTS
      XIS(PT) = XI

*     Looping over four lines for each panel

      DO 10 LN = 1, 4

*     Recognizing parameters for each line

      IF (LN .LE. 2) THEN
            SCHEME = 'CSBAR'
      ELSE
            SCHEME = 'MSBDI'
      END IF
      IF ((LN .EQ. 1) .OR. (LN .EQ. 3)) THEN
!          'HARD'
            PAR(11) = 4.0d0 / 15.0d0
      ELSE
!           'SOFT'
            PAR(11) = 0.0d0
      END IF


*     Calculating two CFFs needed for present line and point ...

      P = 1
      CALL INIT
      QS(1) = Q2
      CALL EVOLC(1,1)
      CALL CFFF 
      
      P = 0
      CALL INIT
      QS(1) = Q2
      CALL EVOLC(1,1)
      CALL CFFF 

*     ... and saving them to array

      POINTS(NDEL, LN, PT) = (ABS(CFF(1)) / ABS(CFF(0)) -
     &        1.0d0) * 100.0d0
      POINTS(2+NDEL, LN, PT) = DCARG(CFF(1) / CFF(0))

 10   CONTINUE

 20   XI = XI + XISTEP

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
