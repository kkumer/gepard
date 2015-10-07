C     ****h* gepard/evolutNS.f
C  FILE DESCRIPTION
C    evolution of {\cal H} non-singlet
C    DVCS form factor
C    Produces data files for Fig. 8 of  KMPKS06b paper
C
C    $Id$
C     *******


C     ****p* evolutNS.f/EVOLUTNS
C  NAME
C     EVOLUTNS  --  Program producing data for Figure evolutNS
C  DESCRIPTION
C    evolution of NLO corrections to {\cal H} singlet
C    DVCS form factor in CSBAR and MSBAR scheme
C    Produces data files for Fig. 8 of  KMPKS06b
C  OUTPUT
C       evolutNS[0-3].dat  --  4 files (one for each panel of Figure
C                            evolutNS) with 4 sets of point coordinates
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
C                 Order is:  1.   LO,        HARD 
C                            2.   LO,        SOFT
C                            3.  NLO, CSBAR, HARD 
C                            4.  NLO, CSBAR, SOFT
C                            5.  NLO, MSBAR, HARD
C                            6.  NLO, MSBAR, SOFT
C
C                P -- approximation order N^{P}LO P=0,1,2
C
C  CHILDREN
C      READPAR, INIT, CFFF, DCARG
C  SYNOPSIS

      PROGRAM EVOLUTNS

C  SOURCE
C
      IMPLICIT NONE
      INTEGER PT, LNPOINTS, LN, NDEL
      DOUBLE PRECISION XISTART, XIEND, XISTEP, MP
      DOUBLE PRECISION DCARG, Q2MEM
      DOUBLE COMPLEX CFFQ, CFF0
      PARAMETER ( LNPOINTS = 40 )
      DOUBLE PRECISION POINTS(0:3, 6, LNPOINTS)
      DOUBLE PRECISION XIS(LNPOINTS)
      PARAMETER ( XISTART = 0.01d0, XIEND = 0.5d0,
     &       XISTEP = (XIEND - XISTART) / (LNPOINTS - 1) )
*     Proton mass:
      PARAMETER (MP = 0.938272d0 )

#include "../header.F"

      PID = 1

      CALL READPAR
    
      DEL2 = -0.25d0

#include "ansatz.F"
      ANSATZ = 'NSFIT'

*     Files that will hold results

      OPEN (UNIT = 10, FILE = "evolutNS0.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 11, FILE = "evolutNS1.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 12, FILE = "evolutNS2.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 13, FILE = "evolutNS3.dat", STATUS = "UNKNOWN")

      DO 5 NDEL = 0, 3
  5         WRITE (10 + NDEL, *) '# Output of evolutNS.f. See prolog of 
     & that program'



*     Looping over two different walues of Q^2

      DO 40 NDEL = 0, 1

      Q2 = 1.0d0
      IF (NDEL .EQ. 1) Q2 = 10.0d0


      XI = XISTART
      DO 20 PT = 1, LNPOINTS
      XIS(PT) = XI

*     Looping over six lines for each panel

      DO 10 LN = 1, 6

*     Recognizing parameters for each line

      IF (LN .LE. 2) THEN
            P = 0
            SCHEME = 'CSBAR'
      ELSEIF  (LN .LE. 4) THEN
            P = 1
            SCHEME = 'CSBAR'
      ELSE
            P = 1
            SCHEME = 'MSBDI'
      END IF
      IF ((LN .EQ. 1) .OR. (LN .EQ. 3) .OR. (LN .EQ. 5)) THEN
!          'HARD'
            PAR(11) = 4.0d0 / 15.0d0
      ELSE
!           'SOFT'
            PAR(11) = 0.0d0
      END IF



*     Calculating two CFFs needed for present line and point ...

      CALL INIT
      CALL CFFF 
      CFFQ = CFF(P)
      Q2MEM = Q2

      Q2 = Q02
      CALL INIT
      CALL CFFF 
      CFF0 = CFF(P)
      Q2 = Q2MEM

*     ... and saving them to array

      POINTS(NDEL, LN, PT) = (ABS(CFFQ) / ABS(CFF0) -
     &        1.0d0) * 100.0d0
      POINTS(2+NDEL, LN, PT) = DCARG(CFFQ / CFF0)

 10   CONTINUE

 20   XI = XI + XISTEP

*     Printing all the results from arrays to files

      DO 40 LN = 1, 6
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
