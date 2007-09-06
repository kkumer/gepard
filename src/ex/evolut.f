C     ****h* gepard/evolut.f
C  FILE DESCRIPTION
C    evolution of {\cal H} non-singlet
C    DVCS form factor
C    Produces data files for Fig. 9 of  KMPKS06b paper
C
C    $Id$
C     *******


C     ****p* evolut.f/EVOLUT
C  NAME
C     EVOLUT  --  Program producing data for Figure evolut
C  SYNOPSIS
      
      PROGRAM EVOLUT

C  DESCRIPTION
C    evolution of NLO corrections to {\cal H} singlet
C    DVCS form factor in CSBAR and MSBAR scheme
C    Produces data files for Fig. 9 of  KMPKS06b
C  OUTPUT
C       evolut[0-3].dat  --  4 files (one for each panel of Figure
C                            evolut) with 4 sets of point coordinates
C                            (\xi, \Delta K) or (\xi, \Delta\varphi),
C                            corresponding to 4 lines on a graph
C  IDENTIFIERS
C          NPOINTS -- number of points on each line
C       LOGXISTART -- log(\xi) where line begins
C         LOGXIEND -- log(\xi) where line ends
C              XIS -- array(NPOINTS), values of XI
C           POINTS -- array(0:3, 6, NPOINTS), holds results for
C                 four panels (Fig 7, in order  a, b, c and d) and
C                 6 K_{\lambda,arg} lines depicted on these panels
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
C  SOURCE
C
      IMPLICIT NONE
      INTEGER PT, NPOINTS, LN, NDEL
      DOUBLE PRECISION LOGXI, LOGXISTART, LOGXIEND, LOGXISTEP, MP
      DOUBLE PRECISION DCARG, Q2MEM
      DOUBLE COMPLEX CFFQ, CFF0
      PARAMETER ( NPOINTS = 40 )
      DOUBLE PRECISION POINTS(0:3, 6, NPOINTS)
      DOUBLE PRECISION XIS(NPOINTS)
      PARAMETER ( LOGXISTART = -5.0d0, LOGXIEND = -0.30103d0,
     &       LOGXISTEP = (LOGXIEND - LOGXISTART) / (NPOINTS - 1) )
*     Proton mass:
      PARAMETER (MP = 0.938272d0 )

      INCLUDE '../header.f'

      PROCESS = 'DVCS'
      FFTYPE = 'SINGLET'

      CALL READPAR
    
      DEL2 = -0.25d0
      NQS = 1

      INCLUDE 'ansatz.f'
      ANSATZ = 'FIT'

*     Files that will hold results

      OPEN (UNIT = 10, FILE = "evolut0.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 11, FILE = "evolut1.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 12, FILE = "evolut2.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 13, FILE = "evolut3.dat", STATUS = "UNKNOWN")

      DO 5 NDEL = 0, 3
  5         WRITE (10 + NDEL, *) '# Output of evolut.f. See prolog of 
     & that program'



*     Looping over two different walues of Q^2

      DO 40 NDEL = 0, 1

      Q2 = 5.0d0
      IF (NDEL .EQ. 1) Q2 = 25.0d0


      LOGXI = LOGXISTART
      DO 20 PT = 1, NPOINTS
      XI = 10**LOGXI
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
            SCHEME = 'MSBND'
      END IF
      IF ((LN .EQ. 1) .OR. (LN .EQ. 3) .OR. (LN .EQ. 5)) THEN
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
      QS(1) = Q2
      CALL EVOLC(1)
      CALL GETMBGPD
      CALL CFFF 
      CFFQ = CFF(P)
      Q2MEM = Q2
      Q2 = Q02

      CALL INIT
      QS(1) = Q2
      CALL EVOLC(1)
      CALL GETMBGPD
      CALL CFFF 
      CFF0 = CFF(P)
      Q2 = Q2MEM

*     ... and saving them to array

      POINTS(NDEL, LN, PT) = (ABS(CFFQ) / ABS(CFF0) -
     &        1.0d0) * 100.0d0
      POINTS(2+NDEL, LN, PT) = DCARG(CFFQ / CFF0)

 10   CONTINUE

 20   LOGXI = LOGXI + LOGXISTEP

*     Printing all the results from arrays to files

      DO 40 LN = 1, 6
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
