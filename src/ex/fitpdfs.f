C     ****h* gepard/fitpdfs.f
C  FILE DESCRIPTION
C    calculation of x * PDF(x) for a given parameters
C    Produces data files for Fig. 14. of  KMPKS06b paper
C
C    $Id$
C     *******


C     ****p* fitpdfs.f/PDF
C  NAME
C     FITPDFS  --  Ploting PDFs resulting from fitted GPDs
C  SYNOPSIS

      PROGRAM FITPDFS

C  DESCRIPTION
C    Produces data files for Fig. 14 of  KMPKS06b
C  OUTPUT
C       fitpdfs[0-1].dat  --  2 files (one for each panel of Figure
C                            fitpdfs) with 4 sets of point coordinates
C                            corresponding to 4 lines on a graph
C                         (Alekhin's PDFs come in another file)
C  IDENTIFIERS
C          NPOINTS -- number of points on each line
C              XS  -- array(NPOINTS), values of x-coordinates 
C           POINTS -- array(2, 4, NPOINTS), holds results for
C                 two panels (Fig 14 a, b),and 4 lines are:
C                    LO, MS-bar NLO, CS-bar NLO, and NNLO
C
C                P -- approximation order N^{P}LO P=0,1,2
C
C  CHILDREN
C      READPAR, INIT, XSPACE
C  SOURCE
C
      IMPLICIT NONE
      INTEGER PT, NPOINTS, LN
      DOUBLE PRECISION HX(2)
      DOUBLE PRECISION LOGXI, LOGXISTART, LOGXIEND, LOGXISTEP, MP
      PARAMETER ( NPOINTS = 20 )
      DOUBLE PRECISION POINTS(0:1, 4, NPOINTS)
      DOUBLE PRECISION XIS(NPOINTS)
      PARAMETER ( LOGXISTART = -4.0d0, LOGXIEND = -2.0d0,
     &       LOGXISTEP = (LOGXIEND - LOGXISTART) / (NPOINTS - 1) )
*     Proton mass:
      PARAMETER (MP = 0.938272d0 )

      INCLUDE '../header.f'

      PROCESS = 'DVCS'
      FFTYPE = 'SINGLET'

      CALL READPAR
    
      SCHEME = 'CSBAR'
      CALL INIT

      INCLUDE 'ansatz.f'

* ----  "Common params, but different than in ansatz.f" ----------------

*   Del M^2 is 0
      PAR(15)   =     0.0
      PAR(25)   =     0.0
*   Removing valence quarks
      PAR(31)   =     0.0d0
      PAR(41)   =     0.0d0

*     Files that will hold results

      OPEN (UNIT = 10, FILE = "fitpdfs0.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 11, FILE = "fitpdfs1.dat", STATUS = "UNKNOWN")

      DO 5 LN = 0, 1
  5         WRITE (10 + LN, *) '# Output of fitpdfs.f. See prolog of 
     & that program'

      LOGXI = LOGXISTART
      DO 20 PT = 1, NPOINTS
      XI = 10**LOGXI
      XIS(PT) = XI

*     Looping over four lines
      DO 10 LN =1, 4

        IF (LN .EQ. 1) THEN
          P = 0
          PAR(11) =  0.157
          PAR(12) =  1.17
          PAR(14) =  0.228
          PAR(21) =  0.527
          PAR(22) =  1.25
          PAR(24) =  0.263
        ELSE IF (LN .EQ. 2) THEN
          SCHEME = 'MSBAR'
          P=1
          PAR(11) =  0.172
          PAR(12) =  1.14
          PAR(14) =  1.93
          PAR(21) =  0.472
          PAR(22) =  1.08
          PAR(24) =  4.45
        ELSE IF (LN .EQ. 3) THEN
          P=1
          PAR(11) =  0.167
          PAR(12) =  1.14
          PAR(14) =  1.34
          PAR(21) =  0.535
          PAR(22) =  1.09
          PAR(24) =  1.59
        ELSE
          P=2
          PAR(11) =  0.167
          PAR(12) =  1.14
          PAR(14) =  1.17
          PAR(21) =  0.571
          PAR(22) =  1.07
          PAR(24) =  1.39
        END IF

*     Calculating PDF needed for present line and point ...

      CALL XSPACE ( HX, XI, 0.0d0)

*     ... and saving them to array

      POINTS(0, LN, PT) = HX(1)
      POINTS(1, LN, PT) = HX(2)

 10   CONTINUE

 20   LOGXI = LOGXI + LOGXISTEP

*     Printing all the results from arrays to files

      DO 40 LN = 1, 4
      DO 30 PT = 1, NPOINTS
         WRITE (UNIT=10,FMT=998) XIS(PT), POINTS(0,LN,PT)
         WRITE (UNIT=11,FMT=998) XIS(PT), POINTS(1,LN,PT)
 30   CONTINUE
*     Empty line is put between the sets
      WRITE (UNIT=10, FMT=999)
      WRITE (UNIT=11, FMT=999)
 40   CONTINUE

998   FORMAT (F12.7,5X,F12.7)
999   FORMAT (1X)

      STOP
      END
C     ****
