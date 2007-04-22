C     ****h* gepard/slope.f
C  FILE DESCRIPTION
C    calculation of slope B for a given parameters
C    Produces data files for Fig. 15. of  KMPKS06b paper
C
C    $Id$
C     *******


C     ****p* slope.f/SLOPE
C  NAME
C     SLOPE  --  Ploting slope B resulting from fitted GPDs
C  DESCRIPTION
C    Produces data files for Fig. 14 of  KMPKS06b
C  OUTPUT
C       slope[0-1].dat  --  2 files (one for each panel of Figure
C                            slope) with 4 sets of point coordinates
C                            corresponding to 4 lines on a graph
C                         (Alekhin's PDFs come in another file)
C  IDENTIFIERS
C          NPOINTS -- number of points on each line
C              XS  -- array(NPOINTS), values of x-coordinates 
C           POINTS -- array(2, 6, NPOINTS), holds results for
C                 two panels (Fig 14 a, b),and 6 lines are:
C                  low-band CS NLO, hi-band CS NLO,
C                  LO, MS-bar NLO, CS-bar NLO, and NNLO
C
C                P -- approximation order N^{P}LO P=0,1,2
C
C  CHILDREN
C      READPAR, INIT, LNHX, DERIV
C  SOURCE
C

      PROGRAM SLOPE

      IMPLICIT NONE
      INTEGER PT, NPOINTS, LN, NEVALS
      DOUBLE PRECISION HX(2)
      DOUBLE PRECISION LOGXI, LOGXISTART, LOGXIEND, LOGXISTEP, MP
      DOUBLE PRECISION LNHXQ, LNHXG, T, H, DRV, ERR
      EXTERNAL LNHXQ, LNHXG
      PARAMETER ( NEVALS = 5, H = 0.41d0 )
      PARAMETER ( NPOINTS = 20 )
      DOUBLE PRECISION POINTS(0:1, 6, NPOINTS)
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

      OPEN (UNIT = 10, FILE = "slope0.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 11, FILE = "slope1.dat", STATUS = "UNKNOWN")

      DO 5 LN = 0, 1
  5         WRITE (10 + LN, *) '# Output of fitpdfs.f. See prolog of 
     & that program'

      LOGXI = LOGXISTART
      DO 20 PT = 1, NPOINTS
      XI = 10**LOGXI
      XIS(PT) = XI

*     Looping over six lines
      DO 10 LN =1, 6

        IF (LN .EQ. 3) THEN
          P = 0
          PAR(11) =  0.157
          PAR(12) =  1.17
          PAR(14) =  0.228
          PAR(21) =  0.527
          PAR(22) =  1.25
          PAR(24) =  0.263
        ELSE IF (LN .EQ. 4) THEN
          SCHEME = 'MSBAR'
          P=1
          PAR(11) =  0.172
          PAR(12) =  1.14
          PAR(14) =  1.93
          PAR(21) =  0.472
          PAR(22) =  1.08
          PAR(24) =  4.45
        ELSE IF (LN .EQ. 5) THEN
          P=1
          PAR(11) =  0.167
          PAR(12) =  1.14
          PAR(14) =  1.34
          PAR(21) =  0.535
          PAR(22) =  1.09
          PAR(24) =  1.59
        ELSE IF (LN .EQ. 1) THEN
          P=1
          PAR(11) =  0.167
          PAR(12) =  1.14
          PAR(14) =  1.23
          PAR(21) =  0.535
          PAR(22) =  1.09
          PAR(24) =  1.19
        ELSE IF (LN .EQ. 2) THEN
          P=1
          PAR(11) =  0.167
          PAR(12) =  1.14
          PAR(14) =  1.45
          PAR(21) =  0.535
          PAR(22) =  1.09
          PAR(24) =  1.99
        ELSE
          P=2
          PAR(11) =  0.167
          PAR(12) =  1.14
          PAR(14) =  1.17
          PAR(21) =  0.571
          PAR(22) =  1.07
          PAR(24) =  1.39
        END IF

*     Calculating slopes needed for present line and point ...
*     ... and saving them to array

      T = 0.0d0

      CALL DERIV(LNHXQ, T, H, NEVALS, DRV, ERR)
!      WRITE (30, *) XI, DRV, ERR
      POINTS(0, LN, PT) = DRV

      T = 0.0d0

      CALL DERIV(LNHXG, T, H, NEVALS, DRV, ERR)
!      WRITE (31, *) XI, DRV, ERR
      POINTS(1, LN, PT) = DRV



 10   CONTINUE

 20   LOGXI = LOGXI + LOGXISTEP

*     Printing all the results from arrays to files

*     First the error band

      DO 25 PT = 1, NPOINTS
         WRITE (UNIT=10,FMT=998) XIS(PT), POINTS(0,1,PT)
         WRITE (UNIT=11,FMT=998) XIS(PT), POINTS(1,1,PT)
 25   CONTINUE
      DO 26 PT = NPOINTS, 1, -1
         WRITE (UNIT=10,FMT=998) XIS(PT), POINTS(0,2,PT)
         WRITE (UNIT=11,FMT=998) XIS(PT), POINTS(1,2,PT)
 26   CONTINUE
*     Empty line is put between the sets
      WRITE (UNIT=10, FMT=999)
      WRITE (UNIT=11, FMT=999)

      DO 40 LN = 3, 6
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



      DOUBLE PRECISION FUNCTION LNHXQ(T)

      IMPLICIT NONE
      DOUBLE PRECISION HX(2), T
      INCLUDE '../header.f'

      CALL XSPACE ( HX, XI, T)

      LNHXQ = LOG( HX(1) )

      RETURN
      END

      DOUBLE PRECISION FUNCTION LNHXG(T)

      IMPLICIT NONE
      DOUBLE PRECISION HX(2), T
      INCLUDE '../header.f'

      CALL XSPACE ( HX, XI, T)

      LNHXG = LOG( HX(2) )

      RETURN
      END
