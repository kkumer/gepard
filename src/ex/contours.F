C     ****h* gepard/contours.f
C  FILE DESCRIPTION
C    calculation of parton density contours
C
C    $Id$
C     *******


C     ****p* contours.f/CONTOURS
C  NAME
C        CONTOURS  --  Ploting parton density contours resulting from fitted GPDs
C  DESCRIPTION
C  OUTPUT
C  IDENTIFIERS
C          LNPOINTS -- number of points on each line
C             DDS  -- array(LNPOINTS), values of Delta-coordinates 
C              XS  -- array(LNPOINTS), values of x-coordinates 
C           POINTS -- array(2, LNPOINTS, LNPOINTS), holds results
C                P -- approximation order N^{P}LO P=0,1,2
C
C  CHILDREN
C      READPAR, INIT, XSPACE
C  SOURCE
C

      PROGRAM CONTOURS

      IMPLICIT NONE
      INTEGER PTD, PTX, LNPOINTS, LN
      DOUBLE PRECISION HX(2)
      DOUBLE PRECISION DD, DDSTART, DDEND, DDSTEP, MP, T
      DOUBLE PRECISION X, LOGX, LOGXSTART, LOGXEND, LOGXSTEP
      PARAMETER ( LNPOINTS = 50 )
      DOUBLE PRECISION POINTS(2, LNPOINTS, LNPOINTS)
      DOUBLE PRECISION DDS(LNPOINTS), XS(LNPOINTS)
      PARAMETER ( DDSTART = 0.0d0, DDEND = 5.0d0,
     &       DDSTEP = (DDEND - DDSTART) / (LNPOINTS - 1) )
      PARAMETER ( LOGXSTART = -5.0d0, LOGXEND = -0.30103d0,
     &       LOGXSTEP = (LOGXEND - LOGXSTART) / (LNPOINTS - 1) )
*     Proton mass:
      PARAMETER (MP = 0.938272d0 )

#include "../header.F"

      PROCESS = 'DVCS'
      FFTYPE = 'SINGLET'

      CALL READPAR
    
      SCHEME = 'CSBAR'
      CALL INIT

#include "ansatz.F"

* ----  "Common params, but different than in ansatz.f" ----------------

*   Del M^2 is 0
      PAR(15)   =     0.0
      PAR(25)   =     0.0
*   Removing valence quarks
      PAR(31)   =     0.0d0
      PAR(41)   =     0.0d0

*     Files that will hold results

      OPEN (UNIT = 11, FILE = "quarks.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 12, FILE = "gluons.dat", STATUS = "UNKNOWN")


C*     LO
C          P = 0
C          PAR(11) =  0.157
C          PAR(12) =  1.17
C          PAR(14) =  0.228
C          PAR(21) =  0.527
C          PAR(22) =  1.25
C          PAR(24) =  0.263
C*     NLO MSBAR
C           SCHEME = 'MSBDI'
C           P=1
C           PAR(11) =  0.172
C           PAR(12) =  1.14
C           PAR(14) =  1.93
C           PAR(21) =  0.472
C           PAR(22) =  1.08
C           PAR(24) =  4.45
C*     NLO CSBAR
           P=1
           PAR(11) =  0.167
           PAR(12) =  1.14
           PAR(14) =  1.34
           PAR(21) =  0.535
           PAR(22) =  1.09
           PAR(24) =  1.59
C*     NNLO CSBAR
C          P=2
C          PAR(11) =  0.167
C          PAR(12) =  1.14
C          PAR(14) =  1.17
C          PAR(21) =  0.571
C          PAR(22) =  1.07
C          PAR(24) =  1.39


      DD = DDSTART
      DO 20 PTD = 1, LNPOINTS
        LOGX = LOGXSTART
        DO 10 PTX = 1, LNPOINTS
          DDS(PTD) = DD
          X = 10**LOGX
          XS(PTX) = X


*     Calculating GPD(x, \Delta)

      T = - DD**2

      CALL XSPACE ( HX, X, T)

      POINTS(1, PTD, PTX) = HX(1)
      POINTS(2, PTD, PTX) = HX(2)


 10   LOGX = LOGX + LOGXSTEP
 20   DD = DD + DDSTEP

*     Printing all the results from arrays to files

      DO 40 PTD = 1, LNPOINTS
      DO 30 PTX = 1, LNPOINTS
         WRITE (UNIT=11,FMT=998) DDS(PTD), XS(PTX), POINTS(1,PTD,PTX)
         WRITE (UNIT=12,FMT=998) DDS(PTD), XS(PTX), POINTS(2,PTD,PTX)
 30   CONTINUE
C     WRITE (UNIT=11, FMT=999)
C     WRITE (UNIT=12, FMT=999)
 40   CONTINUE

998   FORMAT (F12.7,3X,F12.7,5X,F12.7)
999   FORMAT (1X)

      STOP
      END
C     ****
