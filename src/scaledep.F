C     ****h* gepard/scaledep.F
C  FILE DESCRIPTION
C    calculation of relative change of scale dependence of {\cal H}
C    singlet DVCS form factor.
C    Produces data files for Fig. 2. of KMPKS06 Letter
C
C    $Id$
C     *******


C     ****p* scaledep.F/SCALEDEP
C  NAME
C     SCALEDEP  --  Main program
C  DESCRIPTION
C    calculation of relative change of scale dependene of  {\cal H}
C    singlet DVCS form factor.
C    Produces data files for Fig. 2. of KMPKS06 Letter
C  OUTPUT
C       FIG2[AB].dat  --  2 files with 8 sets of point coordinates
C                          (\xi, \dot{K}_\lambda) or (\xi, \dot{K}_arg),
C                          corresponding to 8 lines on a graph
C  IDENTIFIERS
C          NPOINTS -- number of points on each line
C       LOGXISTART -- log(\xi) where line begins
C         LOGXIEND -- log(\xi) where line ends
C              XIS -- array(NPOINTS), values of XI
C           POINTS -- array(2, 8, NPOINTS), holds results for
C                 two panels (Fig 2, in order  a, b) and
C                 eight \dot{K}_{abs,arg} lines depicted on these panels
C                 Order is:  1. NNLO, DEL2=0,    HARD 
C                            2. NNLO, DEL2=-0.5, HARD
C                            3.  NLO, DEL2=0   , HARD
C                            4.  NLO, DEL2=-0.5, HARD
C                            5. NNLO, DEL2=0,    SOFT 
C                            6. NNLO, DEL2=-0.5, SOFT
C                            7.  NLO, DEL2=0   , SOFT
C                            8.  NLO, DEL2=-0.5, SOFT
C
C                P -- approximation order N^{P}LO P=0,1,2
C             DRVU -- derivative of CFF for given P
C             DRVD -- derivative of CFF for P-1
C  CHILDREN
C      READPAR, INIT, LOGABSH, ARGH, DERIV
C  BUGS
C    - Numerical derivation routine DERIV is not perfect. One needs to
C      be careful with choice of NEVALS and H parameters, and check
C      the errors returned in ERR
C  SOURCE
C

      PROGRAM SCALEDEP

      IMPLICIT NONE
      INTEGER PT, NPOINTS, LN, NDEL
      DOUBLE PRECISION QVAR2
      DOUBLE PRECISION LOGXI, LOGXISTART, LOGXIEND, LOGXISTEP
      INTEGER NEVALS
      DOUBLE PRECISION  LOGABSH, ARGH, H, DRVU, DRVD, ERR
      EXTERNAL LOGABSH, ARGH
      PARAMETER ( NPOINTS = 2 )
      DOUBLE PRECISION POINTS(2, 8, NPOINTS)
      DOUBLE PRECISION XIS(NPOINTS)
      PARAMETER ( NEVALS = 5, H = 0.41d0 )
      PARAMETER ( LOGXISTART = -5.0d0, LOGXIEND = -1.30103d0,
     &       LOGXISTEP = (LOGXIEND - LOGXISTART) / (NPOINTS - 1)  )
#include "header.F"

      PROCESS = 'DVCS'
      FFTYPE = 'SINGLET'

      CALL READPAR

*     Files that will hold results

      OPEN (UNIT = 11, FILE = "FIG2A.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 12, FILE = "FIG2B.dat", STATUS = "UNKNOWN")

      DO 5 NDEL = 1, 2
  5         WRITE (10 + NDEL, *) '# Output of scaledep.F. See prolog of 
     & that program'

*     Scales 

      QVAR2 = 4.0d0

      NQS = 1


      SCHEME = 'CSBAR'
      P = 2
      CALL INIT


*     Looping over NPOINTS points for each line

      LOGXI = LOGXISTART
      DO 20 PT = 1, NPOINTS
      XI = 10**LOGXI
      XIS(PT) = XI

*     Looping over eight lines for each panel

      DO 10 LN = 1, 8

*     Recognizing parameters for each line

      IF ((LN .LE. 2) .OR. (LN .EQ. 5) .OR. (LN .EQ. 6)) THEN
            P = 2
      ELSE
            P = 1
      END IF
      IF (DBLE(LN/2) .EQ. LN/2.0d0) THEN
            NDEL = 1
      ELSE
            NDEL = 0
      END IF
      IF (LN .LE. 4) THEN
            ANSATZ = 'HARD'
      ELSE
            ANSATZ = 'SOFT'
      END IF

      CALL INIT
      MTIND = 0
      IF (NDEL .EQ. 1) MTIND = NMTS + 3
      DEL2 = -MTS(MTIND)

*     Calculating derivatives of CFFs needed for present line and point ...
*     ... and saving them to array (note that P changes value)

      CALL DERIV(LOGABSH, QVAR2, H, NEVALS, DRVU, ERR)
      P = P - 1
      CALL DERIV(LOGABSH, QVAR2, H, NEVALS, DRVD, ERR)
      POINTS(1, LN, PT) = DRVU / DRVD

*     Returning P to the correct value 
      P = P + 1
      CALL DERIV(ARGH, QVAR2, H, NEVALS, DRVU, ERR)
      P = P - 1
      CALL DERIV(ARGH, QVAR2, H, NEVALS, DRVD, ERR)
      POINTS(2, LN, PT) = DRVU / DRVD

 10   CONTINUE

 20   LOGXI = LOGXI + LOGXISTEP

*     Printing all the results from arrays to files

      DO 40 LN = 1, 8
      DO 30 PT = 1, NPOINTS
         WRITE (UNIT=11,FMT=998) XIS(PT), POINTS(1, LN, PT)
         WRITE (UNIT=12,FMT=998) XIS(PT), POINTS(2, LN, PT)
 30   CONTINUE
*     Empty line is put between the sets
      WRITE (UNIT=11, FMT='(1X)')
      WRITE (UNIT=12, FMT='(1X)')
 40   CONTINUE

998   FORMAT (F12.7,5X,F12.7)

      STOP
      END
C     ****


C     ****f* scaledep.F/LOGABSH
C  NAME
C    LOGABSH  --  ln(abs( CFF ))
C  SYNOPSIS
C     DOUBLE PRECISION FUNCTION LOGABSH(QQ2)
C
C     DOUBLE PRECISION QQ2
C  CHILDREN
C       CFFF
C  SOURCE
C

      DOUBLE PRECISION FUNCTION LOGABSH(QQ2)

      DOUBLE PRECISION QQ2
#include "header.F"

      Q2 = QQ2
      QS(1) = Q2
      CALL EVOLC(1)

      CALL CFFF 

      LOGABSH = LOG(ABS( CFF(P) ))

      RETURN
      END
C     ****


C     ****f* scaledep.F/ARGH
C  NAME
C    ARGH  --  arg( CFF )
C  SYNOPSIS
C     DOUBLE PRECISION FUNCTION ARGH(QQ2)
C
C     DOUBLE PRECISION QQ2
C  CHILDREN
C       CFFF, DCARG
C  SOURCE
C

      DOUBLE PRECISION FUNCTION ARGH(QQ2)

      DOUBLE PRECISION QQ2
      DOUBLE PRECISION DCARG
#include "header.F"

      Q2 = QQ2
      QS(1) = Q2
      CALL EVOLC(1)

      CALL CFFF 

      ARGH = DCARG( CFF(P) )

      RETURN
      END
C     ****
