C     ****h* gepard/scaledep.f
C  FILE DESCRIPTION
C    calculation of relative change of scale dependence of {\cal H}
C    singlet DVCS form factor.
C    Produces data files for Fig. 2. of KMPKS06 Letter
C
C    $Id: radcorr.f,v 1.1 2006-04-13 19:27:11+02 kkumer Exp kkumer $
C     *******


C     ****p* scaledep.f/SCALEDEP
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
C      LOGABSH, ARGH, DERIV
C  BUGS
C    - Numerical derivation routine DERIV is not perfect. One needs to
C      be careful with choice of NEVALS and H parameters, and check
C      the errors returned in ERR
C    - One should probably control numerical integration
C      parameters (errors, limits, ...) from here (or in some input file), 
C      instead of hard-wiring them in CFFF.
C  SOURCE
C

      PROGRAM SCALEDEP

      IMPLICIT NONE
      INTEGER P, PT, NPOINTS, LN, NDEL
      DOUBLE PRECISION XI, DEL2, Q2, Q02, QVAR2
      DOUBLE COMPLEX CFF(0:2)
      CHARACTER SCHEME*5, ANSATZ*6
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

*     Output common-blocks 

      COMMON / KINEMATICS /  XI, DEL2, Q2, Q02
      COMMON / APPROX     /  P
      COMMON / LABELS     /  SCHEME, ANSATZ
      COMMON / CFF        /  CFF

*     Files that will hold results

      OPEN (UNIT = 11, FILE = "FIG2A.DAT", STATUS = "NEW")
      OPEN (UNIT = 12, FILE = "FIG2B.DAT", STATUS = "NEW")

      DO 5 NDEL = 1, 2
  5         WRITE (10 + NDEL, *) '# Output of scaledep.f. See prolog of 
     & that program'

*     Scales 

      QVAR2 = 4.0d0
      Q02 = 1.0d0

*     Looping over NPOINTS points for each line

      SCHEME = 'CSBAR'
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
            DEL2 = -0.5d0
      ELSE
            DEL2 = 0.0d0
      END IF
      IF (LN .LE. 4) THEN
            ANSATZ = 'HARD'
      ELSE
            ANSATZ = 'SOFT'
      END IF

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


C     ****f* scaledep.f/LOGABSH
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
      DOUBLE PRECISION XI, DEL2, Q2, Q02
      INTEGER P
      DOUBLE COMPLEX CFF(0:2)
      
      COMMON / KINEMATICS /  XI, DEL2, Q2, Q02
      COMMON / APPROX     /  P
      COMMON / CFF        /  CFF

      Q2 = QQ2

      CALL CFFF 

      LOGABSH = LOG(ABS( CFF(P) ))

      RETURN
      END
C     ****


C     ****f* scaledep.f/ARGH
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
      DOUBLE PRECISION XI, DEL2, Q2, Q02, DCARG
      INTEGER P
      DOUBLE COMPLEX CFF(0:2)
      
      COMMON / KINEMATICS /  XI, DEL2, Q2, Q02
      COMMON / APPROX     /  P
      COMMON / CFF        /  CFF

      Q2 = QQ2

      CALL CFFF 

      ARGH = DCARG( CFF(P) )

      RETURN
      END
C     ****
