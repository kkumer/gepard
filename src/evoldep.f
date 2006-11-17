C     ****h* gepard/evoldep.f
C  FILE DESCRIPTION
C    calculation of relative effects of various evolutions
C
C    $Id: radcorr.f 33 2006-07-26 18:42:37Z kuk05260 $
C     *******


C     ****p* evoldep.f/EVOLDEP
C  NAME
C     EVOLDEP -- Main program
C  DESCRIPTION
C    calculation of relative effects of various evolutions
C    in MSBAR scheme: LO, NLO, NLO non-diagonal, ...
C  OUTPUT
C        EVOL.dat  --  file with 2 sets of point coordinates
C                      corresponding to 2 lines on a graph
C  IDENTIFIERS
C          NPOINTS -- number of points on each line
C       LOGXISTART -- log(\xi) where line begins
C         LOGXIEND -- log(\xi) where line ends
C              XIS -- array(NPOINTS), values of XI
C            PREDS -- array(6, NPOINTS), holds results for
C                     each prediction calculated
C                 Order is:  1.  TOTAL predicition (NLO with NLO ND evol.)
C                            2.  LO prediction
C                            3.  LO prediction with LO evolution
C                            4.  NLO correction  with LO evolution
C                            5.  NLO corr.  with LO+NLO diagonal evolution
C                            6.  NLO corr.  with LO+NLO diag.+ NLO non-diagonal evolution
C
C                P -- approximation order N^{P}LO P=0,1,2
C
C  CHILDREN
C      READPAR, INIT, CFFF
C  SOURCE
C

      PROGRAM EVOLDEP

      IMPLICIT NONE
      INTEGER SPEED, ACC, P, NF
      CHARACTER SCHEME*5, ANSATZ*6
      DOUBLE PRECISION XI, DEL2, Q2, Q02
      DOUBLE COMPLEX CFF(0:2)
      INTEGER PT, NPOINTS, LN, NDEL, CZERO
      DOUBLE PRECISION LOGXI, LOGXISTART, LOGXIEND, LOGXISTEP
      PARAMETER ( NPOINTS = 40 )
      DOUBLE COMPLEX PREDS(6, NPOINTS)
      DOUBLE PRECISION XIS(NPOINTS), RES, TOT, C0
      PARAMETER ( LOGXISTART = -5.0d0, LOGXIEND = -0.30103d0,
     &       LOGXISTEP = (LOGXIEND - LOGXISTART) / (NPOINTS - 1)  )

*     Output common-blocks 

      COMMON / PARINT /  SPEED, ACC, P, NF
      COMMON / PARCHR /  SCHEME, ANSATZ

      COMMON / SWITCH /  CZERO

      COMMON / KINEMATICS /  XI, DEL2, Q2, Q02
      COMMON / CFF        /  CFF

      CALL READPAR

*   File that will hold results

      OPEN (UNIT = 11, FILE = "EVOL.DAT", STATUS = "UNKNOWN")

      DO 5 NDEL = 1, 1
  5         WRITE (10 + NDEL, *) '# Output of evoldep.f. See prolog of 
     & that program'

*   Fixed things
      Q2 = 10.0d0
      ANSATZ = 'NSTOY'

*   Looping over points (XI's) on each line

      LOGXI = LOGXISTART
      DO 20 PT = 1, NPOINTS
      XI = 10**LOGXI
      XIS(PT) = XI

*   Looping over five predictions calculated for each XI

      DO 10 LN = 1, 6

*   Recognizing parameters for each line

      Q02 = 1.0d0
      SCHEME = 'MSBAR'
      CZERO = 0
      IF (LN .EQ. 1) THEN
          P = 1
          SCHEME = 'MSBND'
          CZERO = 1
      ELSEIF (LN .EQ. 2) THEN
          P = 0
          Q02 = 10.0d0
          CZERO = 1
      ELSEIF (LN .EQ. 3) THEN
          P = 0
          CZERO = 1
      ELSEIF (LN .EQ. 4) THEN
          P = 1
          SCHEME = 'MSBLO'
      ELSEIF (LN .EQ. 5) THEN
          P = 1
      ELSEIF (LN .EQ. 6) THEN
          P = 1
          SCHEME = 'MSBND'
      ENDIF

*   Doing calculation ...

      CALL INIT
      CALL CFFF 

*   ... and saving it to arrray

      PREDS(LN, PT) = CFF(P)

 10   CONTINUE

 20   LOGXI = LOGXI + LOGXISTEP

*     Calculating plot points from results in PREDS and printing to files

      DO 40 LN = 2, 6
      DO 30 PT = 1, NPOINTS
        C0 =  ABS(PREDS(2, PT))
        TOT =  ABS(PREDS(1, PT))
        IF (LN .EQ. 2) THEN
          RES = TOT / C0
        ELSE IF (LN .EQ. 3) THEN
          RES = ABS(PREDS(LN,PT)) / C0
        ELSE IF (LN .EQ. 4) THEN
          RES = ABS(PREDS(LN,PT)) / C0
        ELSE
          RES = ABS(PREDS(LN,PT)-PREDS(LN-1,PT)) / C0
        ENDIF
        WRITE (UNIT=11, FMT=998) XIS(PT), RES
        IF (LN .EQ. 5) THEN
          WRITE (17, *) XIS(PT), PREDS(1, PT), PREDS(2,PT), PREDS(3, PT)
          WRITE (17, *) '   ', PREDS(4, PT), PREDS(5,PT), PREDS(6,PT)
        ENDIF
 30   CONTINUE
*     Empty line is put between the sets
      WRITE (UNIT=11, FMT=999)
 40   CONTINUE

998   FORMAT (F12.7,5X,F22.7)
999   FORMAT (1X)

      STOP
      END
C     ****
