C     ****h* gepard/anatomy.f
C  FILE DESCRIPTION
C    calculation of relative effects of various evolutions
C
C    $Id: radcorr.f 33 2006-07-26 18:42:37Z kuk05260 $
C     *******


C     ****p* anatomy.f/ANATOMY
C  NAME
C     ANATOMY -- Main program
C  DESCRIPTION
C    calculation of relative effects of various evolutions
C    in MSBAR scheme: LO, NLO, NLO non-diagonal, ... (singlet)
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

      PROGRAM ANATOMY

      IMPLICIT NONE
      INTEGER PT, NPOINTS, LN, NDEL
      DOUBLE PRECISION LOGXI, LOGXISTART, LOGXIEND, LOGXISTEP
      PARAMETER ( NPOINTS = 40 )
      DOUBLE COMPLEX PREDS(6, NPOINTS), TOT
      DOUBLE PRECISION XIS(NPOINTS), MODUL, PHASE, Q2EXP, MP
      DOUBLE PRECISION DCARG
      CHARACTER SUBANSATZ*4
      PARAMETER ( LOGXISTART = -5.0d0, LOGXIEND = -0.30103d0,
     &       LOGXISTEP = (LOGXIEND - LOGXISTART) / (NPOINTS - 1)  )
      PARAMETER (MP = 0.938272d0 )
      INCLUDE '../header.f'

      PROCESS = 'DVCS'
      FFTYPE = 'SINGLET'

      CALL READPAR

      OPEN (UNIT = 10, FILE = "anatomy0.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 11, FILE = "anatomy1.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 12, FILE = "anatomy2.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 13, FILE = "anatomy3.dat", STATUS = "UNKNOWN")

      DO 5 NDEL = 0, 3
  5         WRITE (10 + NDEL, *) '# Output of anatomy.f. See prolog of 
     & that program'

      INCLUDE 'ansatz.f'
      ANSATZ = 'FIT'



      Q02 = 2.5d0
      Q2EXP = 10.0d0
      DEL2 =  -0.25d0
      NQS = 1


*     Looping over two different ansaetze

      DO 40 NDEL = 0, 1
        IF ( NDEL .EQ.  1 ) THEN
!         'HARD'
           PAR(21) = 0.4d0
           PAR(22) = PAR(12) + 0.05d0
        ELSE
!         'SOFT'
          PAR(21) = 0.3d0
          PAR(22) = PAR(12) - 0.2d0
        END IF
        PAR(11) = (2.0d0/3.0d0) - PAR(21)


*   Looping over points (XI's) on each line

      LOGXI = LOGXISTART
      DO 20 PT = 1, NPOINTS
      XI = 10**LOGXI
      XIS(PT) = XI

*   Looping over five predictions calculated for each XI

      DO 10 LN = 1, 6

*   Recognizing parameters for each line

* 1st line is the default total result:
      P = 1
      PAR(1) = Q02
      Q2 = Q2EXP
      SCHEME = 'MSBND'
      CZERO = 1
      IF (LN .EQ. 2) THEN
          P = 0
          PAR(1) = Q2
      ELSEIF (LN .EQ. 3) THEN
          P = 0
      ELSEIF (LN .EQ. 4) THEN
          CZERO = 0
          SCHEME = 'MSBLO'
      ELSEIF (LN .EQ. 5) THEN
          CZERO = 0
          SCHEME = 'MSBAR'
      ELSEIF (LN .EQ. 6) THEN
          CZERO = 0
          SCHEME = 'MSBND'
      ENDIF

*   Doing calculation ...
  
      CALL INIT
      QS(1) = Q2
      CALL EVOLC(1)
      CALL GETMBGPD
      CALL CFFF 

*   ... and saving it to arrray

      PREDS(LN, PT) = CFF(P)

 10   CONTINUE

 20   LOGXI = LOGXI + LOGXISTEP

*     Calculating plot points from results in PREDS and printing to files

      DO 40 LN = 1, 3
      DO 30 PT = 1, NPOINTS
        TOT = PREDS(1, PT)
        IF (LN .EQ. 1) THEN
          MODUL = ABS( PREDS(4,PT) ) / ABS( TOT )
          PHASE = DCARG( (TOT - PREDS(4,PT)) / TOT )
        ELSEIF (LN .EQ. 2) THEN
          MODUL = ABS( PREDS(5,PT)- PREDS(4,PT)) / ABS( TOT )
          PHASE = DCARG( (TOT-PREDS(5,PT)+PREDS(4,PT)) / TOT )
        ELSEIF (LN .EQ. 3) THEN
          MODUL = ABS( PREDS(6,PT)-PREDS(5,PT) ) / ABS( TOT )
          PHASE = DCARG( (TOT-PREDS(6,PT)+PREDS(5,PT)) / TOT )
        ENDIF
        WRITE (UNIT=10+NDEL, FMT=998) XIS(PT), MODUL
        WRITE (UNIT=12+NDEL, FMT=998) XIS(PT), PHASE
 30   CONTINUE
*     Empty line is put between the sets
      WRITE (UNIT=10+NDEL, FMT=999)
      WRITE (UNIT=12+NDEL, FMT=999)
 40   CONTINUE

998   FORMAT (F12.7,5X,F22.7)
999   FORMAT (1X)

      STOP
      END
C     ****
