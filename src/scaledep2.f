C     ****h* gepard/scaledep2.f
C  FILE DESCRIPTION
C    calculation of relative change of scale dependence of {\cal H}
C    singlet DVCS form factor.
C    Produces data files for Fig. 2. of KMPKS06 Letter
C
C    $Id: scaledep.f 33 2006-07-26 18:42:37Z kuk05260 $
C     *******


C     ****p* scaledep2.f/SCALEDEP2
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
      PARAMETER ( NPOINTS = 40 )
      DOUBLE PRECISION POINTS(2, 8, NPOINTS)
      DOUBLE PRECISION XIS(NPOINTS)
      PARAMETER ( NEVALS = 5, H = 0.41d0 )
      PARAMETER ( LOGXISTART = -5.0d0, LOGXIEND = -0.30103d0,
     &       LOGXISTEP = (LOGXIEND - LOGXISTART) / (NPOINTS - 1)  )
      INCLUDE 'header.f'


      CALL READPAR
      ANSATZ = 'FIT'

*       1  'Q02       '   7.0       .1  
        PAR(1) =   1.0    
*       2  'AS0       '   0.0536    .1  
        PAR(2) =   0.05
*       3  'MU02      '   2.5       .1  
        PAR(3) =   2.5    
*       11 'NS        '   0.2       .1  
        PAR(11) =  0.2    
*       12 'AL0S      '   1.1       .1   
        PAR(12) =  1.1    
*       13 'ALPS      '   0.25      .1  
        PAR(13) =  0.25   
*       14 'M02S      '   3.53      .1  
        PAR(14) =  3.53   
*       15 'DELM2S    '   0.88      .1  
        PAR(15) =  0.88   
*       16 'PS        '   3.0       .1  
        PAR(16) =  3.0    
*       21 'NG        '   0.5       .1  
        PAR(21) =  0.5    
*       22 'AL0G      '   1.0       .1   
        PAR(22) =  1.0    
*       23 'ALPG      '   0.25      .1  
        PAR(23) =  0.25   
*       24 'M02G      '   3.53      .1  
        PAR(24) =  3.53   
*       25 'DELM2G    '   0.88      .1  
        PAR(25) =  0.88   
*       26 'PG        '   2.0       .1  
        PAR(26) =  2.0    
*       31 'NU        '   0.0       .1  
        PAR(31) =  2.0    
*       32 'AL0U      '   0.5       .1   
        PAR(32) =  0.5    
*       33 'ALPU      '   1.0       .1  
        PAR(33) =  1.0    
*       34 'M02U      '   3.53      .1  
        PAR(34) =  3.53   
*       35 'DELM2U    '   0.88      .1  
        PAR(35) =  0.88   
*       36 'PU        '   1.0       .1  
        PAR(36) =  1.0    
*       41 'ND        '   0.0       .1  
        PAR(41) =  1.0    
*       42 'AL0D      '   0.5       .1   
        PAR(42) =  0.5    
*       43 'ALPD      '   1.0       .1  
        PAR(43) =  1.0    
*       44 'M02D      '   3.53      .1  
        PAR(44) =  3.53   
*       45 'DELM2D    '   0.88      .1  
        PAR(45) =  0.88   
*       46 'PD        '   1.0       .1  
        PAR(46) =  1.0    


*     Files that will hold results

      OPEN (UNIT = 11, FILE = "FIG2A.DAT", STATUS = "UNKNOWN")
      OPEN (UNIT = 12, FILE = "FIG2B.DAT", STATUS = "UNKNOWN")

      DO 5 NDEL = 1, 2
  5         WRITE (10 + NDEL, *) '# Output of scaledep2.f. See prolog of
     & that program'

*     Scales 

      QVAR2 = 4.0d0


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
            DEL2 = -0.5d0
      ELSE
            DEL2 = 0.0d0
      END IF
      IF (LN .LE. 4) THEN
            PAR(21) = 0.4d0
            PAR(22) = PAR(12) + 0.1d0
!            ANSATZ = 'HARD'
      ELSE
            PAR(21) = 0.3d0
            PAR(22) = PAR(12)
!            ANSATZ = 'SOFT'
      END IF
      PAR(11) = (2.0d0/3.0d0) - PAR(21)


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


C     ****f* scaledep2.f/LOGABSH
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
      INCLUDE 'header.f'

      Q2 = QQ2

      CALL CFFF 

      LOGABSH = LOG(ABS( CFF(P) ))

      RETURN
      END
C     ****


C     ****f* scaledep2.f/ARGH
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
      INCLUDE 'header.f'

      Q2 = QQ2

      CALL CFFF 

      ARGH = DCARG( CFF(P) )

      RETURN
      END
C     ****
