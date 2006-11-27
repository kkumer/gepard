C     ****h* gepard/radcorr2.f
C  FILE DESCRIPTION
C    calculation of relative radiative corrections to {\cal H} singlet
C    DVCS form factor
C    Produces data files for Fig. 1. of  KMPKS06 Letter
C
C    $Id: radcorr.f 33 2006-07-26 18:42:37Z kuk05260 $
C     *******


C     ****p* radcorr2.f/RADCORR2
C  NAME
C     RADCORR  --  Main program
C  DESCRIPTION
C    calculation of relative radiative corrections to {\cal H} singlet
C    DVCS form factor
C    Produces data files for Fig. 1. of  KMPKS06 Letter
C  OUTPUT
C       FIG1[A-D].dat  --  4 files with six sets of point coordinates
C                          (\xi, K_\lambda) or (\xi, K_arg),
C                          corresponding to six lines on a graph
C  IDENTIFIERS
C          NPOINTS -- number of points on each line
C       LOGXISTART -- log(\xi) where line begins
C         LOGXIEND -- log(\xi) where line ends
C              XIS -- array(NPOINTS), values of XI
C           POINTS -- array(4, 6, NPOINTS), holds results for
C                 four panels (Fig 1, in order  a, b, c and d) and
C                 six K_{\lambda,arg} lines depicted on these panels
C                 Order is:  1. NNLO, CSBAR, HARD 
C                            2.  NLO, CSBAR, HARD
C                            3.  NLO, MSBAR, HARD
C                            4. NNLO, CSBAR, SOFT
C                            5.  NLO, CSBAR, SOFT
C                            6.  NLO, MSBAR, SOFT
C
C                P -- approximation order N^{P}LO P=0,1,2
C
C  CHILDREN
C      READPAR, INIT, CFFF, DCARG
C  SOURCE
C

      PROGRAM RADCORR2

      IMPLICIT NONE
      INTEGER PT, NPOINTS, LN, NDEL
      DOUBLE PRECISION LOGXI, LOGXISTART, LOGXIEND, LOGXISTEP
      DOUBLE PRECISION DCARG
      PARAMETER ( NPOINTS = 40 )
      DOUBLE PRECISION POINTS(4, 6, NPOINTS)
      DOUBLE PRECISION XIS(NPOINTS)
      PARAMETER ( LOGXISTART = -5.0d0, LOGXIEND = -0.30103d0,
     &       LOGXISTEP = (LOGXIEND - LOGXISTART) / (NPOINTS - 1)  )
      INCLUDE 'header.f'


      CALL READPAR
    
      ANSATZ = 'FIT'

*       1  'Q02       '   7.0       .1  
        PAR(1) =   2.5    
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

      OPEN (UNIT = 11, FILE = "FIG1A.DAT", STATUS = "UNKNOWN")
      OPEN (UNIT = 12, FILE = "FIG1B.DAT", STATUS = "UNKNOWN")
      OPEN (UNIT = 13, FILE = "FIG1C.DAT", STATUS = "UNKNOWN")
      OPEN (UNIT = 14, FILE = "FIG1D.DAT", STATUS = "UNKNOWN")

      DO 5 NDEL = 1, 4
  5         WRITE (10 + NDEL, *) '# Output of radcorr2.f. See prolog of 
     & that program'

*     Scales 

      Q2 = 2.5d0

*     Looping over two different walues of \Delta^2

      DO 40 NDEL = 0, 1

      DEL2 = 0.0d0
      IF (NDEL .EQ. 1) DEL2 = -0.5d0

      LOGXI = LOGXISTART
      DO 20 PT = 1, NPOINTS
      XI = 10**LOGXI
      XIS(PT) = XI

*     Looping over six lines for each panel

      DO 10 LN = 1, 6

*     Recognizing parameters for each line

      IF ((LN .EQ. 1) .OR. (LN .EQ. 4)) THEN
            P = 2
      ELSE
            P = 1
      END IF
      IF ((LN .EQ. 3) .OR. (LN .EQ. 6)) THEN
            SCHEME = 'MSBAR'
      ELSE
            SCHEME = 'CSBAR'
      END IF
      IF (LN .LE. 3) THEN
            PAR(21) = 0.4d0
            PAR(22) = PAR(12) + 0.1d0
!            ANSATZ = 'HARD'
      ELSE
            PAR(21) = 0.3d0
            PAR(22) = PAR(12)
!            ANSATZ = 'SOFT'
      END IF
      PAR(11) = (2.0d0/3.0d0) - PAR(21)


      CALL INIT

*     Calculating two CFFs needed for present line and point ...

      CALL CFFF 
      P = P - 1
      CALL CFFF 

*     ... and saving them to array

      POINTS(1+NDEL, LN, PT) = LOG(ABS( CFF(P+1) )) / LOG(ABS( CFF(P) ))
      POINTS(3+NDEL, LN, PT) = DCARG( CFF(P+1) ) / DCARG( CFF(P) )

 10   CONTINUE

 20   LOGXI = LOGXI + LOGXISTEP

*     Printing all the results from arrays to files

      DO 40 LN = 1, 6
      DO 30 PT = 1, NPOINTS
         WRITE (UNIT=11+NDEL,FMT=998) XIS(PT), POINTS(1+NDEL, LN, PT)
         WRITE (UNIT=13+NDEL,FMT=998) XIS(PT), POINTS(3+NDEL, LN, PT)
 30   CONTINUE
*     Empty line is put between the sets
      WRITE (UNIT=11+NDEL, FMT=999)
      WRITE (UNIT=13+NDEL, FMT=999)
 40   CONTINUE

998   FORMAT (F12.7,5X,F12.7)
999   FORMAT (1X)

      STOP
      END
C     ****
