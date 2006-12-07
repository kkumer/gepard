C     ****h* gepard/radNLONS.f
C  FILE DESCRIPTION
C    calculation of relative radiative corrections to {\cal H} non-singlet
C    DVCS form factor
C    Produces data files for Fig. 6. of  KMPKS06b paper
C
C    $Id$
C     *******


C     ****p* radNLONS.f/RADNLONS
C  NAME
C     RADNLONS  --  Program producing data for Figure radNLONS
C  DESCRIPTION
C    calculation of relative radiative NLO corrections to {\cal H} non-singlet
C    DVCS form factor in MSBAR and CSBAR schemes
C    Produces data files for Fig. 6 of  KMPKS06b
C  OUTPUT
C       radNLONS[0-3].dat  --  4 files (one for each panel of Figure
C                            radNLONS) with 4 sets of point coordinates
C                            (\xi, K_\lambda) or (\xi, \delta\varphi),
C                            corresponding to 4 lines on a graph
C  IDENTIFIERS
C          NPOINTS -- number of points on each line
C       LOGXISTART -- log(\xi) where line begins
C         LOGXIEND -- log(\xi) where line ends
C              XIS -- array(NPOINTS), values of XI
C           POINTS -- array(0:3, 4, NPOINTS), holds results for
C                 four panels (Fig 7, in order  a, b, c and d) and
C                 4 K_{\lambda,arg} lines depicted on these panels
C                 Order is:  1.  NLO, CSBAR, HARD 
C                            2.  NLO, CSBAR, SOFT
C                            3.  NLO, MSBAR, HARD
C                            4.  NLO, MSBAR, SOFT
C
C                P -- approximation order N^{P}LO P=0,1,2
C
C  CHILDREN
C      READPAR, INIT, CFFF, DCARG
C  SOURCE
C

      PROGRAM RADNLONS

      IMPLICIT NONE
      INTEGER PT, NPOINTS, LN, NDEL
      DOUBLE PRECISION XISTART, XIEND, XISTEP, MP
      DOUBLE PRECISION DCARG
      PARAMETER ( NPOINTS = 40 )
      DOUBLE PRECISION POINTS(0:3, 4, NPOINTS)
      DOUBLE PRECISION XIS(NPOINTS)
      PARAMETER ( XISTART = 0.01d0, XIEND = 0.5d0,
     &       XISTEP = (XIEND - XISTART) / (NPOINTS - 1) )
*     Proton mass:
      PARAMETER (MP = 0.938272d0 )

      INCLUDE '../header.f'


      CALL READPAR
    
      NF = 4
      ANSATZ = 'NSFIT'

*       1  Q02       
        PAR(1) =   2.5d0    
*       2  AS0       
        PAR(2) =   0.05d0
*       3  MU02      
        PAR(3) =   2.5d0 

* ------------ ANSATZ  -------------
* ----  11 NS   --------------------
!        PAR(11) =  see below
*       12 AL0S      
        PAR(12) =  1.1d0 
*       13 ALPS      
        PAR(13) =  0.15d0
*       14 M02S      
        PAR(14) =  (2.0d0 * MP)**2
*       15 DELM2S    
        PAR(15) =  MP**2
*       16 PS        
        PAR(16) =  3.0d0
* ----  21 NG   --------------------
!        PAR(21) =  0.5d0
*       22 AL0G      
!        PAR(22) =  1.0d0 
*       23 ALPG      
        PAR(23) =  0.15d0
*       24 M02G      
        PAR(24) =  (2.0d0 * MP)**2
*       25 DELM2G    
        PAR(25) =  MP**2
*       26 PG        
        PAR(26) =  2.0d0
* ----  31 NU   --------------------
        PAR(31) =  2.0d0 
*       32 AL0U      
        PAR(32) =  0.5d0 
*       33 ALPU      
        PAR(33) =  1.0d0 
*       34 M02U      
        PAR(34) =  (2.0d0 * MP)**2
*       35 DELM2U    
        PAR(35) =  MP**2
*       36 PU        
        PAR(36) =  1.0d0 
* ----  41 ND   --------------------
        PAR(41) =  1.0d0 
*       42 AL0D      
        PAR(42) =  0.5d0 
*       43 ALPD      
        PAR(43) =  1.0d0 
*       44 M02D      
        PAR(44) =  (2.0d0 * MP)**2
*       45 DELM2D    
        PAR(45) =  MP**2
*       46 PD        
        PAR(46) =  1.0d0    


*     Files that will hold results

      OPEN (UNIT = 10, FILE = "radNLONS0.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 11, FILE = "radNLONS1.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 12, FILE = "radNLONS2.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 13, FILE = "radNLONS3.dat", STATUS = "UNKNOWN")

      DO 5 NDEL = 0, 3
  5         WRITE (10 + NDEL, *) '# Output of radNLONS.f. See prolog of 
     & that program'

*     Scales 

      Q2 = 2.5d0

*     Looping over two different walues of \Delta^2

      DO 40 NDEL = 0, 1

      DEL2 = 0.0d0
      IF (NDEL .EQ. 1) DEL2 = -1.0d0

      XI = XISTART
      DO 20 PT = 1, NPOINTS
      XIS(PT) = XI

*     Looping over four lines for each panel

      DO 10 LN = 1, 4

*     Recognizing parameters for each line

      IF (LN .LE. 2) THEN
            SCHEME = 'CSBAR'
      ELSE
            SCHEME = 'MSBAR'
      END IF
      IF ((LN .EQ. 1) .OR. (LN .EQ. 3)) THEN
!          'HARD'
            PAR(11) = 4.0d0 / 15.0d0
      ELSE
!           'SOFT'
            PAR(11) = 0.0d0
      END IF


*     Calculating two CFFs needed for present line and point ...

      P = 1
      CALL INIT
      CALL CFFF 
      P = 0
      CALL CFFF 

*     ... and saving them to array

      POINTS(NDEL, LN, PT) = (ABS(CFF(1)) / ABS(CFF(0)) -
     &        1.0d0) * 100.0d0
      POINTS(2+NDEL, LN, PT) = DCARG(CFF(1)) - DCARG(CFF(0))

 10   CONTINUE

 20   XI = XI + XISTEP

*     Printing all the results from arrays to files

      DO 40 LN = 1, 4
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
