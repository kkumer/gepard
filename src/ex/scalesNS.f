C     ****h* gepard/scalesNS.f
C  FILE DESCRIPTION
C    calculation of dependance on \mu_f and \mu_r of {\cal H} non-singlet
C    DVCS form factor
C    Produces data files for Table 1. of  KMPKS06b paper
C
C    $Id$
C     *******


C     ****p* scalesNS.f/SCALESNS
C  NAME
C     SCALESNS  --  Program producing data for Figure scalesNS
C  DESCRIPTION
C    calculation of dependance on \mu_f and \mu_r of {\cal H} non-singlet
C    DVCS form factor
C    Produces data files for Table 1. of  KMPKS06b
C  OUTPUT
C       scalesNS.dat  --  percent variations of H_NS
C
C  IDENTIFIERS
C
C  CHILDREN
C      READPAR, INIT, CFFF, DCARG
C  SOURCE
C

      PROGRAM SCALESNS

      IMPLICIT NONE
      INTEGER PT, NPOINTS, LN, NDEL
      DOUBLE PRECISION DCARG, MP, FACF, FACR
      PARAMETER ( NPOINTS = 4 )
      DOUBLE PRECISION XIS(NPOINTS), POINTS(0:1, 5, NPOINTS)
      DOUBLE COMPLEX CFFD, CFFM, CFFU
*     Proton mass:
      PARAMETER (MP = 0.938272d0 )

      INCLUDE '../header.f'

      DATA XIS / 0.05D0, 0.1D0, 0.25D0, 0.5D0 /

      CALL READPAR
    
      NF = 4
      ANSATZ = 'NSFIT'
      PROCESS = 'NSDVCS'

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
* ----  41 ND    (valence) --
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


*     File that will hold results

      OPEN (UNIT = 10, FILE = "scalesNS.dat", STATUS = "UNKNOWN")

      WRITE (10, *) '# Output of scalesNS.f. See prolog of that program'

*     Scales 

      Q2 = 4.0d0
      DEL2 = -0.25d0

*     Looping over two different ansaetze

      DO 20 NDEL = 0, 1

      IF (NDEL .EQ. 0) THEN
*     "soft"
            PAR(11) = 0.0d0
      ELSE
*     "hard"
            PAR(11) = 4.0d0 / 15.0d0
      END IF

      DO 10 PT = 1, NPOINTS
        XI = XIS(PT)

*     Looping over five table rows

      DO 10 LN = 1, 5

*     Recognizing parameters for each line

      IF (LN .EQ. 1) THEN
            P = 0
            SCHEME = 'CSBAR'
            FACF = 2.0
            FACR = 1.0
      ELSEIF (LN .EQ. 2) THEN
            P = 1
            SCHEME = 'MSBND'
            FACF = 2.0
            FACR = 1.0
      ELSEIF (LN .EQ. 3) THEN
            P = 1
            SCHEME = 'CSBAR'
            FACF = 2.0
            FACR = 1.0
      ELSEIF (LN .EQ. 4) THEN
            P = 1
            SCHEME = 'MSBND'
            FACF = 1.0
            FACR = 2.0
      ELSE
            P = 1
            SCHEME = 'CSBAR'
            FACF = 1.0
            FACR = 2.0
      END IF


*     Calculating three CFFs needed for present line and point ...

      RR2 = 1.0d0
      RF2 = 1.0d0
      CALL INIT
      CALL CFFF 
      CFFM = CFF(P)

      RR2 = 1.0d0 / FACR
      RF2 = 1.0d0 / FACF
      CALL INIT
      CALL CFFF 
      CFFD = CFF(P)

      RR2 = FACR
      RF2 = FACF
      CALL INIT
      CALL CFFF 
      CFFU = CFF(P)
      
*     ... and saving them to array

      POINTS(NDEL, LN, PT) = ( ABS(CFFD) - ABS(CFFU) ) / 
     &        ABS(CFFM) * 100.0d0

 10   CONTINUE

*     Printing all the results for this table to file

 20   CONTINUE

      DO 30 LN = 1, 5
         WRITE (UNIT=10,FMT=999) ( POINTS(0,LN,PT),
     &  POINTS(1,LN,PT) , PT=1,NPOINTS )
 30   CONTINUE

999   FORMAT (1X, 5(F5.1, ' [', F5.1, ']', 2X))

      STOP
      END
C     ****
