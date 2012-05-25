C     ****h* gepard/radQ.f
C  FILE DESCRIPTION
C    calculation of Q^2 dependence of {\cal H} singlet
C    DVCS form factor
C    Produces data files for Fig. 8. of  KMPKS06b paper
C
C    $Id$
C     *******


C     ****p* radQ.f/RADQ
C  NAME
C     RADQ  --  Program producing data for Figure radQ
C  DESCRIPTION
C    calculation of Q^2 dependence of LO, NLO and NNLO
C    {\cal H} singlet DVCS form factor in CSBAR scheme
C    for collider kinematics
C    Produces data files for Fig. 8? of  KMPKS06b
C  OUTPUT
C       radQ[0-1].dat  --  2 files (one for each panel of Figure
C                            radQ) with 12 sets of point coordinates
C                            (\xi, K_\lambda) or (\xi, \delta\varphi),
C                            corresponding to 12 lines on a graph
C  IDENTIFIERS
C          LNPOINTS -- number of points on each line
C              Q2S -- array(LNPOINTS), values of Q2
C           POINTS -- array(2, 12, LNPOINTS), holds results for
C                 two panels (Fig 8?, in order  a, b)
C                 12 lines depicted on these panels:
C                   1-3 LO, NLO, NNLO 5e-1 (x 0.9 )
C                   3-6 LO, NLO, NNLO 5e-2
C                   6-9 LO, NLO, NNLO 1e-3 (x 1.4)
C                   9-12 LO, NLO, NNLO 1e-5 (x 2)
C
C               PP -- approximation order N^{P}LO P=0,1,2
C
C  CHILDREN
C      READPAR, INIT, CFFF, DCARG
C  SYNOPSIS

      PROGRAM RADQ

C  SOURCE
C
      IMPLICIT NONE
      INTEGER PT, LNPOINTS, LN, NXI, PANEL
      DOUBLE PRECISION Q2START, Q2END, Q2STEP, MP
      DOUBLE PRECISION COLQ2START, COLQ2END
      DOUBLE PRECISION DCARG
      PARAMETER ( LNPOINTS = 40 )
      DOUBLE PRECISION POINTS(0:3, 12, LNPOINTS)
      DOUBLE PRECISION Q2S(LNPOINTS), XIS(4), FACT(4), ADDF(4)
      PARAMETER ( COLQ2START = 2.5d0, COLQ2END =  100.0d0 )
*     Proton mass:
      PARAMETER (MP = 0.938272d0 )

#include "../header.F"

      WRITE (*,*) 'Warning: This one is still broken!'
      PROCESS = 'DVCS'
      FFTYPE = 'SINGLET'

      CALL READPAR
    
#include "ansatz.F"
      ANSATZ = 'FITBP'
* ----  "SOFT"   --------------------
      PAR(21) = 0.3d0
      PAR(22) = PAR(12) - 0.2d0
      PAR(11) = (2.0d0/3.0d0) - PAR(21)


      DEL2 = -0.25d0
      SCHEME = 'CSBAR'
      NQS = 1


*     Files that will hold results

      OPEN (UNIT = 10, FILE = "radQ0.dat", STATUS = "UNKNOWN")
      OPEN (UNIT = 11, FILE = "radQ1.dat", STATUS = "UNKNOWN")

      DO 5 PANEL = 0, 1
  5         WRITE (10 + PANEL, *) '# Output of radQ.f. See prolog of 
     & that program'


*   We do collider kinematics
      Q2START = COLQ2START
      Q2END = COLQ2END
      Q2STEP = (Q2END - Q2START) / (LNPOINTS - 1)

*   \xi points
      XIS(1) = 5.0d-1
      FACT(1) = 0.9d0
      ADDF(1) = 0.45d0
      XIS(2) = 5.0d-2
      FACT(2) = 1.0d0
      ADDF(2) = 0.0d0
      XIS(3) = 1.0d-3
      FACT(3) = 1.4d0
      ADDF(3) = 0.0d0
      XIS(4) = 1.0d-5
      FACT(4) = 2.0d0
      ADDF(4) = 0.25d0


*     Looping over Q2's

      DO 30 PT = 1, LNPOINTS

        Q2 = Q2START + (PT - 1) * Q2STEP
        Q2S(PT) = Q2
        QS(1) = Q2
        P = 2
        CALL INIT


*       Looping over \xi's

        DO 20 NXI = 1, 4
          
          XI = XIS( NXI )

*         Looping over LO, NLO and NNLO ...

          DO 10 P = 0, 2

            CALL EVOLC(1,1)
            CALL CFFF 

*         ... and saving them to array
            
            LN = 3 * (NXI - 1) + P + 1
            POINTS(0, LN, PT) = XIS(NXI) * ABS(CFF(P)) * FACT(NXI)
            POINTS(1, LN, PT) = DCARG(CFF(P)) + ADDF(NXI)

 10       CONTINUE
 20     CONTINUE
 30   CONTINUE

*     Printing all the results from arrays to files

      DO 60 NXI = 1, 4
        DO 50  P = 0, 2
          LN = 3 * (NXI - 1) + P + 1
          DO 40 PT = 1, LNPOINTS
             WRITE (UNIT=10,FMT=998) Q2S(PT), POINTS(0,LN,PT)
             WRITE (UNIT=11,FMT=998) Q2S(PT), POINTS(1,LN,PT)
 40       CONTINUE
*     Empty line is put between the sets
          WRITE (UNIT=10, FMT=999)
          WRITE (UNIT=11, FMT=999)
 50     CONTINUE
 60   CONTINUE

      CLOSE(10)
      CLOSE(11)

998   FORMAT (F6.1,5X,F12.7)
999   FORMAT (1X)

      STOP
      END
C     ****
