
*   "header" file with definitions

*   Parameters
*   ----------

*   1. Constant parameters loaded from 'GEPARD.INI' or fixed by program

      INTEGER SPEED, ACC, P, NF, CZERO
      DOUBLE PRECISION RF2, RR2
      DOUBLE PRECISION C, PHI, CND, PHIND
      CHARACTER SCHEME*5, ANSATZ*6, PROCESS*6, FFTYPE*10

*   2. Parameters from 'MINUIT.CMD' (candidates for fitting parameters)

      INTEGER NPARMAX
      PARAMETER (NPARMAX = 50)
      DOUBLE PRECISION PAR(NPARMAX)

*   3. Kinematics

      DOUBLE PRECISION XI, DEL2, Q2

*   4. Other

      DOUBLE PRECISION PI, CF, CA, TF
      PARAMETER ( PI = 3.1415 92653 58979 D0 )
      PARAMETER ( CF = 4.0d0 / 3.0d0, CA = 3.0d0, TF = 0.5d0 )

*     - QCD beta function
      INTEGER NFMIN, NFMAX
      PARAMETER (NFMIN = 3, NFMAX = 6)
      DOUBLE PRECISION BETA0, BETA1, BETA2, BETA3

*     - Mellin-Barnes integration contour points
      INTEGER NPTS, NPTSMAX
      PARAMETER (NPTSMAX = 768)
      DOUBLE PRECISION Y(NPTSMAX), WG(NPTSMAX)
      DOUBLE COMPLEX N(2,NPTSMAX)

*     - Values on a particular contour point
      DOUBLE COMPLEX S1, S2, S3, S4
      DOUBLE COMPLEX CNS0, CNS1, CNS2

*     - Values on the whole contour
      DOUBLE COMPLEX CDISNS1(NPTSMAX,0:2)
      DOUBLE COMPLEX CDIS1(NPTSMAX,0:2,2), CDIS2(NPTSMAX,0:2,2)
      DOUBLE COMPLEX BIGCNS(NPTSMAX,0:2)
      DOUBLE COMPLEX BIGC(2,NPTSMAX,0:2,2)
      DOUBLE COMPLEX BIGCF2(NPTSMAX,0:2,2)
      DOUBLE COMPLEX GAMNS(NPTSMAX,0:2)
      DOUBLE COMPLEX GAM(2,NPTSMAX,0:2,2,2)

*     - Final observables
      DOUBLE PRECISION F2(0:2)
      DOUBLE COMPLEX CFF(0:2)
        

*   Common blocks
*   -------------
      
*   1. Constant parameters loaded from 'GEPARD.INI' or fixed by program

      COMMON / PARINT /  SPEED, ACC, P, NF, CZERO
      COMMON / PARFLT /  RF2, RR2
      COMMON / MBCONT /  C, PHI, CND, PHIND
      COMMON / PARCHR /  SCHEME, ANSATZ, PROCESS, FFTYPE

*   2. Parameters from 'MINUIT.CMD' (candidates for fitting parameters)

      COMMON / PAR    /  PAR

*   3. Kinematics

      COMMON / KINEMATICS /  XI, DEL2, Q2

*   4. Other

*     - QCD beta function
      COMMON / BETABLK / BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX),
     &                   BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)

*     - Mellin-Barnes integration contour points
      COMMON / CONTOUR  /  NPTS
      COMMON / POINTS   /  Y, WG
      COMMON / NPOINTS  /  N
      
*     - Values on a particular contour point
      COMMON / HARMONIC /  S1, S2, S3, S4

*     - Values on the whole contour
      COMMON / CDISNS   /  CDISNS1
      COMMON / CDIS     /  CDIS1, CDIS2
      COMMON / BIGCNS   /  BIGCNS
      COMMON / BIGC     /  BIGC
      COMMON / GAMNS    /  GAMNS
      COMMON / GAM      /  GAM
      COMMON / BIGCF2   /  BIGCF2

*     - FInal observables
      COMMON / CFF      /  CFF
      COMMON / F2       /  F2
