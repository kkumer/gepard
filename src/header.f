
*   "header" file with definitions

*   Parameters
*   ----------

*   1. Constant parameters loaded from 'GEPARD.INI'

      INTEGER SPEED, ACC, P, NF, CZERO
      DOUBLE PRECISION RF2, RR2
      DOUBLE PRECISION C, PHI, CND, PHIND
      CHARACTER SCHEME*5, ANSATZ*6

*   2. Parameters from 'MINUIT.CMD' (candidates for fitting parameters)

      INTEGER NPARMAX
      PARAMETER (NPARMAX = 50)
      DOUBLE PRECISION PAR(NPARMAX)

*   3. Kinematics

      DOUBLE PRECISION XI, DEL2, Q2

*   4. Other

      DOUBLE PRECISION PI
      PARAMETER ( PI = 3.1415 92653 58979 D0 )

*     - QCD beta function
      INTEGER NFMIN, NFMAX
      PARAMETER (NFMIN = 3, NFMAX = 6)
      DOUBLE PRECISION BETA0, BETA1, BETA2, BETA3

*     - Mellin-Barnes integration contour points
      INTEGER NPTS, NPTSMAX
      PARAMETER (NPTSMAX = 768)
      DOUBLE PRECISION Y(NPTSMAX), WG(NPTSMAX)
      DOUBLE COMPLEX N(NPTSMAX)

*     - Values on a particular contour point
      DOUBLE COMPLEX S1, S2, S3, S4
      DOUBLE COMPLEX GAM0(2,2), GAM1(2,2), GAM2(2,2)
      DOUBLE COMPLEX GAMNS0, GAMNS1, GAMNS2
      DOUBLE COMPLEX C0(2), C1(2), C2(2)

*     - Values on the whole contour
      DOUBLE COMPLEX BIGC(NPTSMAX,0:2,2)
      DOUBLE COMPLEX NGAM(NPTSMAX,0:2,2,2)
      DOUBLE COMPLEX NGAMNS(NPTSMAX,0:2)
      DOUBLE COMPLEX BIGCF2(NPTSMAX,0:2,2)

*     - Final observables
      DOUBLE PRECISION F2(0:2)
      DOUBLE COMPLEX CFF(0:2)
        

*   Common blocks
*   -------------
      
*   1. Constant parameters loaded from 'GEPARD.INI'

      COMMON / PARINT /  SPEED, ACC, P, NF, CZERO
      COMMON / PARFLT /  RF2, RR2
      COMMON / MBCONT /  C, PHI, CND, PHIND
      COMMON / PARCHR /  SCHEME, ANSATZ

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
      COMMON / WGAMMA   /  GAM0, GAM1, GAM2
      COMMON / WGAMMANS /  GAMNS0, GAMNS1, GAMNS2
      COMMON / WC       /  C0, C1, C2

*     - Values on the whole contour
      COMMON / BIGC     /  BIGC
      COMMON / NGAM     /  NGAM
      COMMON / NGAMNS   /  NGAMNS
      COMMON / BIGCF2   /  BIGCF2

*     - FInal observables
      COMMON / CFF      /  CFF
      COMMON / F2       /  F2
