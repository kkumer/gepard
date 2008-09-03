C     ****h* gepard/header.f
C  FILE DESCRIPTION
C    Fortran header file with definition of common blocks
C    $Id:$
C     *******
    
#if defined(__GNUC__)
#define DBLCMPLX(X, Y)  COMPLEX(X, Y)
#elif defined(__INTEL_COMPILER) || defined(INTEL_COMPILER)
#define DBLCMPLX(X, Y)  CMPLX(X, Y, KIND(1.0d0))
#define REALPART(X) REAL(X)
#define IMAGPART(X) IMAG(X)
#endif

*   Parameters
*   ----------

*   1. Constant parameters loaded from 'GEPARD.INI' or fixed by program

C     ****v* header.f/SPEED
C  DESCRIPTION
C     SPEED -- Speed of evaluations: 1 (most accurate) - 3 (fastest).
C              This controls how fine is subdivision of Mellin-Barnes
C              integral.
C  MEMBER OF
C     PARINT
C  SYNOPSIS
      INTEGER SPEED
C     ***

C     ****v* header.f/ACC
C  DESCRIPTION
C     ACC -- Accuracy of numerical integrations: 1 (worst) - 3 (recommended) - 6 (best)
C            Mellin-Barnes integration is done on 2**ACC Gauss-Legendre
C            points per each subinterval.
C            Integration of partial cross-sections over DEL2 are done on
C            2**(ACC-SPEED+1) Gauss-Legendre points on interval (0, -1).
C  MEMBER OF
C     PARINT
C  SYNOPSIS
      INTEGER ACC
C     ***

C     ****v* header.f/P
C  DESCRIPTION
C         P -- approximation order, which is N^{P}LO, P=0, 1 or 2
C  MEMBER OF
C     PARINT
C  SYNOPSIS
      INTEGER P
C     ***

C     ****v* header.f/NF
C  DESCRIPTION
C        NF -- number of active light quark flavours
C  MEMBER OF
C     PARINT
C  SYNOPSIS
      INTEGER NF
C     ***

C     ****v* header.f/CZERO
C  DESCRIPTION
C     CZERO -- This multiplies $\alpha_{s}^{0}$ i. e. LO result.
C       Should normally be 1, but can be put to 0, so one can get NLO
C       correction, instead of prediction, from the routines.
C  MEMBER OF
C     PARINT
C  SYNOPSIS
      INTEGER CZERO
C     ***
      

C     ****v* header.f/MU02
C  DESCRIPTION
C       MU02 -- Starting scale for RGE evolution of alpha_strong.
C  MEMBER OF
C     ASTRONG
C  SYNOPSIS
      DOUBLE PRECISION MU02
C     ***

C     ****v* header.f/ASP
C  DESCRIPTION
C   ASP(P) -- Array holding the input values of alpha_strong/(2 pi) at the
C       scale MU02. This must be specified by the user, and for other
C       scales numerical integration of N^{P}LO order RGE equation will
C       be used (subroutine AS2PF).
C  MEMBER OF
C     ASTRONG
C  SYNOPSIS
      DOUBLE PRECISION ASP(0:2)
C     ***

C     ****v* header.f/Q02
C  DESCRIPTION
C       Q02 -- GPD input scale
C  MEMBER OF
C     PARFLT
C  SYNOPSIS
      DOUBLE PRECISION Q02
C     ***

C     ****v* header.f/RF2
C  DESCRIPTION
C        RF2 --  $ {\cal Q}^2/{\mu_{f}^2}$ i.e. ratio
C               of photon virtuality and factorization scales squared
C  MEMBER OF
C     PARFLT
C  SYNOPSIS
      DOUBLE PRECISION RF2
C     ***

C     ****v* header.f/RR2
C  DESCRIPTION
C        RR2 -- $ {\cal Q}^2/{\mu_{r}^2}$ i.e. ratio
C               of photon virtuality and renormalization scales squared
C  MEMBER OF
C     PARFLT
C  SYNOPSIS
      DOUBLE PRECISION RR2
C     ***

C     ****v* header.f/CHARGEFAC
C  DESCRIPTION
C     CHARGEFAC -- charge factor, cf. Eq. (79) hep-ph/0703179
C  MEMBER OF
C     PARFLT
C  SYNOPSIS
      DOUBLE PRECISION CHARGEFAC
C     ***

C     ****v* header.f/C
C  DESCRIPTION
C          C -- MB integration contour crosses real axis here
C  MEMBER OF
C     MBCONT
C  SYNOPSIS
      DOUBLE PRECISION C
C     ***

C     ****v* header.f/PHI
C  DESCRIPTION
C        PHI -- angle MB contour makes with positive real axis 
C  MEMBER OF
C     MBCONT
C  SYNOPSIS
      DOUBLE PRECISION PHI
C     ***

C     ****v* header.f/CND
C  DESCRIPTION
C        CND -- MB integration contour crosses real axis here
C               (contour for non-diagonal MS-bar evolution) 
C  MEMBER OF
C     MBCONT
C  SYNOPSIS
      DOUBLE PRECISION CND
C     ***

C     ****v* header.f/PHIND
C  DESCRIPTION
C      PHIND -- angle MB contour makes with positive real axis 
C               (contour for non-diagonal MS-bar evolution) 
C  MEMBER OF
C     MBCONT
C  SYNOPSIS
      DOUBLE PRECISION PHIND
C     ***

C     ****v* header.f/SCHEME
C  DESCRIPTION
C      SCHEME -- renormalization scheme: CSBAR, MSBAR, MSBND, MSBLO,
C                EVOLQ, EVOLG
C  MEMBER OF
C     PARCHR
C  SYNOPSIS
      CHARACTER SCHEME*5
C     ***

C     ****v* header.f/ANSATZ
C  DESCRIPTION
C      ANSATZ -- GPD ansatz: SOFT, HARD, FITBP, FIT, SPLICE, MMA, ...
C  MEMBER OF
C     PARCHR
C  SYNOPSIS
      CHARACTER ANSATZ*6
C     ***

C     ****v* header.f/PROCESS
C  DESCRIPTION
C      PROCESS -- process (DVCS or DIS)
C  MEMBER OF
C     PARCHR
C  SYNOPSIS
      CHARACTER PROCESS*6
C     ***

C     ****v* header.f/FFTYPE
C  DESCRIPTION
C      FFTYPE -- form-factor type (SINGLET or NONSINGLET)
C  MEMBER OF
C     PARCHR
C  SYNOPSIS
      CHARACTER FFTYPE*10
C     ***

C     ****v* header.f/DATFILE
C  DESCRIPTION
C      DATFILE -- file with specification of datasets
C                 used for fitting
C  MEMBER OF
C     FILENAMES
C  SYNOPSIS
      CHARACTER DATFILE*20
C     ***

C     ****v* header.f/OUTFILE
C  DESCRIPTION
C      OUTFILE -- file where output goes
C  MEMBER OF
C     FILENAMES
C  SYNOPSIS
      CHARACTER OUTFILE*20
C     ***

C     ****v* header.f/CMDFILE
C  DESCRIPTION
C      CMDFILE -- file with MINUIT batch commands
C  MEMBER OF
C     FILENAMES
C  SYNOPSIS
      CHARACTER CMDFILE*20
C     ***


*   2. Parameters from 'MINUIT.CMD' (candidates for fitting parameters)

C     ****v* header.f/PAR
C  DESCRIPTION
C      PAR -- When MINUIT calls chi-square function FCN values of
C             parameters are written to this array.
C             This array is in the same-name common block.
C  MEMBER OF
C     PAR(block)
C  SYNOPSIS
      INTEGER NPARMAX
      PARAMETER (NPARMAX = 50)
      DOUBLE PRECISION PAR(NPARMAX)
C     ***

*   3. Kinematics

C     ****v* header.f/XI
C  DESCRIPTION
C      XI --  Bjorken-like variable \xi for DVCS, or x_Bj for DIS
C  MEMBER OF
C     KINEMATICS
C  SYNOPSIS
      DOUBLE PRECISION XI
C     ***

C     ****v* header.f/DEL2
C  DESCRIPTION
C      DEL2 --  momentum transfer \Delta^2 = t
C  MEMBER OF
C     KINEMATICS
C  SYNOPSIS
      DOUBLE PRECISION DEL2
C     ***

C     ****v* header.f/Q2
C  DESCRIPTION
C      Q2  --  photon virtuality 
C  MEMBER OF
C     KINEMATICS
C  SYNOPSIS
      DOUBLE PRECISION Q2
C     ***

C     ****v* header.f/NQS
C  DESCRIPTION
C      NQS -- Number of different Q2 values occuring in DVCS
C             data. Set in subroutine INITC.
C  MEMBER OF
C     NQS(block)
C  SYNOPSIS
      INTEGER NQS
C     ***

C     ****v* header.f/NQSDIS
C  DESCRIPTION
C      NQSDIS -- Number of different Q2 values occuring in DIS
C                data. Set in subroutine INITC.
C  MEMBER OF
C     NQS(block)
C  SYNOPSIS
      INTEGER NQSDIS
C     ***

C     ****v* header.f/QS
C  DESCRIPTION
C      QS -- values of Q2 occuring in DVCS data. Set in INITC.
C  MEMBER OF
C     QS(block)
C  SYNOPSIS
      INTEGER QINDMAX
      PARAMETER (QINDMAX = 50)
      DOUBLE PRECISION QS(QINDMAX)
C     ***

C     ****v* header.f/QSDIS
C  DESCRIPTION
C      QSDIS -- values of Q2 occuring in DIS data. Set in INITC.
C  MEMBER OF
C     QS(block)
C  SYNOPSIS
      DOUBLE PRECISION  QSDIS(QINDMAX)
C     ***

C     ****v* header.f/MTIND
C  DESCRIPTION
C      MTIND -- index of position of actual -t=MT in MTS array
C  MEMBER OF
C      MT
C  SYNOPSIS
      INTEGER MTIND
C     ***

C     ****v* header.f/NMTS
C  DESCRIPTION
C      NMTS -- number of points for integration over -t.
C              Set in subroutine INIT.
C  MEMBER OF
C      MT
C  SYNOPSIS
      INTEGER NMTS
C     ***

C     ****v* header.f/NMTSEXP
C  DESCRIPTION
C      NMTSEXP -- number of -t values occuring in data.
C              Set by hand in subroutine INIT (unlike Q2
C              values which are read from data)!
C  MEMBER OF
C      MT
C  SYNOPSIS
      INTEGER NMTSEXP
C     ***

C     ****v* header.f/MTS
C  DESCRIPTION
C      MTS --  Values of -t=MT occuring in data or in
C              integration over -t:
C          MTS(0) = 0    (forward case, DIS)         
C          MTS(1..NMTS) = values for Gauss integration over -t
C          MTS(NMTS+1..NMTS+NMTSEXP) = values of -t in data
C  MEMBER OF
C     MTS(block)
C  SYNOPSIS
      INTEGER MTINDMAX
      PARAMETER (MTINDMAX = 100)
      DOUBLE PRECISION MTS(0:MTINDMAX)
C     ***


C     ****v* header.f/MTWG
C  DESCRIPTION
C      MTWG -- weights for Gauss integration over -t
C  MEMBER OF
C     MTS(block)
C  SYNOPSIS
      DOUBLE PRECISION MTWG(MTINDMAX)
C     ***


C     ****v* header.f/W2
C  DESCRIPTION
C      W2 --  c.m. energy squared
C  MEMBER OF
C     BCAKIN
C  SYNOPSIS
      DOUBLE PRECISION W2
C     ***

C     ****v* header.f/XB
C  DESCRIPTION
C      XB --  X Bjorken
C  MEMBER OF
C     BCAKIN
C  SYNOPSIS
      DOUBLE PRECISION XB
C     ***

C     ****v* header.f/YB
C  DESCRIPTION
C      YB --  "Y Bjorken": lepton energy fraction 
C  MEMBER OF
C     BCAKIN
C  SYNOPSIS
      DOUBLE PRECISION YB
C     ***

C     ****v* header.f/PHIAZ
C  DESCRIPTION
C      PHIAZ --  azymuthal angle between lepton plane and
C                incoming photon-scattered proton plane
C  MEMBER OF
C     BCAKIN
C  SYNOPSIS
      DOUBLE PRECISION PHIAZ
C     ***


*   4. Other

      DOUBLE PRECISION PI, PIHALF, CF, CA, TF, MP2
      PARAMETER ( PI = 3.1415 92653 58979 D0 )
      PARAMETER ( PIHALF = 1.5707 963267 948966 D0 )
      PARAMETER ( CF = 4.0d0 / 3.0d0, CA = 3.0d0, TF = 0.5d0 )
      PARAMETER ( MP2 = 0.880354d0 )
      DOUBLE COMPLEX TWOPII
      PARAMETER ( TWOPII = (0.0d0, 6.2831 853071 795865 D0) )

*     - QCD beta function
      INTEGER NFMIN, NFMAX
      PARAMETER (NFMIN = 3, NFMAX = 6)
      DOUBLE PRECISION BETA0, BETA1, BETA2, BETA3

*     - Mellin-Barnes integration contour points

C     ****v* header.f/NPTS
C  DESCRIPTION
C       NPTS  --  number of points on MB contour
C            Set in subroutine INIT.
C  MEMBER OF
C     CONTOUR
C  SYNOPSIS
      INTEGER NPTS
C     ***


C     ****v* header.f/NPTSMAX
C  DESCRIPTION
C       NPTSMAX  --  Maximal value of NPTS
C  SYNOPSIS
      INTEGER NPTSMAX
      PARAMETER (NPTSMAX = 768)
C     ***

C     ****v* header.f/Y
C  DESCRIPTION
C      Y(K) -- Real coordinate along MB contour
C            Set in subroutine INIT.
C  MEMBER OF
C     POINTS
C  SYNOPSIS
      DOUBLE PRECISION Y(NPTSMAX)
C     ***

C     ****v* header.f/WG
C  DESCRIPTION
C      WG(K) -- Weights for Gauss-Legendre integration over MB contour
C            Set in subroutine INIT.
C  MEMBER OF
C     POINTS
C  SYNOPSIS
      DOUBLE PRECISION  WG(NPTSMAX)
C     ***

C     ****v* header.f/N
C  DESCRIPTION
C      N(K)  --  Complex coordinates of points on MB contour 
C  MEMBER OF
C     NPOINTS
C  SYNOPSIS
      DOUBLE COMPLEX N(NPTSMAX)
C     ***

*     ----- Values on a particular contour point
      DOUBLE COMPLEX S1, S2, S3, S4

*     ----- Values on the whole contour

C     ****v* header.f/CDISNS1
C  DESCRIPTION
C     CDISNS1(K,P) -- DIS F1 nonsinglet Wilson coefficient  
C  MEMBER OF
C     CDISNS(block)
C  SYNOPSIS
      DOUBLE COMPLEX CDISNS1(NPTSMAX,0:2)
C     ***

C     ****v* header.f/CDIS1
C  DESCRIPTION
C     CDIS1(K,P,flavour) -- DIS F1 singlet Wilson coefficients
C  MEMBER OF
C     CDIS(block)
C  SYNOPSIS
      DOUBLE COMPLEX CDIS1(NPTSMAX,0:2,2)
C     ***

C     ****v* header.f/CDIS2
C  DESCRIPTION
C     CDIS2(1,P,flavour) -- DIS F2 singlet Wilson coefficients
C  MEMBER OF
C     CDIS(block)
C  SYNOPSIS
      DOUBLE COMPLEX CDIS2(NPTSMAX,0:2,2)
C     ***


C     ****v* header.f/BIGCNS
C  DESCRIPTION
C     BIGCNS(K,P) -- "Big C" DVCS nonsinglet Wilson coefficient  
C  MEMBER OF
C     BIGCNS(block)
C  SYNOPSIS
      DOUBLE COMPLEX BIGCNS(NPTSMAX,0:2)
C     ***


C     ****v* header.f/BIGC
C  DESCRIPTION
C     BIGC(K,P,flavour) -- "Big C" DVCS singlet Wilson coefficient  
C  MEMBER OF
C     BIGC(block)
C  SYNOPSIS
      DOUBLE COMPLEX BIGC(NPTSMAX,0:2,2)
C     ***

C     ****v* header.f/BIGCF2
C  DESCRIPTION
C     BIGCF2(K,P,flavour) -- "Big C" DIS F2 singlet Wilson coefficient  
C  MEMBER OF
C     BIGCF2(block)
C  SYNOPSIS
      DOUBLE COMPLEX BIGCF2(NPTSMAX,0:2,2)
C     ***


C     ****v* header.f/GAMNS
C  DESCRIPTION
C      GAMNS(K,P) -- non-singlet DIS anomalous dimensions 
C  MEMBER OF
C     GAMNS(block)
C  SYNOPSIS
      DOUBLE COMPLEX GAMNS(NPTSMAX,0:2)
C     ***

C     ****v* header.f/GAM
C  DESCRIPTION
C      GAM(K,P,flavour_i,flavour_j) -- non-singlet DIS anomalous dimensions 
C  MEMBER OF
C     GAMNS(block)
C  SYNOPSIS
      DOUBLE COMPLEX GAM(NPTSMAX,0:2,2,2)
C     ***

C     ****v* header.f/CGRIDNS
C  DESCRIPTION
C      CGRIDNS(QIND,K) -- "Evolved" "big C" nonsinglet Wilson
C      coefficients of DVCS i.e. multiplied by evolution operator for
C      evolution from Q02 to QS(QIND)
C  MEMBER OF
C    CGRIDNS(block)
C  SYNOPSIS
      DOUBLE COMPLEX CGRIDNS(QINDMAX, NPTSMAX)
C     ***


C     ****v* header.f/CGRID
C  DESCRIPTION
C      CGRID(QIND,K,flavour) -- "Evolved" "big C" singlet Wilson
C      coefficients of DVCS i.e. multiplied by evolution operator for
C      evolution from Q02 to QS(QIND)
C  MEMBER OF
C    CGRID(block)
C  SYNOPSIS
      DOUBLE COMPLEX CGRID(QINDMAX, NPTSMAX, 2)
C     ***

C     ****v* header.f/CGRIDDIS
C  DESCRIPTION
C      CGRIDDIS(QIND,K,flavour) -- "Evolved" "big C" singlet Wilson
C      coefficients of DIS F2 i.e. multiplied by evolution operator for
C      evolution from Q02 to QSDIS(QIND)
C  MEMBER OF
C    CGRIDDIS(block)
C  SYNOPSIS
      DOUBLE COMPLEX CGRIDDIS(QINDMAX, NPTSMAX, 2)
C     ***

C     ****v* header.f/HGRID
C  DESCRIPTION
C       HGRID(MTIND, K, flavour) -- values of GPD H
C               obsoleted by MBGPD, should be changed to MBGPD in code
C  MEMBER OF
C    HGRID(block) 
C  SYNOPSIS
      DOUBLE COMPLEX HGRID(0:MTINDMAX, NPTSMAX, 2)
C     ***

C     ****v* header.f/MBGPD
C  DESCRIPTION
C      MBGPD(K, flavour)  --  values of GPD H 
C  MEMBER OF
C    MBGPD(block) 
C  SYNOPSIS
      DOUBLE COMPLEX MBGPD(NPTSMAX, 2)
C     ***

*     - Final observables
C     ****v* header.f/F2
C  DESCRIPTION
C      F2(P)   --  DIS F2 function 
C  MEMBER OF
C    F2(block) 
C  SYNOPSIS
      DOUBLE PRECISION F2(0:2)
C     ***

C     ****v* header.f/CFF
C  DESCRIPTION
C      CFF(P)  --  DVCS Compton form factor  \mathcal{H}
C  MEMBER OF
C    CFF(block) 
C  SYNOPSIS
      DOUBLE COMPLEX CFF(0:2)
C     ***
        
C     ****v* header.f/CFFE
C  DESCRIPTION
C      CFFE(P)  --  DVCS Compton form factor \mathcal{E}
C  MEMBER OF
C    CFFE(block) 
C  SYNOPSIS
      DOUBLE COMPLEX CFFE(0:2)
C     ***
        
C     ****v* header.f/CHISQ
C  DESCRIPTION
C      CHISQ(0:10)   --  Total [CHISQ(0)] and partial chi-squares
C  MEMBER OF
C    CHISQBLK
C  SYNOPSIS
      DOUBLE PRECISION CHISQ(0:10)
C     ***

C     ****v* header.f/NDATAPTS
C  DESCRIPTION
C      NDATAPTS(0:10)   --  Number of data points belonging
C               to total [NDATAPTS(0)] and partial chi-squares
C  MEMBER OF
C    NDATAPTSBLK
C  SYNOPSIS
      INTEGER NDATAPTS(0:10)
C     ***


C     ****v* header.f/CHIN
C  DESCRIPTION
C      CHIN   --  Index of last partial chi-square
C  MEMBER OF
C    CHINBLK
C  SYNOPSIS
      INTEGER CHIN
C     ***


*   Common blocks
*   -------------
      
*   1. Constant parameters loaded from 'GEPARD.INI' or fixed by program


C     ****c* header.f/PARINT
C  DESCRIPTION
C     Integer parameters
C  SYNOPSIS
      COMMON / PARINT /  SPEED, ACC, P, NF, CZERO
C     ***

C     ****c* header.f/ASTRONG
C  DESCRIPTION 
C     Initial value of alpha_strong
C  SYNOPSIS
      COMMON / ASTRONG/  MU02, ASP
C     ***

C     ****c* header.f/PARFLT
C  DESCRIPTION
C     Floating point parameters
C  SYNOPSIS
      COMMON / PARFLT /  Q02, RF2, RR2, CHARGEFAC
C     ***

C     ****c* header.f/MBCONT
C  DESCRIPTION
C     Parameters specifying Mellin-Barnes contours
C  SYNOPSIS
      COMMON / MBCONT /  C, PHI, CND, PHIND
C     ***

C     ****c* header.f/PARCHR
C  DESCRIPTION
C     String parameters
C  SYNOPSIS
      COMMON / PARCHR /  SCHEME, ANSATZ, PROCESS, FFTYPE
C     ***

C     ****c* header.f/FILENAMES
C  DESCRIPTION
C     Filenames for input and output
C  SYNOPSIS
      COMMON / FILENAMES /  DATFILE, OUTFILE, CMDFILE
C     ***


*   2. Parameters from 'MINUIT.CMD' (candidates for fitting parameters)

C     ****c* header.f/PAR(block)
C  DESCRIPTION
C      Parameters from 'MINUIT.CMD' (candidates for fitting parameters)
C      This block holds only array of the same name.
C  SYNOPSIS
      COMMON / PAR    /  PAR
C     ***

*   3. Kinematics

C     ****c* header.f/KINEMATICS
C  DESCRIPTION
C      Kinematics of the process. Experimental parameters.
C  SYNOPSIS
      COMMON / KINEMATICS /  XI, DEL2, Q2
C     ***

C     ****c* header.f/NQS(block)
C  DESCRIPTION
C      Number of different Q2 occuring in the DVCS and DIS data
C  SYNOPSIS
      COMMON / NQS        /  NQS, NQSDIS
C     ***

C     ****c* header.f/QS(block)
C  DESCRIPTION
C      Values of Q2 occuring in the DVCS and DIS data
C  SYNOPSIS
      COMMON / QS         /  QS, QSDIS
C     ***

C     ****c* header.f/MT
C  DESCRIPTION
C      MT (i.e. -t) index and maximal values of this index
C  SYNOPSIS
      COMMON / MT         /  MTIND, NMTS, NMTSEXP
C     ***

C     ****c* header.f/MTS(block)
C  DESCRIPTION
C      Values of -t, and weights for integrals over -t
C  SYNOPSIS
      COMMON / MTS        /  MTS, MTWG
C     ***

C     ****c* header.f/BCAKIN
C  DESCRIPTION
C      Kinematics of the process, relevant for BCA
C  SYNOPSIS
      COMMON / BCAKIN     /  W2, XB,YB, PHIAZ
C     ***


*   4. Other


C     ****c* header.f/BETABLK
C  DESCRIPTION
C     QCD beta function coefficients
C  SYNOPSIS
      COMMON / BETABLK / BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX),
     &                   BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)
C     ***

*     ------- Mellin-Barnes integration contour points

C     ****c* header.f/CONTOUR
C  DESCRIPTION
C      Number of points on Mellin-Barnes contour. Actual points
C      are in common blocks POINTS (real coordinate along the
C      contour) and NPOINTS (complex coordinates).
C  SYNOPSIS
      COMMON / CONTOUR  /  NPTS
C     ***

C     ****c* header.f/POINTS
C  DESCRIPTION
C      Real coordinates along MB contour, and integration weights.
C      Corresponding complex coordinates are in common block NPOINTS.
C  SYNOPSIS
      COMMON / POINTS   /  Y, WG
C     ***

C     ****c* header.f/NPOINTS
C  DESCRIPTION
C      Complex coordinates of MB contour points.
C  SYNOPSIS
      COMMON / NPOINTS  /  N
C     ***
      
*     ----- Values on a particular contour point

C     ****c* header.f/HARMONIC
C  DESCRIPTION
C      Values of harmonic sums at actual MB point
C  SYNOPSIS
      COMMON / HARMONIC /  S1, S2, S3, S4
C     ***

*     ---- Values on the whole contour

C     ****c* header.f/CDISNS
C  DESCRIPTION
C     Nonsinglet DIS Wilson coefficients (small "c")  for F1
C  SYNOPSIS
      COMMON / CDISNS   /  CDISNS1
C     ***

C     ****c* header.f/CDIS
C  DESCRIPTION
C     Singlet DIS Wilson coefficients (small "c") for F1 and F2
C  SYNOPSIS
      COMMON / CDIS     /  CDIS1, CDIS2
C     ***

C     ****c* header.f/BIGCNS(block)
C  DESCRIPTION
C     Nonsinglet Wilson DIS coefficients  (big "C")
C  SYNOPSIS
      COMMON / BIGCNS   /  BIGCNS
C     ***

C     ****c* header.f/BIGC(block)
C  DESCRIPTION
C     Singlet Wilson DIS coefficients  (big "C")
C  SYNOPSIS
      COMMON / BIGC     /  BIGC
C     ***

C     ****c* header.f/GAMNS(block)
C  DESCRIPTION
C     Nonsinglet DIS anomalous dimensions 
C  SYNOPSIS
      COMMON / GAMNS    /  GAMNS
C     ***

C     ****c* header.f/GAM(block)
C  DESCRIPTION
C     Singlet DIS anomalous dimensions 
C  SYNOPSIS
      COMMON / GAM      /  GAM
C     ***

C     ****c* header.f/BIGCF2(block)
C  DESCRIPTION
C     DIS F2 Wilson coefficient 
C  SYNOPSIS
      COMMON / BIGCF2   /  BIGCF2
C     ***

C     ****c* header.f/CGRIDNS(block)
C  DESCRIPTION
C     "Evolved" non-singlet DVCS Wilson coefficients
C  SYNOPSIS
      COMMON / CGRIDNS  /  CGRIDNS
C     ***

C     ****c* header.f/CGRID(block)
C  DESCRIPTION
C     "Evolved" singlet DVCS Wilson coefficients
C  SYNOPSIS
      COMMON / CGRID    /  CGRID
C     ***

C     ****c* header.f/CGRIDDIS(block)
C  DESCRIPTION
C     "Evolved" DIS F2 Wilson coefficients
C  SYNOPSIS
      COMMON / CGRIDDIS /  CGRIDDIS
C     ***

C     ****c* header.f/HGRID(block)
C  DESCRIPTION
C      Values of GPDs on MB contour for all needed -t (in principle)
C      Maybe obsoleted by MBGPD
C  SYNOPSIS
      COMMON / HGRID    /  HGRID
C     ***

C     ****c* header.f/MBGPD(block)
C  DESCRIPTION
C     Values of GPDs on MB contour
C  SYNOPSIS
      COMMON / MBGPD    /  MBGPD
C     ***

*     ---- Final observables

C     ****c* header.f/CFF(block)
C  DESCRIPTION
C     Compton form factor \mathcal{H}
C  SYNOPSIS
      COMMON / CFF      /  CFF, CFFE
C     ***

C     ****c* header.f/F2(block)
C  DESCRIPTION
C     DIS F2 form factor
C  SYNOPSIS
      COMMON / F2       /  F2
C     ***

C     ****c* header.f/CHISQBLK
C  DESCRIPTION
C      Values of chi-square
C  SYNOPSIS
      COMMON / CHISQBLK /  CHISQ
C     ***

C     ****c* header.f/NDATAPTSBLK
C  DESCRIPTION
C      Number of data points used for fit
C  SYNOPSIS
      COMMON / NDATAPTSBLK /  NDATAPTS
C     ***


C     ****c* header.f/CHINBLK
C  DESCRIPTION
C      Index of last partial chi-square in CHISQ
C  SYNOPSIS
      COMMON / CHINBLK /  CHIN
C     ***


