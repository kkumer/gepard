C     ****h* gepard/header.F
C  FILE DESCRIPTION
C    Fortran header file with definition of common blocks
C    $Id:$
C     *******
    
#if defined(__GNUC__)
#define DBLCMPLX(X, Y)  COMPLEX(X, Y)
#define DCTAN  G77TAN
#elif defined(__INTEL_COMPILER) || defined(INTEL_COMPILER)
#define DBLCMPLX(X, Y)  CMPLX(X, Y, KIND(1.0d0))
#define REALPART(X) REAL(X)
#define IMAGPART(X) IMAG(X)
#define DCTAN  TAN
#elif defined(NAGf90Fortran)
#define DBLCMPLX(X, Y)  CMPLX(X, Y, KIND(1.0d0))
#define REALPART(X) REAL(X)
#define IMAGPART(X) AIMAG(X)
#define DCTAN  G77TAN
#define IARGC COMMAND_ARGUMENT_COUNT
#define GETARG GET_COMMAND_ARGUMENT
#endif

*   Parameters
*   ----------

*   1. Constant parameters loaded from 'GEPARD.INI' or fixed by program

C     ****v* header.F/SPEED
C  DESCRIPTION
C     SPEED -- Speed of evaluations: 1 (most accurate) - 3 (fastest).
C              This controls how fine is subdivision of Mellin-Barnes
C              integral.
C  MEMBER OF
C     PARINT
C  SYNOPSIS
      INTEGER SPEED
C     ***

C     ****v* header.F/ACC
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

C     ****v* header.F/P
C  DESCRIPTION
C         P -- approximation order, which is N^{P}LO, P=0, 1 or 2
C  MEMBER OF
C     PARINT
C  SYNOPSIS
      INTEGER P
C     ***

C     ****v* header.F/NF
C  DESCRIPTION
C        NF -- number of active light quark flavours
C  MEMBER OF
C     PARINT
C  SYNOPSIS
      INTEGER NF
C     ***

C     ****v* header.F/CZERO
C  DESCRIPTION
C     CZERO -- This multiplies $\alpha_{s}^{0}$ i. e. LO result.
C       Should normally be 1, but can be put to 0, so one can get NLO
C       correction, instead of prediction, from the routines.
C  MEMBER OF
C     PARINT
C  SYNOPSIS
      INTEGER CZERO
C     ***
      
C     ****v* header.F/PID
C  DESCRIPTION
C       PID -- process ID. with possible values:
C               0 = DIS F2
C               1 = DVCS (CSBAR)
C               2 = DVCS (MSBAR)
C               3 = DVMP (FIXME: just rho presently)
C              -1 = quark GPD at eta = 0
C              -2 = gluon GPD at eta = 0
C              -3 = quark GPD at eta = x
C              -4 = gluon GPD at eta = x
C  MEMBER OF
C     PARINT
C  SYNOPSIS
      INTEGER PID
C     ***
      


C     ****v* header.F/MU02
C  DESCRIPTION
C       MU02 -- Starting scale for RGE evolution of alpha_strong.
C  MEMBER OF
C     ASTRONG
C  SYNOPSIS
      DOUBLE PRECISION MU02
C     ***

C     ****v* header.F/ASP
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

C     ****v* header.F/Q02
C  DESCRIPTION
C       Q02 -- GPD input scale
C  MEMBER OF
C     PARFLT
C  SYNOPSIS
      DOUBLE PRECISION Q02
C     ***

C     ****v* header.F/RF2
C  DESCRIPTION
C        RF2 --  $ {\cal Q}^2/{\mu_{f}^2}$ i.e. ratio
C               of photon virtuality and factorization scales squared
C               (for DVCS)
C  MEMBER OF
C     PARFLT
C  SYNOPSIS
      DOUBLE PRECISION RF2
C     ***

C     ****v* header.F/RGPDF2
C  DESCRIPTION
C        RGPDF2 --  $ {\cal Q}^2/{\mu_{f,GPD}^2}$ i.e. ratio
C               of photon virtuality and GPD factorization scales squared
C               (for DVEM)
C  MEMBER OF
C     PARFLT
C  SYNOPSIS
      DOUBLE PRECISION RGPDF2
C     ***

C     ****v* header.F/RDAF2
C  DESCRIPTION
C        RDAF2 --  $ {\cal Q}^2/{\mu_{f,DA}^2}$ i.e. ratio
C               of photon virtuality and DA factorization scales squared
C               (for DVEM)
C  MEMBER OF
C     PARFLT
C  SYNOPSIS
      DOUBLE PRECISION RDAF2
C     ***


C     ****v* header.F/RR2
C  DESCRIPTION
C        RR2 -- $ {\cal Q}^2/{\mu_{r}^2}$ i.e. ratio
C               of photon virtuality and renormalization scales squared
C  MEMBER OF
C     PARFLT
C  SYNOPSIS
      DOUBLE PRECISION RR2
C     ***

C     ****v* header.F/CHARGEFAC
C  DESCRIPTION
C     CHARGEFAC -- charge factor, cf. Eq. (79) hep-ph/0703179
C  MEMBER OF
C     PARFLT
C  SYNOPSIS
      DOUBLE PRECISION CHARGEFAC
C     ***

C     ****v* header.F/C
C  DESCRIPTION
C          C -- MB integration contour crosses real axis here
C  MEMBER OF
C     MBCONT
C  SYNOPSIS
      DOUBLE PRECISION C
C     ***

C     ****v* header.F/PHI
C  DESCRIPTION
C        PHI -- angle MB contour makes with positive real axis 
C  MEMBER OF
C     MBCONT
C  SYNOPSIS
      DOUBLE PRECISION PHI
C     ***

C     ****v* header.F/CND
C  DESCRIPTION
C        CND -- MB integration contour crosses real axis here
C               (contour for non-diagonal MS-bar evolution) 
C  MEMBER OF
C     MBCONT
C  SYNOPSIS
      DOUBLE PRECISION CND
C     ***

C     ****v* header.F/PHIND
C  DESCRIPTION
C      PHIND -- angle MB contour makes with positive real axis 
C               (contour for non-diagonal MS-bar evolution) 
C  MEMBER OF
C     MBCONT
C  SYNOPSIS
      DOUBLE PRECISION PHIND
C     ***

C     ****v* header.F/SCHEME
C  DESCRIPTION
C      SCHEME -- DVCS renormalization scheme: CSBAR, MSBAR
C  MEMBER OF
C     PARCHR
C  SYNOPSIS
      CHARACTER SCHEME*5
C     ***

C     ****v* header.F/ANSATZ
C  DESCRIPTION
C      ANSATZ -- GPD ansatz: SOFT, HARD, FITBP, FIT, SPLICE, MMA, ...
C  MEMBER OF
C     PARCHR
C  SYNOPSIS
      CHARACTER ANSATZ*6
C     ***

C     ****v* header.F/FFTYPE
C  DESCRIPTION
C      FFTYPE -- form-factor type (SINGLET or NONSINGLET)
C  MEMBER OF
C     PARCHR
C  SYNOPSIS
      CHARACTER FFTYPE*10
C     ***

C     ****v* header.F/DATFILE
C  DESCRIPTION
C      DATFILE -- file with specification of datasets
C                 used for fitting
C  MEMBER OF
C     FILENAMES
C  SYNOPSIS
      CHARACTER DATFILE*20
C     ***

C     ****v* header.F/OUTFILE
C  DESCRIPTION
C      OUTFILE -- file where output goes
C  MEMBER OF
C     FILENAMES
C  SYNOPSIS
      CHARACTER OUTFILE*20
C     ***

C     ****v* header.F/CMDFILE
C  DESCRIPTION
C      CMDFILE -- file with MINUIT batch commands
C  MEMBER OF
C     FILENAMES
C  SYNOPSIS
      CHARACTER CMDFILE*20
C     ***


*   2. Parameters from 'MINUIT.CMD' (candidates for fitting parameters)

C     ****v* header.F/PAR
C  DESCRIPTION
C      PAR -- When MINUIT calls chi-square function FCN values of
C             parameters are written to this array.
C             This array is in the same-name common block.
C  MEMBER OF
C     PAR(block)
C  SYNOPSIS
      INTEGER NPARMAX
      PARAMETER (NPARMAX = 150)
      DOUBLE PRECISION PAR(NPARMAX)
C     ***

*   3. Kinematics

C     ****v* header.F/XI
C  DESCRIPTION
C      XI --  Bjorken-like variable \xi for DVCS, or x_Bj for DIS
C  MEMBER OF
C     KINEMATICS
C  SYNOPSIS
      DOUBLE PRECISION XI
C     ***

C     ****v* header.F/DEL2
C  DESCRIPTION
C      DEL2 --  momentum transfer \Delta^2 = t
C  MEMBER OF
C     KINEMATICS
C  SYNOPSIS
      DOUBLE PRECISION DEL2
C     ***

C     ****v* header.F/Q2
C  DESCRIPTION
C      Q2  --  photon virtuality 
C  MEMBER OF
C     KINEMATICS
C  SYNOPSIS
      DOUBLE PRECISION Q2
C     ***

C     ****v* header.F/NQS
C  DESCRIPTION
C      NQS -- Number of different Q2 values used for given PID
C  MEMBER OF
C     NQS(block)
C  SYNOPSIS
      INTEGER PIDMAX
      PARAMETER (PIDMAX = 3)
      INTEGER NQS(-4:PIDMAX)
C     ***

C     ****v* header.F/QS
C  DESCRIPTION
C      QS -- values of Q2 occuring in DVCS data. Set in INITC.
C  MEMBER OF
C     QS(block)
C  SYNOPSIS
      INTEGER QINDMAX
      PARAMETER (QINDMAX = 150)
      DOUBLE PRECISION QS(-4:PIDMAX, QINDMAX)
C     ***


C     ****v* header.F/MTIND
C  DESCRIPTION
C      MTIND -- index of position of actual -t=MT in MTS array
C  MEMBER OF
C      MT
C  SYNOPSIS
      INTEGER MTIND
C     ***

C     ****v* header.F/NMTS
C  DESCRIPTION
C      NMTS -- number of points for integration over -t.
C              Set in subroutine INIT.
C  MEMBER OF
C      MT
C  SYNOPSIS
      INTEGER NMTS
C     ***

C     ****v* header.F/NMTSEXP
C  DESCRIPTION
C      NMTSEXP -- number of -t values occuring in data.
C              Set by hand in subroutine INIT (unlike Q2
C              values which are read from data)!
C  MEMBER OF
C      MT
C  SYNOPSIS
      INTEGER NMTSEXP
C     ***

C     ****v* header.F/MTS
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


C     ****v* header.F/MTWG
C  DESCRIPTION
C      MTWG -- weights for Gauss integration over -t
C  MEMBER OF
C     MTS(block)
C  SYNOPSIS
      DOUBLE PRECISION MTWG(MTINDMAX)
C     ***


C     ****v* header.F/W2
C  DESCRIPTION
C      W2 --  c.m. energy squared
C  MEMBER OF
C     BCAKIN
C  SYNOPSIS
      DOUBLE PRECISION W2
C     ***

C     ****v* header.F/XB
C  DESCRIPTION
C      XB --  X Bjorken
C  MEMBER OF
C     BCAKIN
C  SYNOPSIS
      DOUBLE PRECISION XB
C     ***

C     ****v* header.F/YB
C  DESCRIPTION
C      YB --  "Y Bjorken": lepton energy fraction 
C  MEMBER OF
C     BCAKIN
C  SYNOPSIS
      DOUBLE PRECISION YB
C     ***

C     ****v* header.F/PHIAZ
C  DESCRIPTION
C      PHIAZ --  azymuthal angle between lepton plane and
C                incoming photon-scattered proton plane
C  MEMBER OF
C     BCAKIN
C  SYNOPSIS
      DOUBLE PRECISION PHIAZ
C     ***


*   4. Other

      DOUBLE PRECISION PI, PIHALF, CF, CA, CG, TF, MP2, SQRT2
      DOUBLE PRECISION EMC, ZETA2, ZETA3, ZETA4
      PARAMETER ( EMC  = 0.5772156649D0, ZETA2 = 1.6449340668D0,
     ,            ZETA3 = 1.2020569032D0, ZETA4 = 1.0823232337D0 )
      PARAMETER ( PI = 3.14159265358979D0 )
      PARAMETER ( PIHALF = 1.5707963267948966D0 )
      PARAMETER ( CF = 4.0d0/3.0d0, CA = 3.0d0, TF = 0.5d0 )
      PARAMETER ( CG = -1.0d0/6.0d0 )
      PARAMETER ( MP2 = 0.880354d0 )
      PARAMETER ( SQRT2 = 1.414213562373095d0 )
      DOUBLE COMPLEX TWOPII
      PARAMETER ( TWOPII = (0.0d0, 6.2831853071795865D0) )

*     - Meson decay constants
      DOUBLE PRECISION FRHO, FOMEGA, FPHI
      PARAMETER (FRHO = 0.209, FOMEGA = 0.195, FPHI = 0.221)

*     - QCD beta function
      INTEGER NFMIN, NFMAX, NC
      PARAMETER (NFMIN = 3, NFMAX = 6, NC = 3)
      DOUBLE PRECISION BETA0, BETA1, BETA2, BETA3

*     - Mellin-Barnes integration contour points

C     ****v* header.F/NPTS
C  DESCRIPTION
C       NPTS  --  number of points on MB contour
C            Set in subroutine INIT.
C  MEMBER OF
C     CONTOUR
C  SYNOPSIS
      INTEGER NPTS
C     ***


C     ****v* header.F/NPTSMAX
C  DESCRIPTION
C       NPTSMAX  --  Maximal value of NPTS
C  SYNOPSIS
      INTEGER NPTSMAX
      PARAMETER (NPTSMAX = 768)
C     ***

C     ****v* header.F/Y
C  DESCRIPTION
C      Y(K) -- Real coordinate along MB contour
C            Set in subroutine INIT.
C  MEMBER OF
C     POINTS
C  SYNOPSIS
      DOUBLE PRECISION Y(NPTSMAX)
C     ***

C     ****v* header.F/WG
C  DESCRIPTION
C      WG(K) -- Weights for Gauss-Legendre integration over MB contour
C            Set in subroutine INIT.
C  MEMBER OF
C     POINTS
C  SYNOPSIS
      DOUBLE PRECISION  WG(NPTSMAX)
C     ***

C     ****v* header.F/N
C  DESCRIPTION
C      N(K)  --  Complex coordinates of points on MB contour 
C  MEMBER OF
C     NPOINTS
C  SYNOPSIS
      INTEGER NPWMAX
      PARAMETER (NPWMAX = 3)
      DOUBLE COMPLEX N(NPWMAX,NPTSMAX)
C     ***

*     ----- Values on a particular contour point
      DOUBLE COMPLEX S1, S2, S3, S4

*     ----- Values on the whole contour

C     ****v* header.F/CDISNS1
C  DESCRIPTION
C     CDISNS1(K,P) -- DIS F1 nonsinglet Wilson coefficient  
C  MEMBER OF
C     CDISNS(block)
C  SYNOPSIS
      DOUBLE COMPLEX CDISNS1(NPTSMAX,0:2)
C     ***

C     ****v* header.F/CDIS1
C  DESCRIPTION
C     CDIS1(K,P,flavour) -- DIS F1 singlet Wilson coefficients
C  MEMBER OF
C     CDIS(block)
C  SYNOPSIS
      DOUBLE COMPLEX CDIS1(NPTSMAX,0:2,2)
C     ***

C     ****v* header.F/CDIS2
C  DESCRIPTION
C     CDIS2(1,P,flavour) -- DIS F2 singlet Wilson coefficients
C  MEMBER OF
C     CDIS(block)
C  SYNOPSIS
      DOUBLE COMPLEX CDIS2(NPTSMAX,0:2,2)
C     ***


C     ****v* header.F/BIGCNS
C  DESCRIPTION
C     BIGCNS(K,P) -- "Big C" DVCS nonsinglet Wilson coefficient  
C  MEMBER OF
C     BIGCNS(block)
C  SYNOPSIS
      DOUBLE COMPLEX BIGCNS(NPTSMAX,0:2)
C     ***


C     ****v* header.F/BIGC
C  DESCRIPTION
C     BIGC(SEC,K,P,flavour) -- "Big C" DVCS singlet Wilson coefficient  
C  MEMBER OF
C     BIGC(block)
C  SYNOPSIS
      DOUBLE COMPLEX BIGC(NPWMAX,NPTSMAX,0:2,2)
C     ***

C     ****v* header.F/WC
C  DESCRIPTION
C     WC(PID, SEC, K, P, L) -- all Wilson coefficients
C  MEMBER OF
C     WC(block)
C  SYNOPSIS
      DOUBLE COMPLEX WC(-4:PIDMAX,NPWMAX,NPTSMAX,0:2,4)
C     ***


C     ****v* header.F/BIGCF2
C  DESCRIPTION
C     BIGCF2(K,P,flavour) -- "Big C" DIS F2 singlet Wilson coefficient  
C  MEMBER OF
C     BIGCF2(block)
C  SYNOPSIS
      DOUBLE COMPLEX BIGCF2(NPTSMAX,0:2,2)
C     ***


C     ****v* header.F/GAMNS
C  DESCRIPTION
C      GAMNS(SEC, K, P, 0) -- non-singlet DIS anomalous dimensions  \gamma^+
C      GAMNS(SEC, K, P, 1) -- non-singlet DIS anomalous dimensions  \gamma^-
C  MEMBER OF
C     GAM(block)
C  SYNOPSIS
      DOUBLE COMPLEX GAMNS(NPWMAX,NPTSMAX,0:2,0:1)
C     ***

C     ****v* header.F/GAM
C  DESCRIPTION
C      GAM(SEC,K,P,flavour_i,flavour_j) -- singlet DIS anomalous dimensions 
C  MEMBER OF
C     GAM(block)
C  SYNOPSIS
      DOUBLE COMPLEX GAM(NPWMAX,NPTSMAX,0:2,2,2)
C     ***

C     ****v* header.F/WCE
C  DESCRIPTION
C      WCE(PID, SEC, QIND, K, L) -- "Evolved"  Wilson
C      coefficients for all processes i.e. multiplied by evolution operator for
C      evolution from Q02 to QS(QIND)
C  MEMBER OF
C    WCE(block)
C  SYNOPSIS
      DOUBLE COMPLEX WCE(-4:PIDMAX, NPWMAX, QINDMAX, NPTSMAX, 4)
C     ***


C     ****v* header.F/GPDH
C  DESCRIPTION
C      GPDH(K, index)  --  values of conformal moments of GPD H
C          index: 1  H_S
C          index: 2  H_G
C          index: 3  H_NS
C          index: 4  H_NS 
C          index: 5  H_NS ..
C      where NS = u, d, s  or (3), (8), (15), or whatever ...
C  MEMBER OF
C    MBGPD(block) 
C  SYNOPSIS
      DOUBLE COMPLEX GPDH(NPTSMAX, 5)
C     ***

C     ****v* header.F/GPDE
C  DESCRIPTION
C      GPDE(K, index)  --  values of conformal moments of GPD E
C          index: [see GPDH]
C  MEMBER OF
C    MBGPD(block) 
C  SYNOPSIS
      DOUBLE COMPLEX GPDE(NPTSMAX, 5)
C     ***


*     - Final observables
C     ****v* header.F/F2
C  DESCRIPTION
C      F2(P)   --  DIS F2 function 
C  MEMBER OF
C    F2(block) 
C  SYNOPSIS
      DOUBLE PRECISION F2(0:2)
C     ***

C     ****v* header.F/CFF
C  DESCRIPTION
C      CFF(P)  --  DVCS Compton form factor  \mathcal{H}
C  MEMBER OF
C    CFF(block) 
C  SYNOPSIS
      DOUBLE COMPLEX CFF(0:2)
C     ***
        
C     ****v* header.F/CFFE
C  DESCRIPTION
C      CFFE(P)  --  DVCS Compton form factor \mathcal{E}
C  MEMBER OF
C    CFF(block) 
C  SYNOPSIS
      DOUBLE COMPLEX CFFE(0:2)
C     ***

C     ****v* header.F/TFFH
C  DESCRIPTION
C      TFFH(P)  --  DVMP transiton form factor  \mathcal{H}_M
C  MEMBER OF
C    TFF(block)
C  SYNOPSIS
      DOUBLE COMPLEX TFFH(0:2)
C     ***
        
        
C     ****v* header.F/CHISQ
C  DESCRIPTION
C      CHISQ(0:10)   --  Total [CHISQ(0)] and partial chi-squares
C  MEMBER OF
C    CHISQBLK
C  SYNOPSIS
      DOUBLE PRECISION CHISQ(0:10)
C     ***

C     ****v* header.F/NDATAPTS
C  DESCRIPTION
C      NDATAPTS(0:10)   --  Number of data points belonging
C               to total [NDATAPTS(0)] and partial chi-squares
C  MEMBER OF
C    NDATAPTSBLK
C  SYNOPSIS
      INTEGER NDATAPTS(0:10)
C     ***


C     ****v* header.F/CHIN
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


C     ****c* header.F/PARINT
C  DESCRIPTION
C     Integer parameters
C  SYNOPSIS
      COMMON / PARINT /  SPEED, ACC, P, NF, CZERO, PID
C     ***

C     ****c* header.F/ASTRONG
C  DESCRIPTION 
C     Initial value of alpha_strong
C  SYNOPSIS
      COMMON / ASTRONG/  MU02, ASP
C     ***

C     ****c* header.F/PARFLT
C  DESCRIPTION
C     Floating point parameters
C  SYNOPSIS
      COMMON / PARFLT /  Q02, RF2, RGPDF2, RDAF2, RR2, CHARGEFAC
C     ***

C     ****c* header.F/MBCONT
C  DESCRIPTION
C     Parameters specifying Mellin-Barnes contours
C  SYNOPSIS
      COMMON / MBCONT /  C, PHI, CND, PHIND
C     ***

C     ****c* header.F/PARCHR
C  DESCRIPTION
C     String parameters
C  SYNOPSIS
      COMMON / PARCHR /  SCHEME, ANSATZ, FFTYPE
C     ***

C     ****c* header.F/FILENAMES
C  DESCRIPTION
C     Filenames for input and output
C  SYNOPSIS
      COMMON / FILENAMES /  DATFILE, OUTFILE, CMDFILE
C     ***


*   2. Parameters from 'MINUIT.CMD' (candidates for fitting parameters)

C     ****c* header.F/PAR(block)
C  DESCRIPTION
C      Parameters from 'MINUIT.CMD' (candidates for fitting parameters)
C      This block holds only array of the same name.
C  SYNOPSIS
      COMMON / PAR    /  PAR
C     ***

*   3. Kinematics

C     ****c* header.F/KINEMATICS
C  DESCRIPTION
C      Kinematics of the process. Experimental parameters.
C  SYNOPSIS
      COMMON / KINEMATICS /  XI, DEL2, Q2
C     ***

C     ****c* header.F/NQS(block)
C  DESCRIPTION
C      Number of different Q2 occuring in the DVCS and DIS data
C  SYNOPSIS
      COMMON / NQS        /  NQS
C     ***

C     ****c* header.F/QS(block)
C  DESCRIPTION
C      Values of Q2 occuring in the DVCS and DIS data
C  SYNOPSIS
      COMMON / QS         /  QS
C     ***

C     ****c* header.F/MT
C  DESCRIPTION
C      MT (i.e. -t) index and maximal values of this index
C  SYNOPSIS
      COMMON / MT         /  MTIND, NMTS, NMTSEXP
C     ***

C     ****c* header.F/MTS(block)
C  DESCRIPTION
C      Values of -t, and weights for integrals over -t
C  SYNOPSIS
      COMMON / MTS        /  MTS, MTWG
C     ***

C     ****c* header.F/BCAKIN
C  DESCRIPTION
C      Kinematics of the process, relevant for BCA
C  SYNOPSIS
      COMMON / BCAKIN     /  W2, XB,YB, PHIAZ
C     ***


*   4. Other


C     ****c* header.F/BETABLK
C  DESCRIPTION
C     QCD beta function coefficients
C  SYNOPSIS
      COMMON / BETABLK / BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX),
     &                   BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)
C     ***

*     ------- Mellin-Barnes integration contour points

C     ****c* header.F/CONTOUR
C  DESCRIPTION
C      Number of points on Mellin-Barnes contour. Actual points
C      are in common blocks POINTS (real coordinate along the
C      contour) and NPOINTS (complex coordinates).
C  SYNOPSIS
      COMMON / CONTOUR  /  NPTS
C     ***

C     ****c* header.F/POINTS
C  DESCRIPTION
C      Real coordinates along MB contour, and integration weights.
C      Corresponding complex coordinates are in common block NPOINTS.
C  SYNOPSIS
      COMMON / POINTS   /  Y, WG
C     ***

C     ****c* header.F/NPOINTS
C  DESCRIPTION
C      Complex coordinates of MB contour points.
C  SYNOPSIS
      COMMON / NPOINTS  /  N
C     ***
      
*     ----- Values on a particular contour point

C     ****c* header.F/HARMONIC
C  DESCRIPTION
C      Values of harmonic sums at actual MB point
C  SYNOPSIS
      COMMON / HARMONIC /  S1, S2, S3, S4
C     ***

*     ---- Values on the whole contour

C     ****c* header.F/CDISNS
C  DESCRIPTION
C     Nonsinglet DIS Wilson coefficients (small "c")  for F1
C  SYNOPSIS
      COMMON / CDISNS   /  CDISNS1
C     ***

C     ****c* header.F/CDIS
C  DESCRIPTION
C     Singlet DIS Wilson coefficients (small "c") for F1 and F2
C  SYNOPSIS
      COMMON / CDIS     /  CDIS1, CDIS2
C     ***

C     ****c* header.F/BIGCNS(block)
C  DESCRIPTION
C     Nonsinglet Wilson DIS coefficients  (big "C")
C  SYNOPSIS
      COMMON / BIGCNS   /  BIGCNS
C     ***

C     ****c* header.F/BIGC(block)
C  DESCRIPTION
C     Singlet Wilson DIS coefficients  (big "C")
C  SYNOPSIS
      COMMON / BIGC     /  BIGC
C     ***

C     ****c* header.F/WC(block)
C  DESCRIPTION
C     All Wilson coefficients
C  SYNOPSIS
      COMMON / WC       /  WC
C     ***

C     ****c* header.F/GAM(block)
C  DESCRIPTION
C     DIS anomalous dimensions 
C  SYNOPSIS
      COMMON / GAM      /  GAM, GAMNS
C     ***

C     ****c* header.F/BIGCF2(block)
C  DESCRIPTION
C     DIS F2 Wilson coefficient 
C  SYNOPSIS
      COMMON / BIGCF2   /  BIGCF2
C     ***


C     ****c* header.F/WCE(block)
C  DESCRIPTION
C     "Evolved" Wilson coefficients
C  SYNOPSIS
      COMMON / WCE      /  WCE
C     ***


C     ****c* header.F/MBGPD(block)
C  DESCRIPTION
C     Values of GPDs on MB contour
C  SYNOPSIS
      COMMON / MBGPD    /  GPDH, GPDE
C     ***

*     ---- Final observables

C     ****c* header.F/CFF(block)
C  DESCRIPTION
C     Compton form factor \mathcal{H}, \mathcal{E}
C  SYNOPSIS
      COMMON / CFF      /  CFF, CFFE
C     ***

C     ****c* header.F/TFF(block)
C  DESCRIPTION
C     DVMP transition form factors \mathcal{H}_M
C  SYNOPSIS
      COMMON / TFF      /  TFFH
C     ***

C     ****c* header.F/F2(block)
C  DESCRIPTION
C     DIS F2 form factor
C  SYNOPSIS
      COMMON / F2       /  F2
C     ***

C     ****c* header.F/CHISQBLK
C  DESCRIPTION
C      Values of chi-square
C  SYNOPSIS
      COMMON / CHISQBLK /  CHISQ
C     ***

C     ****c* header.F/NDATAPTSBLK
C  DESCRIPTION
C      Number of data points used for fit
C  SYNOPSIS
      COMMON / NDATAPTSBLK /  NDATAPTS
C     ***


C     ****c* header.F/CHINBLK
C  DESCRIPTION
C      Index of last partial chi-square in CHISQ
C  SYNOPSIS
      COMMON / CHINBLK /  CHIN
C     ***


