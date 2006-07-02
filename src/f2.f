C     ****h* gepard/f2.f
C  FILE DESCRIPTION
C    calculation of F2 DIS form factor
C
C    $Id$
C     *******


C     ****s* f2.f/F2F
C  NAME
C       F2F  --   singlet F2 DIS form factor
C  DESCRIPTION
C    calculates singlet F2 DIS form factor by numerical Mellin-Barnes 
C    integration of conformal partial waves provided by subroutine PARWAV. 
C    |latex  \begin{eqnarray*}
C    |latex  {{\cal F}^{\rm S}}(\xi,{\cal Q}^2)
C    |latex  &\!\!\!=\!\!\!& \frac{Q_{S}^2}{2 \pi i}\int_{c-i \infty}^{c+ i \infty}\!
C    |latex  dj\,\xi^{-j} 
C    |latex  \mbox{\boldmath $C$}^{\rm DIS}_{j}(Q^2/\mu^2,\alpha_s(\mu)) 
C    |latex  \mbox{\boldmath $H$}_{j} (\xi=0,\Delta^2=0,\mu^2) \\[2ex]
C    |latex  &\!\!\!=\!\!\!& \frac{Q_{S}^2}{\pi} \xi^{-c}  
C    |latex  \int_{0}^{\infty} dy\: \Im\! m\: e^{i\phi} \xi^{- y e^{i\phi}} 
C    |latex  \texttt{FPW}(j) 
C    |latex  \end{eqnarray*} 
C    |latex  where $j=c+ye^{i\phi}$.
C    |latex  This is implemented now with fixed $\phi = \pi/2$, but
C    |latex  implementation with, say, $\phi = 3\pi /4$ should give faster
C    |latex  convergence.
C    Unlike CFF, this is defined with the singlet charge factor 2/9 !
C    All INPUTS and OUTPUT below is via common blocks!
C  SYNOPSIS
C     SUBROUTINE F2F
C
C     INTEGER P
C     DOUBLE PRECISION XI, DEL2, Q2, Q02, C
C     DOUBLE PRECISION F2(0:2)
C     CHARACTER SCHEME*5, ANSATZ*6
C
C     COMMON / KINEMATICS /  XI, DEL2, Q2, Q02
C     COMMON / APPROX     /  P
C     COMMON / LABELS     /  SCHEME, ANSATZ
C     COMMON / C          /  C
C     COMMON / F2         /  F2 
C  INPUTS
C          XI -- Bjorken x
C        DEL2 -- asymmetry parameter (P2-P1)^2. Set equal to 0 below!
C          Q2 -- photon virtuality squared
C         Q02 -- initial scale squared
C           P -- approximation order, which is N^{P}LO
C      SCHEME -- scheme
C      ANSATZ -- label for ansatz for GPDs on input scale
C           C -- intersection point of Mellin- Barnes integration 
C                path with real axis
C  OUTPUT
C     F2(0:2) --  F2(x,Q^2) in LO, NLO and NNLO approximations 
C                 (for SCHEME='CSBAA'; for SCHEME= 'MSBAR' only LO and
C                 NLO exists and CFF(3) is garbage) 
C                 NB!! - only one of three F2(3) array elements, F2(P)
C                 is changed by single call to F2F!
C  IDENTIFIERS
C     A,B,LAM -- Integration path is divided in two segments: [A, B] and [B, LAM]
C  CHILDREN
C     DQAGS    -- subroutine for adaptive numerical Gauss-Kronrod integration
C                 available in SLATEC library
C     F2IMAG    -- (via DQAGS) Mellin-Barnes integrand 
C     BETAF
C  SOURCE
C

      SUBROUTINE F2F

      IMPLICIT NONE
      INTEGER P
      DOUBLE PRECISION XI, DEL2, Q2, Q02, C
      DOUBLE PRECISION F2(0:2)
      CHARACTER SCHEME*5, ANSATZ*6
      DOUBLE PRECISION NG, NSEA, MG, MSEA, ALPHA0G, ALPHA0SEA
      DOUBLE PRECISION FITPAR(10)
*   Integration parameters:
      INTEGER IER, IWORK, LAST, LENW, LIMIT, NEVAL
      DOUBLE PRECISION A, B, LAM, ABSERR, EPSABS, EPSREL, WORK
      DOUBLE PRECISION F2IMAG, RES, AUX
      PARAMETER (LIMIT = 400, LENW = 4 * LIMIT)
      DIMENSION IWORK(LIMIT), WORK(LENW)
      EXTERNAL F2IMAG

*   Input common-blocks 

      COMMON / FITPARAMS  /  FITPAR

*   Output common-blocks

      COMMON / APPROX     /  P
      COMMON / KINEMATICS /  XI, DEL2, Q2, Q02
      COMMON / LABELS     /  SCHEME, ANSATZ
      COMMON / C          /  C
      COMMON / GPD        /  NG, NSEA, MG, MSEA, ALPHA0G, ALPHA0SEA
      COMMON / F2P        /  F2

      C = 0.5d0

*   Integration limits and errors

      A = 0.0d0
      B = 1.0d0
      LAM = 10.0d0
      EPSABS = 1.0d-9
      EPSREL = 1.0d-2

*     We do fitting at NLO in CSBAR scheme. It is here that this
*     is specified and written in common blocks.

      ANSATZ = 'FIT'
      P = 1
      SCHEME = 'CSBAR'

*     Scales and kinematics

      Q02 = FITPAR(7)
*      (Fixing forward kinematics:)
      DEL2 = 0.0d0

*     GPD parameters

      NG = FITPAR(1)
      NSEA = FITPAR(2)
      MG = FITPAR(3)
      MSEA = FITPAR(4)
      ALPHA0G = FITPAR(5)
      ALPHA0SEA = FITPAR(6)

*   Initialization of QCD beta function coefficients

      CALL BETAF

*   Integration

      CALL DQAGS(F2IMAG, A, B, EPSABS, EPSREL, AUX, ABSERR,
     &   NEVAL, IER, LIMIT, LENW, LAST, IWORK, WORK)
      CALL DQAGS(F2IMAG, B, LAM, EPSABS, EPSREL*50, RES, ABSERR,
     &   NEVAL, IER, LIMIT, LENW, LAST, IWORK, WORK)
      RES = RES + AUX

*     0.0707... = (2/9) / (Pi)

      F2(P) = (1.0d0/XI)**(C) * RES * 0.07073553026306459d0

      RETURN
      END
C     ***



C     ****f* cff.f/F2IMAG
C  NAME
C    F2IMAG
C  DESCRIPTION
C    Mellin-Barnes integrand for F2. "IMAG" in the 
C    name is because it is by optical theorem related
C    to imaginary part of forward Compton scattering
C  SYNOPSIS
C    DOUBLE PRECISION FUNCTION F2IMAG(Y)
C
C    DOUBLE PRECISION Y
C  PARENTS
C       F2F (via DQAGS)
C  CHILDREN
C       PARWAVF
C  SOURCE
C
      DOUBLE PRECISION FUNCTION F2IMAG(Y)
      IMPLICIT NONE
      DOUBLE PRECISION Y
      INTEGER P
      DOUBLE PRECISION XI, DEL2, Q2, Q02, C
      CHARACTER SCHEME*5, ANSATZ*6
      DOUBLE COMPLEX J, FPW

      COMMON / KINEMATICS /  XI, DEL2, Q2, Q02
      COMMON / APPROX     /  P
      COMMON / LABELS     /  SCHEME, ANSATZ
      COMMON / C          /  C

      J = COMPLEX(C, Y)
      CALL PARWAVF (J, XI, DEL2, Q2, Q02, P, SCHEME, ANSATZ, FPW)
      F2IMAG = REALPART( (1.0d0/XI)**COMPLEX(0.0d0,Y) * FPW )

      RETURN
      END
C     ***

