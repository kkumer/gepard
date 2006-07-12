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
      INTEGER SPEED, P, NF
      CHARACTER SCHEME*5, ANSATZ*6
      INTEGER NPTSMAX, K
      DOUBLE PRECISION XI, DEL2, Q2, Q02, C, PHI
      DOUBLE PRECISION F2(0:2)
      DOUBLE PRECISION NG, NSEA, MG, MSEA, ALPHA0G, ALPHA0SEA
      DOUBLE PRECISION FITPAR(10), RES, F2IMAG
      DOUBLE COMPLEX EPH, FPW, J
      PARAMETER (NPTSMAX = 64)
      DOUBLE PRECISION Y(NPTSMAX), WG(NPTSMAX)
      DOUBLE COMPLEX N(NPTSMAX)

*   Input common-blocks 

      COMMON / PARINT /  SPEED, P, NF
      COMMON / PARCHR /  SCHEME, ANSATZ

      COMMON / FITPARAMS  /  FITPAR
      COMMON / CPHI       /  C, PHI
      COMMON / POINTS     /  Y, WG
      COMMON / NPOINTS    /  N

*   Output common-blocks

      COMMON / KINEMATICS /  XI, DEL2, Q2, Q02
      COMMON / GPD        /  NG, NSEA, MG, MSEA, ALPHA0G, ALPHA0SEA
      COMMON / F2P        /  F2


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

      EPH = EXP ( COMPLEX(0.0d0, PHI) )
      RES = 0.0d0

      DO 123 K = 1, NPTSMAX/SPEED
      J = N(K) - 1
      CALL PARWAVF (K, FPW, 'DIS')
      F2IMAG = IMAGPART(EPH * (1.0d0/XI)**(J-C) * FPW )
      RES = RES + WG(K)*F2IMAG
 123  CONTINUE

*     0.0707... = (2/9) / (Pi)

      F2(P) = (1.0d0/XI)**(C) * RES * 0.07073553026306459d0

      RETURN
      END
C     ***
