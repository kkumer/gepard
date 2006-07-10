C     ****h* gepard/cff.f
C  FILE DESCRIPTION
C    calculation of Compton form factors
C
C    $Id$
C  BUGS
C     Some complex arithmetic (REALPART, COMPLEX, ...) could be
C     specific to GNU G77 compiler.
C     *******


C     ****s* cff.f/CFFF
C  NAME
C     CFFF  --   Compton form factors
C  DESCRIPTION
C    calculates singlet Compton form factor H by numerical Mellin-Barnes 
C    integration of conformal partial waves provided by subroutine PARWAV. 
C    |latex  \begin{eqnarray*}
C    |latex  {{\cal H}^{\rm S}}(\xi,\Delta^2,{\cal Q}^2)
C    |latex  &\!\!\!=\!\!\!& \frac{1}{2i}\int_{c-i \infty}^{c+ i \infty}\!
C    |latex  dj\,\xi^{-j-1} \left[i + \tan \left(\frac{\pi j}{2}\right) \right]
C    |latex  \mbox{\boldmath $C$}_{j}({\cal Q}^2/\mu^2,\alpha_s(\mu)) 
C    |latex  \mbox{\boldmath $H$}_{j} (\xi,\Delta^2,\mu^2) \\[2ex]
C    |latex  &\!\!\!=\!\!\!& \frac{(\xi/2)^{-c-1}}{\Gamma(3/2)}\Big\{i \Im\! m\: e^{i\phi} 
C    |latex  \int_{0}^{\infty} dy  (\xi/2)^{- y e^{i\phi}} 
C    |latex  \frac{\Gamma(5/2 + c + y e^{i\phi})}{\Gamma(3 + c + y e^{i\phi})}
C    |latex  \texttt{FPW}(j) \\[2ex]
C    |latex  & & + \Im\! m\:  e^{i\phi} 
C    |latex  \int_{0}^{\infty} dy  (\xi/2)^{- y e^{i\phi}}
C    |latex  \frac{\Gamma(5/2 + c + y e^{i\phi})}{\Gamma(3 + c + y e^{i\phi})}
C    |latex  \tan\big(\frac{\pi(c + y e^{i\phi})}{2}\big)
C    |latex  \texttt{FPW}(j)
C    |latex  \end{eqnarray*} 
C    |latex  where $j=c+ye^{i\phi}$.
C    |latex  This is implemented now with fixed $\phi = \pi/2$, but
C    |latex  implementation with, say, $\phi = 3\pi /4$ should give faster
C    |latex  convergence.
C    All INPUTS and OUTPUT below is via common blocks!
C  SYNOPSIS
C     SUBROUTINE CFFF
C
C     INTEGER P
C     DOUBLE PRECISION XI, DEL2, Q2, Q02, C
C     DOUBLE COMPLEX CFF(0:2)
C     CHARACTER SCHEME*5, ANSATZ*6
C
C     COMMON / KINEMATICS /  XI, DEL2, Q2, Q02
C     COMMON / APPROX     /  P
C     COMMON / LABELS     /  SCHEME, ANSATZ
C     COMMON / C          /  C
C     COMMON / CFF        /  CFF
C  INPUTS
C          XI -- DVCS scaling parameter
C        DEL2 -- DVCS asymmetry parameter (P2-P1)^2
C          Q2 -- photon virtuality squared
C         Q02 -- initial scale squared
C           P -- approximation order, which is N^{P}LO
C      SCHEME -- scheme
C      ANSATZ -- label for ansatz for GPDs on input scale
C           C -- intersection point of Mellin- Barnes integration 
C                path with real axis
C  OUTPUT
C    CFF(0:2) --  Compton form factor in LO, NLO and NNLO approximations 
C                 (for ANSATZ='CSBAA'; for ANSATZ= 'MSBAR' only LO and
C                 NLO exists and CFF(3) is garbage) 
C                 NB!! - only one of three CFF(3) array elements, CFF(P)
C                 is changed by single call to CFFF!
C  IDENTIFIERS
C     A,B,LAM -- Integration path is divided in two segments: [A, B] and [B, LAM]
C  CHILDREN
C     DQAGS    -- subroutine for adaptive numerical Gauss-Kronrod integration
C                 available in SLATEC library
C     FIMAG    -- (via DQAGS) integrand for imaginary part
C     FREAL    -- (via DQAGS) integrand for real part
C     BETAF
C  SOURCE
C

      SUBROUTINE CFFF

      IMPLICIT NONE
      INTEGER P, PMAX, NPTSMAX, SPEED, K
      DOUBLE PRECISION XI, DEL2, Q2, Q02, C, PHI
      DOUBLE COMPLEX CFF(0:2), EPH, PIHALF, J, FPW
      DOUBLE COMPLEX DCTAN
      CHARACTER SCHEME*5, ANSATZ*6
      DOUBLE PRECISION RESREAL, RESIMAG
      DOUBLE PRECISION FREAL, FIMAG
      PARAMETER (NPTSMAX = 64)
      DOUBLE PRECISION Y(NPTSMAX), WG(NPTSMAX)
      DOUBLE COMPLEX N(NPTSMAX)

*   Input common-blocks 

      COMMON / INITPAR    /  SPEED, PMAX
      COMMON / APPROX     /  P
      COMMON / KINEMATICS /  XI, DEL2, Q2, Q02
      COMMON / LABELS     /  SCHEME, ANSATZ
      COMMON / CPHI       /  C, PHI
      COMMON / POINTS     /  Y, WG
      COMMON / NPOINTS    /  N

*   Output common-blocks

      COMMON / CFF        /  CFF

      PIHALF = (1.5707963267948966d0, 0.0d0)
      EPH = EXP ( COMPLEX(0.0d0, PHI) )

      RESREAL = 0.0d0
      RESIMAG = 0.0d0

*   Integration ...

      DO 123 K = 1, NPTSMAX/SPEED

      J = N(K) - 1

      CALL PARWAVF (K, XI, DEL2, Q2, Q02, P, SCHEME, ANSATZ, FPW)

      FREAL = IMAGPART(EPH * (2.0d0/XI)**(J-C) * DCTAN(PIHALF*J) * FPW)

      FIMAG = IMAGPART(EPH * (2.0d0/XI)**(J-C) * FPW)

      RESREAL = RESREAL + WG(K)*FREAL
      RESIMAG = RESIMAG + WG(K)*FIMAG

 123  CONTINUE

      CFF(P) = (2.0d0/XI)**(C + 1.0d0) * 
     &      COMPLEX(RESREAL, RESIMAG) / 0.886226925452758014d0

      RETURN
      END
C     ***
