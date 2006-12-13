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
C  OUTPUT
C    CFF(0:2) --  Compton form factor in LO, NLO and NNLO approximations 
C                 (for SCHEME='CSBAA'; for SCHEME= 'MSBAR' only LO and
C                 NLO exists and CFF(3) is garbage) 
C                 NB!! - only one of three CFF(3) array elements, CFF(P)
C                 is changed by single call to CFFF!
C  IDENTIFIERS
C          XI -- DVCS scaling parameter
C           K -- Mellin-Barnes integration point index
C           J -- conformal moment
C          XI -- DVCS scaling parameter
C           C -- intersection point of Mellin- Barnes integration 
C                path with real axis
C         PHI -- angle between Mellin-Barnes contour and Re(J) axis
C  CHILDREN
C    PARWAV, DCTAN
C  SOURCE
C

      SUBROUTINE CFFF

      IMPLICIT NONE
      INTEGER K
      DOUBLE COMPLEX EPH, PIHALF, J, FPW
      DOUBLE COMPLEX DCTAN
      DOUBLE PRECISION RESREAL, RESIMAG
      DOUBLE PRECISION FREAL, FIMAG
      INCLUDE 'header.f'

      EPH = EXP ( COMPLEX(0.0d0, PHI) )

      RESREAL = 0.0d0
      RESIMAG = 0.0d0

*   Integration ...

      DO 123 K = 1, NPTS

      J = N(K) - 1

      CALL PARWAVF(K, FPW)

      PIHALF = (1.5707963267948966d0, 0.0d0)
      FREAL = IMAGPART(EPH * (2.0d0/XI)**(J-C) * DCTAN(PIHALF*J) * FPW)

      FIMAG = IMAGPART(EPH * (2.0d0/XI)**(J-C) * FPW)

      RESREAL = RESREAL + WG(K)*FREAL
      RESIMAG = RESIMAG + WG(K)*FIMAG

 123  CONTINUE
      
*   Gamma(3/2) = 0.8862....

      CFF(P) = (2.0d0/XI)**(C + 1.0d0) * 
     &      COMPLEX(RESREAL, RESIMAG) / 0.886226925452758014d0

      RETURN
      END
C     ***
