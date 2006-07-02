C     ****h* gepard/cff.f
C  FILE DESCRIPTION
C    calculation of Compton form factors
C
C    $Id: cff.f,v 1.1 2006-04-12 15:46:02+02 kkumer Exp kkumer $
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
      INTEGER P
      DOUBLE PRECISION XI, DEL2, Q2, Q02, C
      DOUBLE COMPLEX CFF(0:2)
      CHARACTER SCHEME*5, ANSATZ*6
*   Integration parameters:
      INTEGER IER, IWORK, LAST, LENW, LIMIT, NEVAL
      DOUBLE PRECISION A, B, LAM, ABSERR, EPSABS, EPSREL, WORK
      DOUBLE PRECISION FREAL, RESREAL, FIMAG, RESIMAG, AUX
      PARAMETER (LIMIT = 400, LENW = 4 * LIMIT)
      DIMENSION IWORK(LIMIT), WORK(LENW)
      EXTERNAL FREAL, FIMAG

*   Input common-blocks 

C      SCHEME -- scheme
      COMMON / APPROX     /  P
      COMMON / KINEMATICS /  XI, DEL2, Q2, Q02
      COMMON / LABELS     /  SCHEME, ANSATZ

*   Output common-blocks

      COMMON / C          /  C
      COMMON / CFF        /  CFF

      C = 0.5d0

*   Integration limits and errors

      A = 0.0d0
      B = 1.0d0
      LAM = 10.0d0
      EPSABS = 1.0d-9
      EPSREL = 1.0d-2

*   Initialization of QCD beta function coefficients

      CALL BETAF

*   Integration

      CALL DQAGS(FREAL, A, B, EPSABS, EPSREL, AUX, ABSERR,
     &   NEVAL, IER, LIMIT,  LENW, LAST, IWORK, WORK)
      CALL DQAGS(FREAL, B, LAM, EPSABS, EPSREL*50, RESREAL, ABSERR,
     &   NEVAL, IER, LIMIT, LENW, LAST, IWORK, WORK)
      RESREAL = RESREAL + AUX

      CALL DQAGS(FIMAG, A, B, EPSABS, EPSREL, AUX, ABSERR,
     &   NEVAL, IER, LIMIT, LENW, LAST, IWORK, WORK)
      CALL DQAGS(FIMAG, B, LAM, EPSABS, EPSREL*50, RESIMAG, ABSERR,
     &   NEVAL, IER, LIMIT, LENW, LAST, IWORK, WORK)
      RESIMAG = RESIMAG + AUX

      CFF(P) = (2.0d0/XI)**(C + 1.0d0) * 
     &      COMPLEX(RESREAL, RESIMAG) / 0.886226925452758014d0

      RETURN
      END
C     ***


C     ****f* cff.f/FREAL
C  NAME
C    FREAL
C  DESCRIPTION
C    real part of the Mellin-Barnes integrand
C  SYNOPSIS
C    DOUBLE PRECISION FUNCTION FREAL(Y)
C
C    DOUBLE PRECISION Y
C  PARENTS
C       CFFF (via DQAGS)
C  CHILDREN
C       PARWAVF, CLNGAMMA, DCTAN
C  SOURCE
C

      DOUBLE PRECISION FUNCTION FREAL(Y)
      IMPLICIT NONE
      DOUBLE PRECISION Y
      INTEGER P
      DOUBLE PRECISION XI, DEL2, Q2, Q02, C
      CHARACTER SCHEME*5, ANSATZ*6
      DOUBLE COMPLEX J, FPW, PIHALF, CLNGAMMA, DCTAN

      COMMON / KINEMATICS /  XI, DEL2, Q2, Q02
      COMMON / APPROX     /  P
      COMMON / LABELS     /  SCHEME, ANSATZ
      COMMON / C          /  C

      PIHALF = (1.5707963267948966d0, 0.0d0)
      J = COMPLEX(C, Y)
      CALL PARWAVF (J, XI, DEL2, Q2, Q02, P, SCHEME, ANSATZ, FPW)
      FREAL = REALPART( (2.0d0/XI)**COMPLEX(0.0d0,Y) * DCTAN(PIHALF*J) *
     &  FPW * EXP(CLNGAMMA(2.5d0 + J)) / EXP(CLNGAMMA(3.0d0 + J)) )

      RETURN
      END
C     ***


C     ****f* cff.f/FIMAG
C  NAME
C    FIMAG
C  DESCRIPTION
C    imaginary part of the Mellin-Barnes integrand
C  SYNOPSIS
C    DOUBLE PRECISION FUNCTION FIMAG(Y)
C
C    DOUBLE PRECISION Y
C  PARENTS
C       CFFF (via DQAGS)
C  CHILDREN
C       PARWAVF, CLNGAMMA
C  SOURCE
C
      DOUBLE PRECISION FUNCTION FIMAG(Y)
      IMPLICIT NONE
      DOUBLE PRECISION Y
      INTEGER P
      DOUBLE PRECISION XI, DEL2, Q2, Q02, C
      CHARACTER SCHEME*5, ANSATZ*6
      DOUBLE COMPLEX J, FPW, CLNGAMMA

      COMMON / KINEMATICS /  XI, DEL2, Q2, Q02
      COMMON / APPROX     /  P
      COMMON / LABELS     /  SCHEME, ANSATZ
      COMMON / C          /  C

      J = COMPLEX(C, Y)
      CALL PARWAVF (J, XI, DEL2, Q2, Q02, P, SCHEME, ANSATZ, FPW)
      FIMAG = REALPART( (2.0d0/XI)**COMPLEX(0.0d0,Y) * FPW *
     &      EXP(CLNGAMMA(2.5d0 + J)) / EXP(CLNGAMMA(3.0d0 + J)) )

      RETURN
      END
C     ***

