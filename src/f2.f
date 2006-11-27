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
C    |latex  {{F_2}^{\rm S}}(\xi,{\cal Q}^2)
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
C  OUTPUT
C     F2(0:2) --  F2(x,Q^2) in LO, NLO and NNLO approximations 
C                 (for SCHEME='CSBAA'; for SCHEME= 'MSBAR' only LO and
C                 NLO exists and CFF(3) is garbage) 
C                 NB!! - only one of three F2(3) array elements, F2(P)
C                 is changed by single call to F2F!
C  IDENTIFIERS
C          XI -- Bjorken x
C        DEL2 -- asymmetry parameter (P2-P1)^2. Set equal to 0 below!
C          Q2 -- photon virtuality squared
C         Q02 -- initial scale squared
C           P -- approximation order, which is N^{P}LO
C      SCHEME -- scheme
C      ANSATZ -- label for ansatz for GPDs on input scale
C           C -- intersection point of Mellin- Barnes integration 
C                path with real axis
C         PHI -- angle between Mellin-Barnes contour and Re(J) axis
C  PARENTS
C     PROCDATA, TEST
C  CHILDREN
C     PARWAVF
C  SOURCE
C

      SUBROUTINE F2F

      IMPLICIT NONE
      INTEGER K
      DOUBLE PRECISION RES, F2IMAG
      DOUBLE COMPLEX EPH, FPW, J
      INCLUDE 'header.f'



*     Fixing forward kinematics:
      DEL2 = 0.0d0

      EPH = EXP ( COMPLEX(0.0d0, PHI) )
      RES = 0.0d0

      DO 123 K = 1, NPTS
      J = N(K) - 1
      CALL PARWAVF (K, FPW, 'DIS')
      F2IMAG = IMAGPART(EPH * (1.0d0/XI)**(J-C) * FPW )
      RES = RES + WG(K)*F2IMAG
 123  CONTINUE

      IF (NF .EQ. 3) THEN
         F2(P) = (1.0d0/XI)**(C) * RES * 2.0D0 / 9.0D0 / PI
      ELSE IF (NF .EQ. 4) THEN
         F2(P) = (1.0d0/XI)**(C) * RES * 5.0D0 / 18.0D0 / PI
      ELSE
         CALL ERROR ('GeParD', '  CFF',
     &   'NF is not integer equal to 3 or 4!                          ',
     &   5, 3)
      END IF

      RETURN
      END
C     ***
