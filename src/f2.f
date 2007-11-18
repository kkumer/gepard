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
C    integration 
C    Unlike CFF, this is defined with the singlet charge factor Q_{S}^2
C    included
C         NB!! - only one of three F2(3) array elements, F2(P)
C         is changed by single call to F2F!
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
C  SYNOPSIS

      SUBROUTINE F2F

C  PARENTS
C     PROCDATA
C  CHILDREN
C     LOOKUPQ
C  SOURCE
C

      IMPLICIT NONE
      INTEGER K, QIND
      DOUBLE PRECISION RES, F2IMAG
      DOUBLE COMPLEX EPH, FPW, J, FCM(2)
      INCLUDE 'header.f'


      CALL LOOKUPQ(QIND)

*     Fixing forward kinematics:
      MTIND = 0

      EPH = EXP ( COMPLEX(0.0d0, PHI) )
      RES = 0.0d0

      DO 123 K = 1, NPTS
      J = N(K) - 1
      FCM(1) = HGRID(MTIND, K, 1)
      FCM(2) = HGRID(MTIND, K, 2)
      FPW = CGRIDDIS(QIND,K,1) * FCM(1) + CGRIDDIS(QIND,K,2) * FCM(2)
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
