C     ****h* gepard/cff.f
C  FILE DESCRIPTION
C    calculation of Compton form factors
C
C    $Id$
C     *******


C     ****s* cff.f/CFFF
C  NAME
C     CFFF  --   Compton form factor \mathcal{H}
C  DESCRIPTION
C    calculates singlet Compton form factor H by numerical Mellin-Barnes 
C    integration according to formula Eq. (58) from [Kumericki:2007sa]
C      NB!! - only one of three CFF(0:2) array elements, CFF(P)
C             is changed by single call to CFFF!
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
C  SYNOPSIS

      SUBROUTINE CFFF

C  PARENTS
C    PARSIGMA, cffHInternal, cffEInternal
C  CHILDREN
C    GETMBGPD, GETMBGPDMMA, LOOKUPQ, DCTAN
C  SOURCE
C


      IMPLICIT NONE
      INTEGER K, QIND
      DOUBLE COMPLEX EPH, J, CFAC, FPWH, FPWE
      DOUBLE COMPLEX DCTAN, FCM(2)
      DOUBLE PRECISION RESREALH, RESIMAGH
      DOUBLE PRECISION FREALH, FIMAGH
      DOUBLE PRECISION RESREALE, RESIMAGE
      DOUBLE PRECISION FREALE, FIMAGE
      INCLUDE 'header.f'

      EPH = EXP ( COMPLEX(0.0d0, PHI) )

      RESREALH = 0.0d0
      RESIMAGH = 0.0d0
      RESREALE = 0.0d0
      RESIMAGE = 0.0d0

      CALL LOOKUPQ(QIND)
#ifdef MMA
      IF (ANSATZ(:3) .EQ. 'MMA') THEN
        CALL GETMBGPDMMA
      ELSE
        CALL GETMBGPD
      END IF
#else
      IF (ANSATZ(:3) .EQ. 'MMA') THEN
        WRITE(*,*) 'ANSATZ=MMA unavailable. Recompile with -D MMA'
      ELSE
        CALL GETMBGPD
      END IF
#endif

*   Integration ...

      DO 123 K = 1, NPTS

      J = N(K) - 1
      CFAC = EPH / XI**(J+1.d0)

      IF ( FFTYPE(:7) .EQ. 'NONSING' ) THEN
        FPWH = CGRIDNS(QIND,K)*MBGPD(K,1)        
      ELSE
        FPWH = CGRID(QIND,K,1)*MBGPD(K,1) + CGRID(QIND,K,2)*MBGPD(K,2)
        FPWE = CGRID(QIND,K,1)*PAR(18)*MBGPD(K,1) 
     &       + CGRID(QIND,K,2)*PAR(28)*MBGPD(K,2)
      END IF

      FREALH = IMAGPART(CFAC * DCTAN(PIHALF*J) * FPWH)
      FREALE = IMAGPART(CFAC * DCTAN(PIHALF*J) * FPWE)

      FIMAGH = IMAGPART(CFAC * FPWH)
      FIMAGE = IMAGPART(CFAC * FPWE)

      RESREALH = RESREALH + WG(K)*FREALH
      RESIMAGH = RESIMAGH + WG(K)*FIMAGH
      
      RESREALE = RESREALE + WG(K)*FREALE
      RESIMAGE = RESIMAGE + WG(K)*FIMAGE

 123  CONTINUE
      
*     Multiplying with charge factor

      CFF(P) = CHARGEFAC * COMPLEX(RESREALH, RESIMAGH)
      CFFE(P) = CHARGEFAC * COMPLEX(RESREALE, RESIMAGE) 

      RETURN
      END
C     ***

