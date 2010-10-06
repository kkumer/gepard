
      DOUBLE PRECISION FUNCTION ImHoptR(p, t, x)

      IMPLICIT NONE
      DOUBLE PRECISION p(0:19), t, x
      DOUBLE PRECISION twox, onex, val, sea, PI
Cf2py intent(in) p
Cf2py intent(in) t, x
      PARAMETER(PI=3.1415926535897932385d0)
      twox = 2.*x / (1.+x)
      onex = (1.-x) / (1.+x)
      val = ( p(6) * p(10) * twox**(-p(7)-p(8)*t) * 
     &         onex**p(11) / (1. - onex*t/(p(9)**2))  )
      sea = ( 2./9. * p(0) * p(4) * twox**(-p(1)-p(2)*t) * 
     &         onex**p(5) / (1. - onex*t/(p(3)**2))**2 )
      ImHoptR = PI * (val + sea) / (1.+x)
*     write (22,*)  x, val, sea, ImHoptR
*     write (22,*)  p(0), p(1), p(4)
      return 
      end


      DOUBLE PRECISION FUNCTION ImHtoptR(p, t, x)

      IMPLICIT NONE
      DOUBLE PRECISION p(0:19), t, x
      DOUBLE PRECISION twox, onex, val, PI
Cf2py intent(in) p
Cf2py intent(in) t, x
      PARAMETER(PI=3.1415926535897932385d0)
      twox = 2.*x / (1.+x)
      onex = (1.-x) / (1.+x)
      val = ( p(14) * p(16) * twox**(-p(7)-p(8)*t) * 
     &         onex**p(17) / (1. - onex*t/(p(15)**2))  )
      ImHtoptR = PI * val / (1.+x)
      return 
      end


      DOUBLE PRECISION FUNCTION subtractionR(p, t)

      IMPLICIT NONE
      DOUBLE PRECISION p(0:19), t
Cf2py intent(in) p
Cf2py intent(in) t
      DOUBLE PRECISION PI
      PARAMETER(PI=3.1415926535897932385d0)
      subtractionR = p(12)/(1.-t/p(13)**2)**2
      return 
      end

      SUBROUTINE PVquadrature(fname, p, t, xi, res)

      IMPLICIT NONE
      CHARACTER*2 fname
      DOUBLE PRECISION p(0:19)
      DOUBLE PRECISION t, xi, res
Cf2py intent(in) fname
Cf2py intent(in) p
Cf2py intent(in) t, xi
Cf2py intent(out) res
      INTEGER K, NMAX
      PARAMETER(NMAX=18)
      DOUBLE PRECISION a, b, ga, y, u, intg, PI
      DOUBLE PRECISION fu, fxi, ImHoptR, ImHtoptR
      DOUBLE PRECISION roots(NMAX), weights(NMAX)
      PARAMETER(PI=3.1415926535897932385d0)
      DATA (roots(K), K=1, NMAX) /
     & -0.99156517, -0.95582395, -0.89260247, -0.80370496, -0.69168704,
     & -0.55977083, -0.41175116, -0.25188623, -0.08477501,  0.08477501,
     &  0.25188623,  0.41175116,  0.55977083,  0.69168704,  0.80370496,
     &  0.89260247,  0.95582395,  0.99156517 /
      DATA (weights(K), K=1, NMAX) /
     &  0.02161601,  0.04971455,  0.07642573,  0.10094204,  0.12255521,
     &  0.14064291,  0.15468468,  0.16427648,  0.16914238,  0.16914238,
     &  0.16427648,  0.15468468,  0.14064291,  0.12255521,  0.10094204,
     &  0.07642573,  0.04971455,  0.02161601 /
      DATA a, b, ga /0.d0, 1.0d0, 0.9d0/

      res = 0.0d0
      IF (fname .EQ. 'H') THEN
        fxi = ImHoptR(p, t, xi)
      ELSE IF (fname .EQ. 'E') THEN
        fxi = 0
      ELSE IF (fname .EQ. 'Ht') THEN
        fxi = ImHtoptR(p, t, xi)
      ELSE  ! (fname .EQ. 'Et')
        fxi =  0
      ENDIF
      DO 20 K = 1, NMAX
        y = (b-a)*(roots(K)+1)/2.0 + a
        u = y**(1./(1.-ga))
        IF (fname .EQ. 'H') THEN
          fu = ImHoptR(p, t, u)
        ELSE IF (fname .EQ. 'E') THEN
          fu = 0
        ELSE IF (fname .EQ. 'Ht') THEN
          fu = ImHtoptR(p, t, u)
        ELSE  ! (fname .EQ. 'Et')
          fu =  0
        ENDIF
        intg = u**ga * ( fu - fxi )
*       write (21,*)  u, fu
        IF ((fname .EQ. 'Ht') .OR. (fname .EQ. 'Et')) THEN
          intg = (2.*xi) / (xi**2 - u**2) * intg / (1.-ga)
        ELSE   ! H or E
          intg = (2.*u) / (xi**2 - u**2) * intg / (1.-ga)
        ENDIF
        res = res + weights(K)*intg
20    CONTINUE
      res = (b-a)/2.0*res
      IF ((fname .EQ. 'Ht') .OR. (fname .EQ. 'Et')) THEN
        res = res + log((1.+xi)/(1.-xi)) * fxi
      ELSE
        res = res + log(xi**2 / (1.-xi**2)) * fxi
      ENDIF
      res = res/PI
      return
      end





*        0  -    NS 
*        1  -   alS 
*        2  -  alpS 
*        3  -    MS 
*        4  -    rS 
*        5  -    bS 
*        6  -    Nv 
*        7  -   alv 
*        8  -  alpv 
*        9  -    Mv 
*       10  -    rv 
*       11  -    bv 
*       12  -     C 
*       13  -    MC 
*       14  -   tNv 
*       15  -   tMv 
*       16  -   trv 
*       17  -   tbv

