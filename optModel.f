
      SUBROUTINE ImHopt(p, t, x, res)

      IMPLICIT NONE
      DOUBLE PRECISION p(0:17), t, x(18), res(18)
      DOUBLE PRECISION twox, onex, val, sea, PI
Cf2py intent(in) p, x
Cf2py intent(in) t
Cf2py intent(out) res
      INTEGER K, NMAX
      PARAMETER(NMAX=18)
      PARAMETER(PI=3.1415926535897932385d0)
      DO 10 K = 1, NMAX
        twox = 2.*x(K) / (1.+x(K))
        onex = (1.-x(K)) / (1.+x(K))
        val = ( p(6) * p(10) * twox**(-p(7)-p(8)*t) * 
     &         onex**p(11) / (1. - onex*t/(p(9)**2))  )
        sea = ( 2./9. * p(0) * p(4) * twox**(-p(1)-p(2)*t) * 
     &         onex**p(5) / (1. - onex*t/(p(3)**2))**2 )
        res(K) = PI * (val + sea) / (1.+x(K))
10       CONTINUE
      return
      end


      DOUBLE PRECISION FUNCTION ImHoptR(p, t, x)

      IMPLICIT NONE
      DOUBLE PRECISION p(0:17), t, x
      DOUBLE PRECISION twox, onex, val, sea, PI
Cf2py intent(in) p
Cf2py intent(in) t, x
Cf2py intent(out) res
      PARAMETER(PI=3.1415926535897932385d0)
      twox = 2.*x / (1.+x)
      onex = (1.-x) / (1.+x)
      val = ( p(6) * p(10) * twox**(-p(7)-p(8)*t) * 
     &         onex**p(11) / (1. - onex*t/(p(9)**2))  )
      sea = ( 2./9. * p(0) * p(4) * twox**(-p(1)-p(2)*t) * 
     &         onex**p(5) / (1. - onex*t/(p(3)**2))**2 )
      ImHoptR = PI * (val + sea) / (1.+x)
      return 
      end



      SUBROUTINE PVquadratureH(p, t, xi, res)

      IMPLICIT NONE
      DOUBLE PRECISION p(0:17)
      DOUBLE PRECISION t, xi, res
Cf2py intent(in) p
Cf2py intent(in) t, xi
Cf2py intent(out) res
      INTEGER K, NMAX
      PARAMETER(NMAX=18)
      DOUBLE PRECISION a, b, ga, y, u, intg, ImHoptR
      DOUBLE PRECISION roots(NMAX), weights(NMAX)
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
      DO 20 K = 1, NMAX
        y = (b-a)*(roots(K)+1)/2.0 + a

        u = y**(1./(1.-ga))
        intg = u**ga * ( ImHoptR(p, t, u) - ImHoptR(p, t, xi) )
        intg = (2.*u) / (xi**2 - u**2) * intg / (1.-ga)

        res = res + weights(K)*intg
20    CONTINUE
      res = (b-a)/2.0*res
      res = res + log(xi**2 / (1.-xi**2)) * ImHoptR(p, t, xi)
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

