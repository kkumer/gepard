*    
       SUBROUTINE WgammaVGG2F (NF, n, res)
*
       IMPLICIT none
       DOUBLE COMPLEX n, res
       DOUBLE COMPLEX S1, S2, S3, S4
       DOUBLE COMPLEX HS1, HS2, HS3, HS4
       DOUBLE COMPLEX PSI, DPSI
       DOUBLE PRECISION NF, CF, CA, TF
       DOUBLE PRECISION EMC, ZETA2, ZETA3, ZETA4, LOG2
       PARAMETER ( EMC  = 0.5772156649D0, ZETA2 = 1.6449340668D0,
     ,            ZETA3 = 1.2020569032D0, ZETA4 = 1.0823232337D0,
     ,            LOG2 = 0.6931471806D0 )
       PARAMETER ( CF=4./3.d0, CA=3.d0, TF=0.5d0 )
*
*      ... Harmonic sum initialization ...
*
      COMMON /HARMONIC/ S1, S2, S3, S4
*
*     ... Function itself spliced in by Mathematica ...
               res=2.5d-1*(-4.425894d3+2.6758000000000006d3/(-1.d0+n)**
     &  2-1.4213999999999998d4/(-1.d0+n)+3.4560000000000004d3/n*
     &  *5+4.32d2/n**4+1.4942d4/n**3+2.7439999999999998d2/n**2+2
     &  .0852d4/n-3.968d3/(1.d0+n)+3.363d3/(2.d0+n)-4.848d3/(3.d
     &  0+n)+(3.589d3*S1)/n+2.6435210000000002d3*(-(1/n)+S1)-7.3
     &  05000000000001d3*(S1/n**2+(S2-ZETA2)/n)-(1.7513999999999
     &  998d4*(-(S1/n**2)-S2/n-S3+ZETA2/n+ZETA3))/n-NF*(-5.28722
     &  9999999999d2-1.5727000000000002d2/(-1.d0+n)**2+1.8296000
     &  000000001d2/(-1.d0+n)+4.551111111111111d2/n**5-5.5466666
     &  66666667d2/n**4+9.826d2/n**3-1.541d3/n**2-3.502000000000
     &  0002d2/n+7.557d2/(1.d0+n)-7.138d2/(2.d0+n)+5.59299999999
     &  9999d2/(3.d0+n)+(3.2d2*S1)/n+4.121720000000001d2*(-(1/n)
     &  +S1)+2.615d1*(S1/n**2+(S2-ZETA2)/n)-(1.6174000000000002d
     &  3*(-(S1/n**2)-S2/n-S3+ZETA2/n+ZETA3))/n)-(6.463d0-2.7983
     &  539094650207d0/(-1.d0+n)+7.111111111111111d0/n**4+1.936d
     &  1/n**3+3.422d0/n**2-1.3878000000000001d1/n+1.534d2/(1.d0
     &  +n)-1.8769999999999998d2/(2.d0+n)+5.2749999999999995d1/(
     &  3.d0+n)+1.7777777777777777d0*(-(1/n)+S1)-1.156d2*(S1/n**
     &  2+(S2-ZETA2)/n)+8.525d1*((1/(1.d0+n)+S1)/(1.d0+n)**2+((1
     &  .d0+n)**(-2)+S2-ZETA2)/(1.d0+n))-(1.2646d2*(-(S1/n**2)-S
     &  2/n-S3+ZETA2/n+ZETA3))/n)*(NF*NF))
*
       RETURN
       END
