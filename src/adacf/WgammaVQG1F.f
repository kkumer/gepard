*    
       SUBROUTINE WgammaVQG1F (NF, n, res)
*
       IMPLICIT none
       DOUBLE COMPLEX n, res
       DOUBLE COMPLEX S1, S2, S3, S4
       DOUBLE COMPLEX HS1, HS2, HS3, HS4
       DOUBLE COMPLEX PSI, DPSI
       DOUBLE PRECISION NF, CF, CA, TF
       DOUBLE PRECISION EMC, ZETA2, ZETA3, ZETA4, LOG2
       PARAMETER ( EMC  = 0.57721 56649D0, ZETA2 = 1.64493 40668D0,
     ,            ZETA3 = 1.20205 69032D0, ZETA4 = 1.08232 32337D0,
     ,            LOG2 = 0.69314 71806D0 )
       PARAMETER ( CF=4./3.d0, CA=3.d0, TF=0.5d0 )
*
*      ... Harmonic sum initialization ...
*
      COMMON /HARMONIC/ S1, S2, S3, S4
*
*     ... Function itself spliced in by Mathematica ...
               res=2.5d-1*(-8.d0*CF*NF*TF*((-4.d0*S1)/n**2+(4.d0+8.d0*n
     &  +2.6d1*n**3+1.1d1*n**4+1.5d1*(n*n))/(n**3*(1.d0+n)**3*(2
     &  .d0+n))+((2.d0+n+n*n)*(5.d0-2.d0*S2+2.d0*(S1*S1)))/(n*(1
     &  .d0+n)*(2.d0+n)))-8.d0*CA*NF*TF*((8.d0*(3.d0+2.d0*n)*S1)
     &  /((1.d0+n)**2*(2.d0+n)**2)+(2.d0*(1.6d1+6.4d1*n+1.28d2*n
     &  **3+8.5d1*n**4+3.5999999999999996d1*n**5+2.5d1*n**6+1.5d
     &  1*n**7+6.d0*n**8+n**9+1.04d2*(n*n)))/((-1.d0+n)*n**3*(1.
     &  d0+n)**3*(2.d0+n)**3)+((2.d0+n+n*n)*(2.d0*S2-2.d0*(S1*S1
     &  )-2.d0*HS2(5.d-1*n)))/(n*(1.d0+n)*(2.d0+n))))
*
       RETURN
       END
