*    
       SUBROUTINE WgammaVGQ1F (NF, n, res)
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
               res=2.5d-1*(-1.0666666666666667d1*CF*NF*TF*((1.d0+n)**(-
     &  2)+((-2.6666666666666665d0+S1)*(2.d0+n+n*n))/((-1.d0+n)*
     &  n*(1.d0+n)))-4.d0*(CF*CF)*((-4.d0*S1)/(1.d0+n)**2-(-4.d0
     &  -1.2d1*n+2.8000000000000003d1*n**3+4.3d1*n**4+3.d1*n**5+
     &  1.2d1*n**6-n*n)/((-1.d0+n)*n**3*(1.d0+n)**3)+((2.d0+n+n*
     &  n)*(1.d1*S1-2.d0*S2-2.d0*(S1*S1)))/((-1.d0+n)*n*(1.d0+n)
     &  ))-8.d0*CF*CA*((1.1111111111111112d-1*(1.44d2+4.32d2*n-1
     &  .3039999999999998d3*n**3-1.031d3*n**4+6.949999999999999d
     &  2*n**5+1.678d3*n**6+1.4000000000000001d3*n**7+6.21d2*n**
     &  8+1.09d2*n**9-1.52d2*(n*n)))/((-1.d0+n)**2*n**3*(1.d0+n)
     &  **3*(2.d0+n)**2)-(3.333333333333333d-1*S1*(-1.2d1-2.2d1*
     &  n+1.7000000000000002d1*n**4+4.1d1*(n*n)))/((-1.d0+n)**2*
     &  n**2*(1.d0+n))+((2.d0+n+n*n)*(S2+S1*S1-HS2(5.d-1*n)))/((
     &  -1.d0+n)*n*(1.d0+n))))
*
       RETURN
       END
