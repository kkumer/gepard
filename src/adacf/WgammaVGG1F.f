*    
       SUBROUTINE WgammaVGG1F (NF, n, res)
*
       IMPLICIT none
       DOUBLE COMPLEX n, res
       DOUBLE COMPLEX S1, S2, S3, S4
       DOUBLE COMPLEX HS1, HS2, HS3, HS4
       DOUBLE COMPLEX PSI, DPSI
       DOUBLE PRECISION NF, CF, CA, TF
       DOUBLE PRECISION EMC, ZETA2, ZETA3, ZETA4, LOG2
       DOUBLE COMPLEX MELLINF2
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
               res=2.5d-1*(CF*NF*TF*(8.d0+(1.6d1*(-4.d0-4.d0*n-1.d1*n**
     &  3+n**4+4.d0*n**5+2.d0*n**6-5.d0*(n*n)))/((-1.d0+n)*n**3*
     &  (1.d0+n)**3*(2.d0+n)))+CA*NF*TF*(1.0666666666666667d1-1.
     &  777777777777778d1*S1+(1.7777777777777777d0*(1.2d1+5.6000
     &  000000000005d1*n+7.6d1*n**3+3.8d1*n**4+9.399999999999999
     &  d1*(n*n)))/((-1.d0+n)*n**2*(1.d0+n)**2*(2.d0+n)))+CA*CA*
     &  (-2.1333333333333333d1+5.955555555555556d1*S1+(6.4d1*S1*
     &  (-2.d0-2.d0*n+8.d0*n**3+5.d0*n**4+2.d0*n**5+7.d0*(n*n)))
     &  /((-1.d0+n)**2*n**2*(1.d0+n)**2*(2.d0+n)**2)-(4.44444444
     &  4444445d-1*(5.76d2+1.488d3*n-1.6320000000000001d3*n**3-2
     &  .344d3*n**4+1.5670000000000002d3*n**5+6.098d3*n**6+6.04d
     &  3*n**7+2.742d3*n**8+4.57d2*n**9+5.6000000000000005d2*(n*
     &  n)))/((-1.d0+n)**2*n**3*(1.d0+n)**3*(2.d0+n)**3)-1.6d1*S
     &  1*HS2(5.d-1*n)+(3.2d1*(1.d0+n+n*n)*HS2(5.d-1*n))/((-1.d0
     &  +n)*n*(1.d0+n)*(2.d0+n))-4.d0*HS3(5.d-1*n)+3.2d1*(S1/n**
     &  2-6.25d-1*ZETA3+MELLINF2(
     &  n)-5.d-1*ZETA2*(-PSI(5.d-1*n)+PSI(5.d-1*(1.d0+n)))
     &  )))
*
       RETURN
       END
