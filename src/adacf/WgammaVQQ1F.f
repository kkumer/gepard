*    
       SUBROUTINE WgammaVQQ1F (NF, n, res)
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
               res=2.5d-1*((-1.6d1*CF*NF*TF*(8.d0+2.8000000000000003d1*
     &  n+4.9d1*n**3+3.2d1*n**4+5.d0*n**5+3.8d1*(n*n)))/((-1.d0+
     &  n)*n**3*(1.d0+n)**3*(2.d0+n)**2)+CF*NF*TF*(1.33333333333
     &  33333d0-1.777777777777778d1*S1+1.0666666666666667d1*S2+(
     &  1.7777777777777777d0*(-3.d0+5.d0*n+1.1d1*(n*n)))/(n**2*(
     &  1.d0+n)**2))+CF*CA*(-7.166666666666667d0+(5.955555555555
     &  556d1+(8.d0*(1.d0+2.d0*n))/(n**2*(1.d0+n)**2))*S1+(-1.73
     &  33333333333332d1+8.d0/(n*(1.d0+n)))*S2-1.6d1*S1*S2-(4.44
     &  4444444444445d-1*(9.d0+3.d0*n+2.63d2*n**3+1.51d2*n**4+9.
     &  7d1*(n*n)))/(n**3*(1.d0+n)**3))+(-5.d-1*CF*CA+CF*CF)*(-3
     &  .d0+(1.6d1*(1.d0+2.d0*n)*S1)/(n**2*(1.d0+n)**2)+2.4d1*S2
     &  -(8.d0*(-1.d0+3.d0*n**3+n*n))/(n**3*(1.d0+n)**3)-(1.6d1*
     &  (1.d0+2.d0*n+2.d0*(n*n)))/(n**3*(1.d0+n)**3)+1.6d1*(-(1/
     &  (n*(1.d0+n)))+2.d0*S1)*(S2-HS2(5.d-1*n))-8.d0*HS3(5.d-1*
     &  n)+6.4d1*(S1/n**2-6.25d-1*ZETA3+MellinF2(
     &  n)-5.d-1*ZETA2*(-PSI(5.d-1*n)+PSI(
     &  5.d-1*(1.d0+n))))))
*
       RETURN
       END
