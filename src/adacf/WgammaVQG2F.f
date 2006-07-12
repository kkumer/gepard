*    
       SUBROUTINE WgammaVQG2F (NF, n, res)
*
       IMPLICIT none
       DOUBLE COMPLEX n, res
       DOUBLE COMPLEX S1, S2, S3, S4
       DOUBLE COMPLEX HS1, HS2, HS3, HS4
       DOUBLE COMPLEX PSI, DPSI
       DOUBLE PRECISION NF, CF, CG, TF
       DOUBLE PRECISION EMC, ZETA2, ZETA3, ZETA4, LOG2
       PARAMETER ( EMC  = 0.57721 56649D0, ZETA2 = 1.64493 40668D0,
     ,            ZETA3 = 1.20205 69032D0, ZETA4 = 1.08232 32337D0,
     ,            LOG2 = 0.69314 71806D0 )
       PARAMETER ( CF=4./3.d0, CG=3.d0, TF=0.5d0 )
*
*      ... Harmonic sum initialization ...
*
      COMMON /HARMONIC/ S1, S2, S3, S4
*
*     ... Function itself spliced in by Mathematica ...
               res=2.5d-1*(-(NF*NF*(4.576131687242798d0/(-1.d0+n)-4.266
     &  666666666667d1/n**5+8.355555555555556d1/n**4-1.815999999
     &  9999998d2/n**3+2.54d2/n**2-2.52d2/n-7.019999999999999d1/
     &  (1.d0+n)**4-1.9613999999999998d2/(1.d0+n)**3+1.58d2/(1.d
     &  0+n)+1.454d2/(2.d0+n)-1.3928000000000003d2/(3.d0+n)+(5.4
     &  96d0*S1)/n-(7.4074074074074066d-1*(S1**3+3.d0*S1*S2+2.d0
     &  *S3))/n-5.309d1*(S1/n**2+(S2-ZETA2)/n)-(1.61232d2*(-(S1/
     &  n**2)-S2/n-S3+ZETA2/n+ZETA3))/n+(7.4074074074074066d0*(S
     &  2+S1*S1))/n))-NF*(2.986666666666667d2/(-1.d0+n)**2-1.268
     &  3d3/(-1.d0+n)+4.764444444444445d2/n**5+8.8d1/n**4+1.7630
     &  000000000001d3/n**3-4.249d2/n**2+2.522d3/n+1.515d3/(1.d0
     &  +n)**4-3.316d3/(1.d0+n)+2.1260000000000003d3/(2.d0+n)-(1
     &  .0442d2*S1)/n+(7.777777777777778d0*(S1**3+3.d0*S1*S2+2.d
     &  0*S3))/n+1.823d3*(S1/n**2+(S2-ZETA2)/n)-(5.044d1*(-(S1/n
     &  **2)-S2/n-S3+ZETA2/n+ZETA3))/n-(1.205d2*(S2+S1*S1))/n+3.
     &  7037037037037033d0*(S1**4/n+(8.d0*S1*S3)/n+(6.d0*S4)/n+(
     &  6.d0*S2*(S1*S1))/n+(3.d0*(S2*S2))/n)))
*
       RETURN
       END
