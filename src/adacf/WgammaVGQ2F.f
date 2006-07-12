*    
       SUBROUTINE WgammaVGQ2F (NF, n, res)
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
               res=2.5d-1*(1.1893d3/(-1.d0+n)**2-6.1631d3/(-1.d0+n)+1.2
     &  705185185185186d3/n**5+1.0453333333333332d3/n**4+3.588d3
     &  /n**3+4.0329999999999995d3/n**2+4.307d3/n+1.9458d3/(1.d0
     &  +n)**3-4.893d2/(1.d0+n)-1.452d3/(2.d0+n)-1.46d2/(3.d0+n)
     &  +(2.193d3*S1)/n+(8.148148148148147d1*(S1**3+3.d0*S1*S2+2
     &  .d0*S3))/n+(8.946000000000002d2*(-(S1/n**2)-S2/n-S3+ZETA
     &  2/n+ZETA3))/n-(6.063d2*(S2+S1*S1))/n-NF*(-7.108199999999
     &  999d1/(-1.d0+n)**2-4.641d1/(-1.d0+n)+1.1377777777777778d
     &  2/n**5-5.214814814814814d1/n**4+4.078d1/n**3-1.748000000
     &  0000002d2/n**2-1.838d2/n+2.1719999999999997d2/(1.d0+n)**
     &  3+3.335d1/(1.d0+n)-2.779d2/(2.d0+n)+(2.9669999999999996d
     &  2*S1)/n+(4.938271604938271d0*(S1**3+3.d0*S1*S2+2.d0*S3))
     &  /n-4.968d1*(S1/n**2+(S2-ZETA2)/n)-(6.806900000000001d1*(
     &  S2+S1*S1))/n)-NF*NF*(-2.3703703703703702d0*(1/(-1.d0+n)-
     &  1/n-2.d0/(1.d0+n))+1.1851851851851851d1*(S1/n-(-(1/n)+S1
     &  )/(-1.d0+n)-(8.d-1*(1/(1.d0+n)+S1))/(1.d0+n))+3.55555555
     &  55555554d0*((-n**(-2)+(-(1/n)+S1)**2+S2)/(-1.d0+n)+(5.d-
     &  1*((1.d0+n)**(-2)+(1/(1.d0+n)+S1)**2+S2))/(1.d0+n)-(S2+S
     &  1*S1)/n))-4.938271604938271d0*(S1**4/n+(8.d0*S1*S3)/n+(6
     &  .d0*S4)/n+(6.d0*S2*(S1*S1))/n+(3.d0*(S2*S2))/n))
*
       RETURN
       END
