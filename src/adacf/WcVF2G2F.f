*    
       SUBROUTINE WcVF2G2F (NF, n, res)
*
       IMPLICIT none
       DOUBLE COMPLEX n, res
       DOUBLE COMPLEX S1, S2, S3, S4
       DOUBLE PRECISION NF, CF, CG, TF
       PARAMETER ( CF=4./3.d0, CG=3.d0, TF=0.5d0 )
*
*      ... Harmonic sum initialization ...
*
      COMMON / HARMONIC /  S1, S2, S3, S4
*
*     ... Function itself spliced in by Mathematica ...
               res=2.5d-1*NF*(-2.71d-1+1.5059d3/(-1.d0+n)-3.19140000000
     &  00002d1/n**4-1.1896d2/n**3-2.0974d3/n**2-4.93834d3/n+1.2
     &  564d3/(1.d0+n)**4-(1.494d3/(-1.d0+n)-1.4482d3/n**3+1.385
     &  12d3/n-1.2564d3/(1.d0+n)**3)*S1-(2.15845d2/n-2.094d2/(1.
     &  d0+n))*S1**3+(2.3200000000000003d3/n**2-2.4d1/n+6.282000
     &  000000001d2/(1.d0+n)**2)*S2+(1.0960699999999999d3/n+6.28
     &  2000000000001d2/(1.d0+n))*S1*S2+(2.76011d3/n+4.188d2/(1.
     &  d0+n))*S3+(8.718d2/n**2-2.4d1/n+6.282000000000001d2/(1.d
     &  0+n)**2)*(S1*S1))
*
       RETURN
       END
