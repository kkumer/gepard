*    
       SUBROUTINE WcVFLG2F (NF, n, res)
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
               res=2.5d-1*NF*(-5.333d0/(-1.d0+n)-3.9659999999999993d1/n
     &  **2+5.333d0/n+2.1542399999999997d3/(1.d0+n)**3+1.00286d3
     &  /(1.d0+n)**2-1.9097680000000001d3/(1.d0+n)+9.84000000000
     &  0002d1/(2.d0+n)**3-9.840000000000002d1/(2.d0+n)**2+(-8.6
     &  48d2/n+8.7312d2/(1.d0+n)**2+9.632000000000001d2/(1.d0+n)
     &  +9.840000000000002d1/(2.d0+n)**2-9.840000000000002d1/(2.
     &  d0+n))*S1+(9.473999999999998d1/n+1.0170599999999999d3/(1
     &  .d0+n)+4.920000000000001d1/(2.d0+n))*S2+(9.4739999999999
     &  98d1/n-1.4393999999999998d2/(1.d0+n)+4.920000000000001d1
     &  /(2.d0+n))*(S1*S1))
*
       RETURN
       END
