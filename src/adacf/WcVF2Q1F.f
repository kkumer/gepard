*    
       SUBROUTINE WcVF2Q1F (NF, n, res)
*
       IMPLICIT none
       DOUBLE COMPLEX n, res
       DOUBLE COMPLEX S1, S2, S3, S4
       DOUBLE PRECISION NF, CF, CA, TF
       PARAMETER ( CF=4./3.d0, CA=3.d0, TF=0.5d0 )
*
*      ... Harmonic sum initialization ...
*
      COMMON / HARMONIC /  S1, S2, S3, S4
*
*     ... Function itself spliced in by Mathematica ...
               res=5.d-1*CF*(-9.d0+2.d0/n**2+3.d0/n+4.d0/(1.d0+n)+(3.d0
     &  -2.d0/(n*(1.d0+n)))*S1-2.d0*S2+2.d0*(S1*S1))
*
       RETURN
       END
