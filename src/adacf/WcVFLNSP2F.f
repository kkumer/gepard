*    
       SUBROUTINE WcVFLNSP2F (NF, n, res)
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
               res=2.5d-1*(-1.6400000000000001d-1-6.2d-2/n**3+1.4849999
     &  999999999d1/n**2+2.07153d2/n+5.312d1/(1.d0+n)**3+9.748d1
     &  /(1.d0+n)+(-1.505d2/n**2+5.579d1/n)*S1+5.925925925925926
     &  d-1*NF*(6.d0/n+6.d0/(1.d0+n)**2-2.5d1/(1.d0+n)-(6.d0*S1)
     &  /(1.d0+n))-(1.3688d2*S2)/n+(1.3619999999999999d1*(S1*S1)
     &  )/n)
*
       RETURN
       END
