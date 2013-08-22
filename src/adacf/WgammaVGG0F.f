*    
       SUBROUTINE WgammaVGG0F (NF, n, res)
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
               res=5.d-1*(-7.333333333333333d0*CA-8.d0*CA*(1/((-1.d0+n)
     &  *n)+1/((1.d0+n)*(2.d0+n))-S1)+2.6666666666666665d0*NF*TF
     &  )
*
       RETURN
       END
