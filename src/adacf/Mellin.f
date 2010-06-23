C     According to
C        Eq. (33) in Bluemlein and Kurth, hep-ph/9708388
     
      DOUBLE COMPLEX FUNCTION MELLINF2 (N)
*
      IMPLICIT NONE
      INTEGER K
      DOUBLE COMPLEX N, PSI, PSITMP, MF2
      DOUBLE PRECISION LOG2, EMC, ZETA2, ABK(8)
      PARAMETER ( LOG2 = 0.69314 71806D0, EMC  = 0.57721 56649D0,
     &            ZETA2 = 1.64493 40668D0 )
      DATA (ABK(K), K=1,8) / 0.9999964239d0, -0.4998741238d0, 
     &   0.3317990258d0, -0.2407338084d0, 0.1676540711d0, 
     &  -0.0953293897d0, 0.0360884937d0, -0.0064535442d0 /
*
*
      PSITMP = PSI(N)
      MF2 = (0.0d0, 0.0d0)
*
      DO 10 K = 1, 8
         PSITMP = PSITMP + 1.0d0 / (N + K - 1)
         MF2 = MF2 + ABK(K) * ( (N - 1) * ( ZETA2 / (N + K - 1) -
     &         (PSITMP + EMC) / (N + K - 1)**2 ) +
     &         (PSITMP + EMC) / (N + K - 1) )
 10   CONTINUE
*
      MELLINF2 =  ZETA2 * LOG2  -  MF2
*
      RETURN
      END
