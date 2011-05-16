*  The code in this file is mostly taken from the source of paper
*        A. Vogt,  hep-ph/0408244
*
* =================================================================av==
*
* .. Harmonic sums ..
*
* =================================================================av==

       DOUBLE COMPLEX FUNCTION HS1 (Z)
*
       DOUBLE COMPLEX Z, PSI
       DOUBLE PRECISION EMC, ZETA2, ZETA3, ZETA4
       PARAMETER ( EMC  = 0.57721 56649D0, ZETA2 = 1.64493 40668D0,
     ,            ZETA3 = 1.20205 69032D0, ZETA4 = 1.08232 32337D0 )
*
*
       HS1 =  EMC  + PSI (Z+1.) 
*
       RETURN
       END

* =================================================================av==

       DOUBLE COMPLEX FUNCTION HS2 (Z)
*
       DOUBLE COMPLEX Z, DPSI
       DOUBLE PRECISION EMC, ZETA2, ZETA3, ZETA4
       PARAMETER ( EMC  = 0.57721 56649D0, ZETA2 = 1.64493 40668D0,
     ,            ZETA3 = 1.20205 69032D0, ZETA4 = 1.08232 32337D0 )
*

*
       HS2 = ZETA2 - DPSI (Z+1.,1)

       RETURN
       END

* =================================================================av==

       DOUBLE COMPLEX FUNCTION HS3 (Z)
*
       DOUBLE COMPLEX Z, DPSI
       DOUBLE PRECISION EMC, ZETA2, ZETA3, ZETA4
       PARAMETER ( EMC  = 0.57721 56649D0, ZETA2 = 1.64493 40668D0,
     ,            ZETA3 = 1.20205 69032D0, ZETA4 = 1.08232 32337D0 )
*

*
       HS3 = ZETA3 + 0.5 * DPSI (Z+1.,2)

       RETURN
       END

* =================================================================av==


       DOUBLE COMPLEX FUNCTION HS4 (Z)
*
       DOUBLE COMPLEX Z, DPSI
       DOUBLE PRECISION EMC, ZETA2, ZETA3, ZETA4
       PARAMETER ( EMC  = 0.57721 56649D0, ZETA2 = 1.64493 40668D0,
     ,            ZETA3 = 1.20205 69032D0, ZETA4 = 1.08232 32337D0 )
*

*
       HS4 = ZETA4 - 1./6.D0 * DPSI (Z+1.,3)

       RETURN
       END

* =================================================================kk==

       DOUBLE COMPLEX FUNCTION HSM1 (Z, K)
*
      IMPLICIT NONE
       INTEGER K
       DOUBLE COMPLEX Z, PSI
       DOUBLE PRECISION EMC, LOG2
       PARAMETER ( EMC  = 0.57721 56649D0, LOG2 = 0.69314 718056D0)
*
       
        HSM1 = - LOG2 + 0.5d0 * (1-2*MOD(K,2)) * ( PSI((Z+2.0d0)/2.0d0)
     &    -PSI((Z+1.d0)/2.0d0) )

       RETURN
       END

* =================================================================kk==

       DOUBLE COMPLEX FUNCTION HSM2 (Z, K)
*
      IMPLICIT NONE
       INTEGER K
       DOUBLE COMPLEX Z, DPSI
       DOUBLE PRECISION EMC, ZETA2, ZETA3, ZETA4
       PARAMETER ( EMC  = 0.57721 56649D0, ZETA2 = 1.64493 40668D0,
     ,            ZETA3 = 1.20205 69032D0, ZETA4 = 1.08232 32337D0 )
*
       
        HSM2 = - ZETA2/2.0d0 + 0.25d0 * (1-2*MOD(K+1,2)) * ( 
     &        DPSI((Z+2.0d0)/2.0d0, 1) - DPSI((Z+1.d0)/2.0d0, 1) )
              

       RETURN
       END

* =================================================================kk==

       DOUBLE COMPLEX FUNCTION HSM3 (Z, K)
*
      IMPLICIT NONE
       INTEGER K
       DOUBLE COMPLEX Z, DPSI
       DOUBLE PRECISION EMC, ZETA2, ZETA3, ZETA4
       PARAMETER ( EMC  = 0.57721 56649D0, ZETA2 = 1.64493 40668D0,
     ,            ZETA3 = 1.20205 69032D0, ZETA4 = 1.08232 32337D0 )
*
       
        HSM3 = - 0.75d0*ZETA3 + (1-2*MOD(K,2)) * ( 
     &    DPSI((Z+2.0d0)/2.0d0, 2) - DPSI((Z+1.d0)/2.0d0, 2) )/16.0d0
              

       RETURN
       END

* =================================================================kk==

      DOUBLE COMPLEX FUNCTION MDILOG (J)

*      Approximate Mellin moment x^J of dilog(x)/(1+x) 
*      (about 1e-3 accuracy)

      IMPLICIT NONE
      DOUBLE COMPLEX J
      INTEGER I, IMAX
      PARAMETER (IMAX=5)
      DOUBLE PRECISION COEF(IMAX)
       
*     coefs of odd powers of x
*     dilog(x)/(1+x) = c1*x + c2*x^3 + ... + ck*x^(2*k-1)
      DATA (COEF(I), I = 1, IMAX)
     & /  0.9191873351383635 D0, -1.1104046996355397 D0,
     &    2.895284481438232 D0, - 3.5219362238384186 D0,
     &    1.640336140321476 D0 /
*
       MDILOG = 0.0d0
       DO 10 I = 1, IMAX
 10      MDILOG = MDILOG + COEF(I)/(J + 2*I)

       RETURN
       END

* =================================================================kk==

       DOUBLE COMPLEX FUNCTION HSM21 (Z, K)
*
       IMPLICIT NONE
       INTEGER K
       DOUBLE COMPLEX Z, MDILOG, HSM1
       DOUBLE PRECISION EMC, ZETA2, ZETA3, LOG2
       EXTERNAL HSM1, MDILOG
       PARAMETER ( EMC  = 0.57721 56649D0, ZETA2 = 1.64493 40668D0,
     ,            ZETA3 = 1.20205 69032D0, LOG2 = 0.69314 718056D0 )
*

        HSM21 = ZETA2*HSM1(Z,K) - 5.0d0*ZETA3/8.0d0 + ZETA2*LOG2
     &     + (1-2*MOD(K+1,2)) * MDILOG(Z)
              

       RETURN
       END

* =================================================================kk==

       DOUBLE COMPLEX FUNCTION HS1M2 (Z, K)
*
       INTEGER K
       DOUBLE COMPLEX Z, HSM21, HSM2, HSM3, HS1
       EXTERNAL HSM21, HSM2, HSM3, HS1
*

       HS1M2 = - HSM21(Z, K) + HSM2(Z, K) * HS1(Z) + HSM3(Z, K)

       RETURN
       END



* =================================================================av==
*
* =====================================================================
* ..The complex psi function,  PZI(Z),  and its m'th derivatives, 
*    DPZI(Z,M),  calculated from the asymtotic expansions. The 
*    functional equations are used for  |Im(Z)| < 10  to shift the 
*    argument to  Re(Z) >= 10  before applying the expansions.
* =====================================================================

       FUNCTION PSI (Z)
*
       IMPLICIT DOUBLE COMPLEX (A-Z)
       SUB = 0.D0
       ZZ = Z
*
* ..Shift of the argument using the functional equation
*
       IF (DABS(DIMAG(ZZ)) .LT. 10.D0) THEN
*
  1      CONTINUE
         IF (DBLE(ZZ) .LT. 10.D0) THEN
           SUB = SUB - 1./ ZZ
           ZZ = ZZ + 1.
           GOTO 1
         END IF
*
       END IF
*
* ..Use of the asymtotic expansion (at the shifted argument)
*
       RZ = 1./ ZZ
       DZ = RZ * RZ
       PSI = SUB + LOG(ZZ) - 0.5 * RZ - DZ/5040.D0 * ( 420.+ DZ *
     1       ( - 42. + DZ * (20. - 21. * DZ) ) )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
       FUNCTION DPSI (Z,M)
*
       IMPLICIT DOUBLE COMPLEX (A-Z)
       INTEGER M, K1, K2
       SUB = 0.D0
       ZZ = Z
*
* ..Shift of the argument using the functional equations
*
       IF (DABS(DIMAG(ZZ)) .LT. 10.D0) THEN
*
  1      CONTINUE
         SUBM = -1./ZZ
         DO 10 K1 = 1, M
           SUBM = - SUBM * K1 / ZZ
 10      CONTINUE
*
         IF (DBLE(ZZ) .LT. 10.D0) THEN
           SUB = SUB + SUBM
           ZZ = ZZ + 1.
           GOTO 1
         END IF
*
       END IF
*
* ..Expansion coefficients for the first derivative
*
       A1 =  1.D0
       A2 =  1./2.D0
       A3 =  1./6.D0
       A4 = -1./30.D0
       A5 =  1./42.D0
       A6 = -1./30.D0
       A7 =  5./66.D0
*
* ..Expansion coefficients for the higher derivatives
*
       IF (M .EQ. 1) GO TO 2
       DO 11 K2 = 2, M
         A1 = A1 * (K2-1.)
         A2 = A2 *  K2
         A3 = A3 * (K2+1.)
         A4 = A4 * (K2+3.)
         A5 = A5 * (K2+5.)
         A6 = A6 * (K2+7.)
         A7 = A7 * (K2+9.)
  11   CONTINUE
  2    CONTINUE 
*
* ..Use of the asymtotic expansion (at the shifted argument)
*
       RZ = 1./ ZZ
       DZ = RZ * RZ
       DPSI = SUB + (-1)**(M+1) * RZ**M * ( A1 + RZ * (A2 + RZ * 
     1        (A3 + DZ * (A4 + DZ * (A5 + DZ * (A6 + A7 * DZ ))))) )
*
       RETURN
       END
*
* =================================================================av==
