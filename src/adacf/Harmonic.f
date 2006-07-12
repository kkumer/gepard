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
