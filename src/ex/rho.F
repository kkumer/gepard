C     ****h* gepard/rho.F
C  FILE DESCRIPTION
C    Simple program for testing DVEM - rho production
C
C    $Id:$
C     *******

C     ****p* rho.F/RHO
C  NAME
C    RHO  --  DVrhoP for low-t H1 kinematics
C            
C  DESCRIPTION
C             For tests and debugging DVEM
C
C  CHILDREN
C      READPAR, INIT, TFF
C
C  SOURCE
C


      PROGRAM TEST

      IMPLICIT NONE
      INTEGER K
      DOUBLE COMPLEX CAUX
      DOUBLE PRECISION XS
#include "header.F"
      DOUBLE PRECISION Q2EXP(5), XSEXP(5)
      DATA Q2EXP / 3.3d0, 6.6d0, 11.5d0, 17.4d0, 33.d0 /
      DATA XSEXP / 2156.d0, 604.d0, 181.d0, 51.1d0, 6.5d0 /

      FFTYPE  = 'SINGLET'

      CALL READPAR
*   modifications of GEPARD.INI
      SPEED = 3
      NF = 4
      Q02 = 4.0
      ASP(0) = 0.0606
      ASP(1) = 0.0518
      ASP(2) = 0.0488
      C = 0.35
      PHI = 1.57079632
      CND = -0.25
      PHIND = 1.57
      ANSATZ = 'FIT'
      SCHEME = 'MSBAR'
      PROCESS  = 'DVMP'


      RGPDF2 = 1.d0
      RDAF2 = 1.d0
      RR2 = 1.d0

*   GPD ansatz KM10b

      PAR(11) = 0.152039d0
      PAR(12) = 1.1575506d0
      PAR(13) = 0.15d0
      PAR(14) = 0.4975d0
      PAR(15) = 0.0d0
      PAR(16) = 2.0d0

      PAR(22) = 1.247317d0
      PAR(23) = 0.15d0
      PAR(24) = 0.7d0
      PAR(25) = 0.0d0
      PAR(26) = 2.0d0

*  second PW
      PAR(17) = -0.46d0
      PAR(27) = -2.515d0

c     PAR(17) = 0
c     PAR(27) = 0

*  third PW
      PAR(37) = 0.0935d0
      PAR(47) = 0.89158d0

c     PAR(37) = 0
c     PAR(47) = 0

*   turn off GPD E
      PAR(18) = 0.0d0
      PAR(28) = 0.0d0

*   turn off xi-dependence
      PAR(19) = 0.0d0
      PAR(29) = 0.0d0



      DEL2 = -0.025d0

      WRITE (*, *) " ---  Test 1: TFF (no evolution) ---- "
      P = 0
      CALL INIT
      W2 = 75.0d0**2
      Q2 = 4.0d0
      NQS = 1
      QS(1) = Q2
      CALL EVOLC(1,1)
      CALL EVOLC(2,1)
      CALL EVOLC(3,1)
      XI = Q2 / ( 2.0d0 * W2 + Q2)
      MTIND = 0
      CALL TFFF
      CAUX = TFFH(0)
      WRITE (*, *) CAUX

      P = 1
      CALL INIT

      WRITE (*, *) " ---  Test 2: NLO TFF (no evolution) ---- "
      W2 = 75.0d0**2
      Q2 = 4.0d0
      NQS = 1
      QS(1) = Q2
      CALL EVOLC(1,1)
      CALL EVOLC(2,1)
      CALL EVOLC(3,1)
      XI = Q2 / ( 2.0d0 * W2 + Q2)
      MTIND = 0
      CALL TFFF
      CAUX = TFFH(1)
      WRITE (*, *) CAUX

      WRITE (*, *)
      WRITE (*, *) " ---  Test 3: dsigma_rho_L / dt ---- "
      P = 0
      CALL INIT
      DO 10 K=1,5
        Q2 = Q2EXP(K)
        NQS = 1
        QS(1) = Q2
        CALL EVOLC(1,1)
        CALL EVOLC(2,1)
        CALL EVOLC(3,1)
        XI = Q2 / ( 2.0d0 * W2 + Q2)
        CALL TFFF
        XS = 112175.5 * (2.d0*XI)**2 
     &       *  ABS(TFFH(P))**2 / Q2**2
c     WRITE (*,*) TFFH(P)
10    WRITE (*, 999) Q2EXP(K), XS, XSEXP(K)

      WRITE (*, *)
      WRITE (*, *) " ---  Test 4: NLO dsigma_rho_L / dt ---- "
      P = 1
      CALL INIT
      DO 20 K=1,5
        Q2 = Q2EXP(K)
        NQS = 1
        QS(1) = Q2
        CALL EVOLC(1,1)
        CALL EVOLC(2,1)
        CALL EVOLC(3,1)
        XI = Q2 / ( 2.0d0 * W2 + Q2)
        CALL TFFF
        XS = 112175.5 * (2.d0*XI)**2 
     &       *  ABS(TFFH(P))**2 / Q2**2
c     WRITE (*,*) TFFH(P)
20    WRITE (*, 999) Q2EXP(K), XS, XSEXP(K)

999   FORMAT (1X,'Q2 = ', F4.1, 5X, 'XS = ', F6.1, 5X, 'vs. ', F6.1)
      STOP
      END
C     ****
