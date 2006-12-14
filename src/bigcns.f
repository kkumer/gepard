C     ****h* gepard/bigcns.f
C  FILE DESCRIPTION
C    calculation of  Wilson coefficients for DVCS in CSbar scheme
C    according to KMKPS06 paper
C
C    $Id$
C     *******


C     ****s* bigcns.f/BIGCNSF
C  NAME
C     BIGCNSF  --  "big C" Wilson coefficients (non-singlet)
C  DESCRIPTION
C    calculates Wilson coefficients for DVCS in CSbar scheme
C    according to KMKPS06 paper.
C  SYNOPSIS
C     SUBROUTINE BIGCNSF (K)
C
C  INPUTS
C           J -- conformal moment
C  OUTPUT
C       BIGC0 -- vector of two C^{(0)} (quark and gluon) Wilson 
C                coefficients (trivially equal to (1,0))
C       BIGC1 -- vector  C^{(1)} 
C       BIGC2 -- vector  C^{(2)} 
C  IDENTIFIERS
C       BETABLK, WGAMMA, WC -- common blocks with beta function coefficients
C       of QCD, and J-th moments of anomalous dimensions and Wilson 
C       coefficients of DIS, initialized by subroutine INIT
C  PARENTS
C      INIT
C  CHILDREN
C      VECMAT, HS1, HS2
C  SOURCE
C

      SUBROUTINE BIGCNSF ( K )

      IMPLICIT NONE
      INTEGER K
      DOUBLE PRECISION LRF2, LRR2
      DOUBLE COMPLEX J, HS1, HS2
      DOUBLE COMPLEX SHIFT1, SHIFT2
      INCLUDE 'header.f'

      J = N(K) - 1

      LRF2 = LOG(RF2)
      LRR2 = LOG(RR2)

      IF ( PROCESS(:3) .EQ. 'DVC' ) THEN
        SHIFT1 = HS1(J + 1.5d0) - HS1(J + 2.0d0) + 2.0d0 * LOG(2.0d0)
     &        - LRF2
        SHIFT2 = SHIFT1*SHIFT1 - HS2(J + 1.5d0) + HS2(J + 2.0d0)
      ELSE
        SHIFT1 =  - LRF2
        SHIFT2 = SHIFT1*SHIFT1
      END IF

      BIGCNS(K, 0) = (1.0d0, 0.0d0)

      IF ( SCHEME .EQ. 'CSBAR' ) THEN

*      -CSBAR-
      IF (P .GE. 1) THEN
        BIGCNS(K, 1) = CDISNS1(K, 1) + 0.5d0 * SHIFT1 * CDISNS1(K, 0) *
     &                  GAMNS(K, 0)
        IF (P .GE. 2) THEN
          BIGCNS(K, 2) = CDISNS1(K, 2) + 0.5d0 * SHIFT1 * (CDISNS1(K, 0)
     &               * GAMNS(K, 1) + CDISNS1(K, 1) * GAMNS(K, 0)) +
     &      0.125d0 * SHIFT2 * CDISNS1(K, 0) * GAMNS(K, 0)**2 + 
     &      0.5d0 * BETA0(NF) * ( BIGCNS(K, 1) * LRR2 + 0.25d0 * 
     &        CDISNS1(K, 0) * GAMNS(K, 0) * LRF2**2 )
        END IF
      END IF

      ELSE
              
*      -MSBAR-
        BIGCNS(K, 1) =  CF * (2.0d0 * HS1(J+1)**2 - (9.0d0/2.0d0)
     &     + (5.0d0 - 4.0d0 * HS1(J+1)) / (2.0d0 * (J+1) * (J+2))
     &     + 1.0d0 / ( (J+1)**2 * (J+2)**2 ) ) - GAMNS(K,0)*LRF2/ 2.0d0

      END IF

      RETURN
      END
C     *******

