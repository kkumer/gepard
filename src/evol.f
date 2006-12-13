C     ****h* gepard/evol.f
C  FILE DESCRIPTION
C     calculates evolution operator
C
C    $Id$
C     *******

C     ****s* evol.f/EVOLNSF
C  NAME
C     EVOLNSF  --  non-singlet evolution operator
C  DESCRIPTION
C    calculates non-singlet evolution operator
C  SYNOPSIS
C     SUBROUTINE EVOLNSF (K, R, EVOLNSA)
C
C     INTEGER K
C     DOUBLE PRECISION R
C     DOUBLE COMPLEX EVOLNSA(0:2), 
C  INPUTS
C           K -- Mellin-Barnes integration point index
C           R -- ratio of astrong(mu)/astrong(mu0)
C  OUTPUT
C       EVOLNSA -- three numbers representing three coefficients 
C                 {\cal A}^{0,1,2}(a_s/a0_s)^{-\gamma^0/beta_0} in
C                 alpha_s/(2\pi) expansion of evolution operator 
C                 \Epsilon(mu, mu0)
C  PARENTS
C      PARWAVF
C  CHILDREN
C      ERFUNCNSF, RNNLONSF
C  SOURCE
C

      SUBROUTINE EVOLNSF (K, R, EVOLNSA)

      IMPLICIT NONE
      INTEGER K
      DOUBLE PRECISION R
      DOUBLE COMPLEX EVOLNSA(0:2), NDINT
      INTEGER ORD
      DOUBLE COMPLEX R1, R2, AUX(0:2)
      INCLUDE 'header.f'


      CALL RNNLONSF (K, R1, R2)

      
      AUX(0) = (1.0d0, 0.0d0)

      IF (SCHEME(4:5) .EQ. 'LO') THEN 
*       Just LO evolution ...
        AUX(1) = (0.0d0, 0.0d0)
        AUX(2) = (0.0d0, 0.0d0)
      ELSE
*       ... or normal (N)NLO one (diagonal part)
        AUX(1) = - ( 1.0d0 - 1.0d0 / R ) * R1
        IF (P .GE. 2) THEN
          AUX(2) = 0.5d0 * ( AUX(1)**2 - (1.0d0 - 1.0d0 / R**2) * R2 )
        END IF
      ENDIF

      DO 10 ORD = 0, P
 10   EVOLNSA(ORD) = AUX(ORD) * R**( - GAMNS(K,0) / BETA0(NF) )

*   Adding MSBAR non-diagonal NLO evolution (time consuming)
      IF ( (P .EQ. 1) .AND. (SCHEME .EQ. 'MSBND') ) THEN
        CALL NDINTF(K, R, NDINT, 0, 0)
        EVOLNSA(1) = EVOLNSA(1) + NDINT
      ENDIF


      RETURN
      END
C     *****


C     ****s* evol.f/EVOLF
C  NAME
C     EVOLF  --   singlet evolution operator
C  DESCRIPTION
C    calculates singlet evolution operator
C    according to KMKPS06 paper
C  SYNOPSIS
C     SUBROUTINE EVOLF (K, R, EVOLA)
C
C     INTEGER K
C     DOUBLE PRECISION R
C     DOUBLE COMPLEX EVOLA(0:2,2,2), 
C  INPUTS
C           K -- Mellin-Barnes integration point index
C           R -- ratio of astrong(mu)/astrong(mu0)
C  OUTPUT
C       EVOLA -- three 2x2 matrices representing three coefficients 
C                 {\cal A}^{0,1,2}(a_s/a0_s)^{-\gamma^0/beta_0} in
C                 alpha_s/(2\pi) expansion of evolution operator 
C                 \Epsilon(mu, mu0)
C  PARENTS
C      PARWAVF
C  CHILDREN
C      ERFUNCF, RNNLOF, KRONECKER
C  SOURCE
C

      SUBROUTINE EVOLF (K, R, EVOLA)

      IMPLICIT NONE
      INTEGER K
      DOUBLE PRECISION R
      DOUBLE COMPLEX EVOLA(0:2,2,2)
      INTEGER A, B, CE, I, J, K1, ORD
      DOUBLE PRECISION RINV
      DOUBLE COMPLEX LAM(2), PR(2,2,2)
      DOUBLE COMPLEX ERFUNC1(2,2), ERFUNC2(2,2)
      DOUBLE COMPLEX R1PROJ(2,2,2,2), R2PROJ(2,2,2,2)
      DOUBLE COMPLEX AUX, DINV, IJAUX(0:2), R1R1(2,2), NDINT
      INTEGER KRONECKER
      INCLUDE 'header.f'

      CALL LAMBDAF(K, LAM)
      CALL ERFUNCF (R, LAM, LAM, ERFUNC1, ERFUNC2)
      CALL RNNLOF (K, LAM, PR, R1PROJ, R2PROJ)

      RINV = 1.0d0 / R

      DO 10 I = 1, 2
      DO 10 J = 1, 2
      DO 10 ORD = 0, P
 10      EVOLA(ORD,I,J) = (0.0d0, 0.0d0)

      DO 40 I = 1, 2
      DO 40 J = 1, 2

      DO 40 B = 1, 2
      DO 15 ORD = 0, P
 15      IJAUX(ORD) = (0.0d0, 0.0d0)
      DO 30 A = 1, 2
         IJAUX(0) = IJAUX(0) + KRONECKER(A,B) * PR(A,I,J)
         IJAUX(1) = IJAUX(1) - ERFUNC1(A,B) * R1PROJ(A,B,I,J)
      IF (P. GE. 2) THEN
         AUX = (0.0d0, 0.0d0)
         DO 20 CE = 1, 2
            DINV = 1.0d0 + (LAM(CE) - LAM(B)) / BETA0(NF)
            R1R1(I,J) = (0.0d0, 0.0d0)
            DO 17 K1 = 1, 2
 17         R1R1(I,J) = R1R1(I,J) + R1PROJ(A,CE,I,K1)*R1PROJ(CE,B,K1,J)
 20         AUX = AUX + (ERFUNC2(A,B) - ERFUNC1(A,CE) * RINV**DINV)/DINV
     &            * R1R1(I,J)
         IJAUX(2) = IJAUX(2) + AUX - ERFUNC2(A,B) * R2PROJ(A,B,I,J)
      END IF
 30   CONTINUE
      DO 35 ORD = 0, P
        IF (SCHEME(4:5) .EQ. 'LO') THEN 
*         Just LO evolution ...
          IJAUX(1) = (0.0d0, 0.0d0)
          IJAUX(2) = (0.0d0, 0.0d0)
        ENDIF
 35   EVOLA(ORD,I,J) = EVOLA(ORD,I,J) 
     &                 + IJAUX(ORD) * R**( -LAM(B) / BETA0(NF) )
 40   CONTINUE

*   Adding MSBAR non-diagonal NLO evolution (time consuming)
      IF ( (P .EQ. 1) .AND. (SCHEME .EQ. 'MSBND') ) THEN
        DO 50 I = 1, 2
        DO 50 J = 1, 2
          CALL NDINTF(K, R, NDINT, I, J)
          EVOLA(1,I,J) = EVOLA(1,I,J) + NDINT
 50     CONTINUE
      ENDIF

      RETURN
      END
C     *****


C     ****f* evol.f/KRONECKER
C  NAME
C    KRONECKER -- kronecker delta symbol
C  SYNOPSIS
C    INTEGER FUNCTION KRONECKER (A, B)
C
C    INTEGER A, B
C  PARENTS
C       EVOLF
C  SOURCE
C

      INTEGER FUNCTION KRONECKER (A, B)

      IMPLICIT NONE
      INTEGER A, B

      IF (A .EQ. B) THEN
            KRONECKER = 1
      ELSE
            KRONECKER = 0
      END IF

      RETURN
      END
C     *****
