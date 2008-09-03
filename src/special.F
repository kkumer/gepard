C     ****h* gepard/special.f
C  FILE DESCRIPTION
C    Double complex special functions (Gamma and Beta)
C
C    $Id$
C  NOTES
C     - Using algorithm of P. Godfrey from http://my.fit.edu/~gabdo/gamma.txt
C     *******


C     ****f* special.f/CLNGAMMA
C  NAME
C     CLNGAMMA  --  logarithm of gamma function of complex argument
C  DESCRIPTION
C     Uses Lanczos algorithm
C  SYNOPSIS

      DOUBLE COMPLEX FUNCTION CLNGAMMA(Z)

      IMPLICIT NONE
      DOUBLE COMPLEX Z

C  INPUTS
C     Z -- argument
C  RETURN VALUE
C     CLNGAMMA -- ln(Gamma(Z))
C  NOTES
C     Error for P=1 is strictly smaller than 2x10^(-10), and in practice is often
C     less then 10^(-14). This is original Lanczos parameterization.
C  PARENTS
C     INIT, CBETA, HJ
C  SOURCE
C

      INTEGER SPD
#include "header.F"
      INTEGER K
      INTEGER NMAX, NP
      PARAMETER (NMAX = 7, NP = 4)
      INTEGER NN(NP)
      DOUBLE PRECISION G(NP), COEF(NP, 0:NMAX-1)
      DOUBLE COMPLEX TMP, X, YY, S
      SAVE NN, G, COEF
      DATA NN            /7, 3, 2, 1/
      DATA G            /5.0d0, 2.60404494991116d0, 
     &                   1.49614999d0, 0.328498d0/
      DATA COEF(4,0)     /1.8113464d0/
      DATA COEF(3,0), COEF(3,1) /0.5613831815d0,0.6055625385d0/
      DATA COEF(2,0), COEF(2,1), COEF(2,2) /
     &                    0.1854248319548753,0.8797922715613486,
     &                                      -0.2588333483439386/
      DATA COEF(1,0), COEF(1,1), COEF(1,2), COEF(1,3), COEF(1,4), 
     &     COEF(1,5), COEF(1,6)     /0.0168895284640819927d0, 
     &                   1.2866458274168037d0,
     &                  -1.46103406972059703d0, 
     &                   0.405586795700707467d0,
     &                  -0.0208035005652801098, 
     &                   0.0000204135450223743664,
     &                  -9.11230491453873551e-8/


      SPD = 3
      IF (SPEED .EQ. 1) SPD = 1

      X = Z
      YY = X - (1.d0, 0.0d0)
      TMP = Z - (0.5d0, 0.0d0)
      TMP = TMP * LOG( (TMP + G(SPD)) ) - TMP
      S = COEF(SPD, 0)
      DO 10 K = 1, NN(SPD)-1
      YY = YY + (1.d0, 0.0d0)
 10   S = S + COEF(SPD, K) / YY

      CLNGAMMA = TMP + LOG(S)

*     Taking principal value of logarithm: -Pi < Im(log()) <= Pi
      DO WHILE ( IMAGPART(CLNGAMMA) .GT. PI )
        CLNGAMMA = CLNGAMMA - TWOPII
      END DO
      DO WHILE ( IMAGPART(CLNGAMMA) .LE. (-PI) )
        CLNGAMMA = CLNGAMMA + TWOPII
      END DO

      RETURN
      END
C     ***

C     ****f* special.f/MULTIGAMMA
C  NAME
C     MULTIGAMMA  --  Ratio of two gamma functions
C  DESCRIPTION
C     Uses Lanczos algorithm
C  SYNOPSIS

      DOUBLE COMPLEX FUNCTION MULTIGAMMA(ZN1,ZN2,ZN3,ZD1,ZD2,ZD3,ZDSQ)

      IMPLICIT NONE
      DOUBLE COMPLEX ZN1,ZN2,ZN3,ZD1,ZD2,ZD3,ZDSQ

C  INPUTS
C     Z... -- arguments
C  RETURN VALUE
C     MULTIGAMMA -- Gamma(ZN1)Gamma(ZN2)Gamma(ZN3)/(
C           Gamma(ZD1)Gamma(ZD2)Gamma(ZD3)Gamma(ZDSQ))
C  PARENTS
C     HJ
C  SOURCE
C

      INTEGER SPD
#include "header.F"
      INTEGER K, L
      INTEGER NMAX, NP
      PARAMETER (NMAX = 7, NP = 4)
      INTEGER NN(NP)
      DOUBLE PRECISION INVEE, G(NP), COEF(NP, 0:NMAX-1)
      PARAMETER (INVEE = 0.3678794411714423215955d0)
      DOUBLE COMPLEX ZN(3), ZD(3), TMP, X, YY, S
      SAVE NN, G, COEF
      DATA NN            /7, 3, 2, 1/
      DATA G            /5.0d0, 2.60404494991116d0, 
     &                   1.49614999d0, 0.328498d0/
      DATA COEF(4,0)     /1.8113464d0/
      DATA COEF(3,0), COEF(3,1) /0.5613831815d0,0.6055625385d0/
      DATA COEF(2,0), COEF(2,1), COEF(2,2) /
     &                    0.1854248319548753,0.8797922715613486,
     &                                      -0.2588333483439386/
      DATA COEF(1,0), COEF(1,1), COEF(1,2), COEF(1,3), COEF(1,4), 
     &     COEF(1,5), COEF(1,6)     /0.0168895284640819927d0, 
     &                   1.2866458274168037d0,
     &                  -1.46103406972059703d0, 
     &                   0.405586795700707467d0,
     &                  -0.0208035005652801098, 
     &                   0.0000204135450223743664,
     &                  -9.11230491453873551e-8/


      ZN(1) = ZN1
      ZN(2) = ZN2
      ZN(3) = ZN3

      ZD(1) = ZD1
      ZD(2) = ZD2
      ZD(3) = ZD3

      SPD = 3
      IF (SPEED .EQ. 1) SPD = 1

      MULTIGAMMA = (1.0d0, 0.0d0)

      DO 40 L = 1, 3
      X = ZN(L)
      YY = X - (1.d0, 0.0d0)
      TMP = ZN(L) - (0.5d0, 0.0d0)
      TMP = ( (TMP + G(SPD)) * INVEE )**TMP
      S = COEF(SPD, 0)
      DO 20 K = 1, NN(SPD)-1
      YY = YY + (1.d0, 0.0d0)
 20   S = S + COEF(SPD, K) / YY
      MULTIGAMMA = MULTIGAMMA * (TMP * S)

      X = ZD(L)
      YY = X - (1.d0, 0.0d0)
      TMP = ZD(L) - (0.5d0, 0.0d0)
      TMP = ( (TMP + G(SPD)) * INVEE )**TMP
      S = COEF(SPD, 0)
      DO 30 K = 1, NN(SPD)-1
      YY = YY + (1.d0, 0.0d0)
 30   S = S + COEF(SPD, K) / YY
      MULTIGAMMA = MULTIGAMMA / (TMP * S)
 40   CONTINUE

      X = ZDSQ
      YY = X - (1.d0, 0.0d0)
      TMP = ZDSQ - (0.5d0, 0.0d0)
      TMP = ( (TMP + G(SPD)) * INVEE )**TMP
      S = COEF(SPD, 0)
      DO 50 K = 1, NN(SPD)-1
      YY = YY + (1.d0, 0.0d0)
 50   S = S + COEF(SPD, K) / YY
      MULTIGAMMA = MULTIGAMMA / (TMP * S * TMP * S)

      RETURN
      END
C     ***

C     ****f* special.f/RATGAMMA
C  NAME
C     RATGAMMA  --  Ratio of two gamma functions
C  DESCRIPTION
C     Uses Lanczos algorithm
C  SYNOPSIS

      DOUBLE COMPLEX FUNCTION RATGAMMA(Z1,Z2)

      IMPLICIT NONE
      DOUBLE COMPLEX Z1, Z2

C  INPUTS
C     Z1, Z2 -- arguments
C  RETURN VALUE
C     RATGAMMA -- Gamma(Z1)/Gamma(Z2)
C  PARENTS
C     HJ
C  SOURCE
C

      INTEGER SPD
#include "header.F"
      INTEGER K
      INTEGER NMAX, NP
      PARAMETER (NMAX = 7, NP = 4)
      INTEGER NN(NP)
      DOUBLE PRECISION INVEE, G(NP), COEF(NP, 0:NMAX-1)
      PARAMETER (INVEE = 0.3678794411714423215955d0)
      DOUBLE COMPLEX TMP, X, YY, S
      SAVE NN, G, COEF
      DATA NN            /7, 3, 2, 1/
      DATA G            /5.0d0, 2.60404494991116d0, 
     &                   1.49614999d0, 0.328498d0/
      DATA COEF(4,0)     /1.8113464d0/
      DATA COEF(3,0), COEF(3,1) /0.5613831815d0,0.6055625385d0/
      DATA COEF(2,0), COEF(2,1), COEF(2,2) /
     &                    0.1854248319548753,0.8797922715613486,
     &                                      -0.2588333483439386/
      DATA COEF(1,0), COEF(1,1), COEF(1,2), COEF(1,3), COEF(1,4), 
     &     COEF(1,5), COEF(1,6)     /0.0168895284640819927d0, 
     &                   1.2866458274168037d0,
     &                  -1.46103406972059703d0, 
     &                   0.405586795700707467d0,
     &                  -0.0208035005652801098, 
     &                   0.0000204135450223743664,
     &                  -9.11230491453873551e-8/


      SPD = 3
      IF (SPEED .EQ. 1) SPD = 1

      X = Z1
      YY = X - (1.d0, 0.0d0)
      TMP = Z1 - (0.5d0, 0.0d0)
      TMP = ( (TMP + G(SPD)) * INVEE )**TMP
      S = COEF(SPD, 0)
      DO 60 K = 1, NN(SPD)-1
      YY = YY + (1.d0, 0.0d0)
 60   S = S + COEF(SPD, K) / YY
      RATGAMMA = TMP * S

      X = Z2
      YY = X - (1.d0, 0.0d0)
      TMP = Z2 - (0.5d0, 0.0d0)
      TMP = ( (TMP + G(SPD)) * INVEE )**TMP
      S = COEF(SPD, 0)
      DO 70 K = 1, NN(SPD)-1
      YY = YY + (1.d0, 0.0d0)
 70   S = S + COEF(SPD, K) / YY
      RATGAMMA = RATGAMMA / (TMP * S)

      RETURN
      END
C     ***


** C     ****f* special.f/CLNGAMMA
** C  NAME
** C     CLNGAMMA  --  logarithm of gamma function of complex argument
** C  DESCRIPTION
** C     Uses GSL library function gsl_sf_lngamma_complex_e via C wrapper clg_
** C  SYNOPSIS
** 
**       DOUBLE COMPLEX FUNCTION CLNGAMMA(XX)
** 
**       IMPLICIT NONE
**       DOUBLE COMPLEX XX
** 
** C  INPUTS
** C     XX -- argument
** C  RETURN VALUE
** C     CLNGAMMA -- ln(Gamma(z))
** C  PARENTS
** C     INIT, CBETA, HJ
** C  CHILDREN
** C     clg_
** C  SOURCE
** C
**       DOUBLE PRECISION REZ, IMZ, RELG, IMLG
**       DOUBLE COMPLEX Z
** 
**       REZ = REALPART(XX)
**       IMZ = IMAGPART(XX)
** 
**       CALL CLG(REZ, IMZ, RELG, IMLG)
** 
**       CLNGAMMA = DBLCMPLX(RELG, IMLG)
** 
**       RETURN
**       END
** C     ***


C     ****f* special.f/CBETA
C  NAME
C     CBETA  --  beta function of complex argument
C  SYNOPSIS

      DOUBLE COMPLEX FUNCTION CBETA(Z, W)

      DOUBLE COMPLEX Z, W 

C  INPUTS
C     Z, W -- arguments
C  RETURN VALUE
C     CBETA -- Beta(Z, W)
C  PARENTS
C     FCN, HJ
C  CHILDREN
C     CLNGAMMA  -- log(Gamma(z)) function
C  SOURCE
C

      DOUBLE COMPLEX CLNGAMMA

      CBETA = EXP( CLNGAMMA(Z) + CLNGAMMA(W) - CLNGAMMA(Z + W) )

      RETURN
      END
C     ***
      
C     ****f* special.f/POCHHAMMER
C  NAME
C     POCHHAMMER  --  Pochhammer symbol
C  SYNOPSIS

      DOUBLE COMPLEX FUNCTION POCHHAMMER(Z, M)

      INTEGER M
      DOUBLE COMPLEX Z

C  INPUTS
C     Z, M -- arguments
C  RETURN VALUE
C     POCHHAMMER -- (Z)_M
C  PARENTS
C     HJ
C  SOURCE
C

      INTEGER K

      POCHHAMMER = Z
      DO 100 K = 1, M - 1
 100     POCHHAMMER = POCHHAMMER * (Z + K)  

      RETURN
      END
C     ***
