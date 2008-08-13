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
C  NOTES
C     Error for P=1 is strictly smaller than 2x10^(-10), and in practice is often
C     less then 10^(-14). This is original Lanczos parameterization.
C  PARENTS
C     INIT, CBETA, HJ
C  BUGS
C     It should be programmed to calculate lngamma directly, and not by
C     taking ln at the end. For this, restriction to P.V. of ln shold be
C     implemented.
C  SOURCE
C

      INTEGER SPD
      INCLUDE 'header.f'
      INTEGER K
      INTEGER NMAX, NP
      PARAMETER (NMAX = 7, NP = 4)
      INTEGER NN(NP)
      DOUBLE PRECISION EE, G(NP), COEF(NP, 0:NMAX-1)
      DOUBLE COMPLEX TMP, X, YY, S
      SAVE NN, EE, G, COEF
      DATA NN            /7, 3, 2, 1/
      DATA EE           /2.71828182845904524d0/
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
      TMP = ( (TMP + G(SPD)) / EE )**TMP
      S = COEF(SPD, 0)
      DO 10 K = 1, NN(SPD)-1
      YY = YY + (1.d0, 0.0d0)
 10   S = S + COEF(SPD, K) / YY

      CLNGAMMA = LOG(TMP * S)

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
**       CLNGAMMA = COMPLEX(RELG, IMLG)
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
      DO 10 K = 1, M - 1
 10     POCHHAMMER = POCHHAMMER * (Z + K)  

      RETURN
      END
C     ***
