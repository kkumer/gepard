C     ****h* gepard/special.f
C  FILE DESCRIPTION
C    Double complex special functions (Gamma and Beta)
C
C    $Id$
C  NOTES
C     - Adapted from "Numerical recipes" Sect. 6.1
C     *******


C     ****f* special.f/CLNGAMMA
C  NAME
C     CLNGAMMA  --  logarithm of gamma function of complex argument
C  DESCRIPTION
C     Uses Lanczos algorithm
C  SYNOPSIS
C     DOUBLE COMPLEX FUNCTION CLNGAMMA(Z)
C
C     DOUBLE COMPLEX Z
C  INPUTS
C     Z -- argument
C  RETURN VALUE
C     CLNGAMMA -- ln(Gamma(z))
C  NOTES
C     Error is strictly smaller than 2x10^(-10), and in practice is often
C     less then 10^(-14)
C  PARENTS
C     INIT, PARWAVF
C  SOURCE
C

      DOUBLE COMPLEX FUNCTION CLNGAMMA(XX)

      INTEGER J
      DOUBLE PRECISION STP, COF(6)
      DOUBLE COMPLEX XX, SER, TMP, X, Y
      SAVE COF, STP 
      DATA COF, STP /76.18009172947146d0,-86.50532032941677d0,
     &  24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     &    -.5395239384953d-5,2.5066282746310005d0/

      X=XX
      Y=X
      TMP=X+(5.5d0,0.0d0)
      TMP=(X+(0.5d0,0.0d0))*LOG(TMP)-TMP
      ser=(1.000000000190015d0, 0.0d0)
      DO 10 J=1, 6
      Y = Y + (1.d0, 0.0d0)
 10   SER = SER + COF(J) / Y
      CLNGAMMA = TMP + LOG( STP * SER / X )

      RETURN
      END
C     ***


C     ****f* special.f/CBETA
C  NAME
C     CBETA  --  beta function of complex argument
C  SYNOPSIS
C     DOUBLE COMPLEX FUNCTION CBETA(Z, W)
C
C     DOUBLE COMPLEX Z
C     DOUBLE COMPLEX W
C  INPUTS
C     Z, W -- arguments
C  RETURN VALUE
C     CBETA -- Beta(Z, W)
C  PARENTS
C     HJ
C  CHILDREN
C     CLNGAMMA  -- log(Gamma(z)) function
C  SOURCE
C

      DOUBLE COMPLEX FUNCTION CBETA(Z, W)

      DOUBLE COMPLEX Z, W, CLNGAMMA

      CBETA = EXP( CLNGAMMA(Z) + CLNGAMMA(W) - CLNGAMMA(Z + W) )

      RETURN
      END
C     ***
      
C     ****f* special.f/POCHHAMMER
C  NAME
C     POCHHAMMER  --  Pochhammer symbol
C  SYNOPSIS
C     DOUBLE COMPLEX FUNCTION POCHHAMMER(Z, M)
C
C     INTEGER M
C     DOUBLE COMPLEX Z
C  INPUTS
C     Z, M -- arguments
C  RETURN VALUE
C     POCHAMMER -- (Z)_M
C  PARENTS
C     HJ
C  SOURCE
C

      DOUBLE COMPLEX FUNCTION POCHHAMMER(Z, M)

      INTEGER M, K
      DOUBLE COMPLEX Z

      POCHHAMMER = Z
      DO 10 K = 1, M - 1
 10     POCHHAMMER = POCHHAMMER * (Z + K)  

      RETURN
      END
C     ***
