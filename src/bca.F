C     ****h* gepard/bca.F
C  FILE DESCRIPTION
C    calculation of DVCS beam charge asymmetry (BCA)
C    NB: Approximation of small XB is used!!!
C
C    $Id:$
C     *******


C     ****f* bca.F/BCA
C  NAME
C       HERABCA  --  DVCS beam charge asymmetry for HERA s
C  SYNOPSIS

      SUBROUTINE BCA (FITFLAG, WAVG, Q2AVG, PHIIN, HERABCA)

      IMPLICIT NONE
      INTEGER FITFLAG
      DOUBLE PRECISION WAVG, Q2AVG, PHIIN, HERABCA

C  PARENTS
C     FCN
C  CHILDREN
C     INTEGRAF, TINT, TDVCS2, TBH2
C  SOURCE
C

      INTEGER K, NBNDMAX, NMINTS, EFFACC
      PARAMETER ( NBNDMAX = 1 )
      DOUBLE PRECISION SHERA
      PARAMETER ( SHERA = 101568.0d0 )
      DOUBLE PRECISION TINT, TDVCS2, TBH2
      DOUBLE PRECISION NUMER, DENOM
      DOUBLE PRECISION BND(NBNDMAX+1), MINTS(64), MINTWG(64)
#include "header.F"

*     Integration over t similar to how it is done for total
*     cross section. H1 have cut-off on  t=-0.05 GeV^2
      DATA BND /0.05D0, 1.0D0/ 
      EFFACC = ACC - SPEED + 1
      NMINTS = 2**EFFACC
      CALL INTEGRAF(EFFACC, 1, BND, NBNDMAX, MINTS, MINTWG)

      PROCESS = 'DVCS'
      PHIAZ = PHIIN

      Q2 = Q2AVG
      IF ( FITFLAG .NE. 1) THEN
*     We are not fitting, should set QS()
        NQS = 1
        QS(1) = Q2
        CALL EVOLC(1,1)
        CALL EVOLC(2,1)
      END IF

      W2 = WAVG**2
      XB = 2.0d0 * Q2AVG / (2.0d0*WAVG**2 + Q2AVG)
      YB = (W2 + Q2) / SHERA
      XI = XB / 2.0d0

*     Integrating over t
      NUMER = 0.0d0
      DENOM = 0.0d0
      DO 10 K = 1, NMINTS
          DEL2 = - MINTS(K)
          NUMER = NUMER + MINTWG(K) * TINT()
          DENOM = DENOM + MINTWG(K) * ( TDVCS2() + TBH2() )
 10   CONTINUE

      HERABCA = - NUMER / DENOM

      RETURN
      END
C     ***



C     ****f* bca.F/TBH2
C  NAME
C       TBH2  --  Bethe-Heitler amplitude squared
C  DESCRIPTION
C        Eq. (25) from BMK [hep-ph/0112108]
C  SYNOPSIS

      DOUBLE PRECISION FUNCTION TBH2 ()

C  PARENTS
C     BCA
C  CHILDREN
C     F1NUCL, F2NUCL, P1, P2, KKIN
C  SOURCE
C

      IMPLICIT NONE
      DOUBLE PRECISION C0BH
      DOUBLE PRECISION F1NUCL, F2NUCL, KKIN, P1, P2
#include "header.F"

      C0BH = 16.d0 * KKIN()**2 * (Q2/DEL2) *
     &  ( F1NUCL(DEL2)**2 - (DEL2/(4.0d0*MP2)) * F2NUCL(DEL2)**2 ) + 
     &        8.d0 * (2.0d0 - YB)**2 *
     &  ( F1NUCL(DEL2)**2 - (DEL2/(4.0d0*MP2)) * F2NUCL(DEL2)**2 )

      TBH2 =  C0BH / (XB**2 * YB**2 * DEL2 * P1() * P2()) 

      RETURN
      END
C     ***



C     ****f* bca.F/TDVCS2
C  NAME
C       TDVCS2  --  DVCS amplitude squared
C  DESCRIPTION
C        Eqs. (26), (43) from BMK [hep-ph/0112108] in small-XB approx.
C  SYNOPSIS

      DOUBLE PRECISION FUNCTION TDVCS2 ()

C  PARENTS
C     BCA
C  CHILDREN
C     .
C  SOURCE
C

      IMPLICIT NONE
      DOUBLE PRECISION CCALDVCS, C0DVCS 
#include "header.F"

      CALL CFFF

      CCALDVCS = ABS(CFF(P))**2 - (DEL2/(4.0d0*MP2)) * ABS(CFFE(P))**2
       
      C0DVCS = 2.0d0 * (2.0d0 - 2.0d0*YB + YB**2) * CCALDVCS

      TDVCS2 =  C0DVCS / ( YB**2 * Q2 )

      RETURN
      END
C     ***


C     ****f* bca.F/TINT
C  NAME
C       TINT  --  DVCS-BH interference amplitude
C  DESCRIPTION
C        Eqs. FIXME  from BMK [hep-ph/0112108] in small-XB approx.
C  SYNOPSIS

      DOUBLE PRECISION FUNCTION TINT ()

C  PARENTS
C     BCA
C  CHILDREN
C     CFFF, F1NUCL, F2NUCL, P1, P2, KKIN
C  SOURCE
C

      IMPLICIT NONE
      DOUBLE COMPLEX  CCALINT
      DOUBLE PRECISION C0INT, C1INT
      DOUBLE PRECISION F1NUCL, F2NUCL, KKIN, P1, P2
#include "header.F"

      CALL CFFF

      CCALINT = F1NUCL(DEL2) * CFF(P) - 
     &             (DEL2/(4.0d0*MP2)) * F2NUCL(DEL2) * CFFE(P)

      C0INT = -8.d0 * (2.d0 - YB) * (2.d0 - 2.d0 * YB + YB**2) *
     &  (-DEL2/Q2) * REALPART(CCALINT)

      C1INT = -8.d0 * (2.d0 - 2.d0 * YB + YB**2) *
     &      KKIN() * REALPART(CCALINT)

      TINT =  ( C0INT + C1INT * COS(PHIAZ) ) /
     &               ( XB * YB**3 * DEL2 * P1() *P2() )
       

      RETURN
      END
C     ***


C     ****f* bca.F/P1
C  NAME
C       P1  --  Bethe-Heitler propagator factor
C  DESCRIPTION
C        Eq. 32 from BMK [hep-ph/0112108]
C  SYNOPSIS

      DOUBLE PRECISION FUNCTION P1 ()

C  PARENTS
C     TBH2, TINT
C  CHILDREN
C     KKIN
C  SOURCE
C

      IMPLICIT NONE
      DOUBLE PRECISION JJ
      DOUBLE PRECISION KKIN
#include "header.F"


      JJ = (1.0d0 - YB) * (1.0d0 + DEL2/Q2) - (2.0d0 - YB) * DEL2 / Q2 

      P1 = - ( JJ + 2.0d0 * KKIN() * COS(PHIAZ) ) / YB
 
      RETURN
      END
C     ***


C     ****f* bca.F/P2
C  NAME
C       P2  --  Bethe-Heitler propagator factor
C  DESCRIPTION
C        Eq. 32 from BMK [hep-ph/0112108]
C  SYNOPSIS

      DOUBLE PRECISION FUNCTION P2 ()

C  PARENTS
C     TBH2, TINT
C  CHILDREN
C     KKIN
C  SOURCE
C

      IMPLICIT NONE
      DOUBLE PRECISION JJ
      DOUBLE PRECISION KKIN
#include "header.F"


      JJ = (1.0d0 - YB) * (1.0d0 + DEL2/Q2) - (2.0d0 - YB) * DEL2 / Q2 

      P2 = 1.0d0 + DEL2/Q2 + ( JJ + 2.0d0 * KKIN() * COS(PHIAZ) ) / YB
 
      RETURN
      END
C     ***


C     ****f* bca.F/KKIN
C  NAME
C       KKIN  --  Kinematical function K
C  DESCRIPTION
C        Eq. 30 from BMK [hep-ph/0112108]
C  SYNOPSIS

      DOUBLE PRECISION FUNCTION KKIN ()

C  PARENTS
C     TBH2, TINT, P1, P2
C  SOURCE
C

      IMPLICIT NONE
#include "header.F"

      KKIN = SQRT( -DEL2 * (1.d0 - YB) / Q2 )
 
      RETURN
      END
C     ***


C     ****f* bca.F/F1NUCL
C  NAME
C       F1NUCL  --  Nucleon elastic form factor F1
C  DESCRIPTION
C        Given to me by DM
C  SYNOPSIS

      DOUBLE PRECISION FUNCTION F1NUCL (T)

      IMPLICIT NONE
      DOUBLE PRECISION T

C  PARENTS
C     TBH2, TINT
C  SOURCE
C

       
      F1NUCL = (1.41d0 * (1.26d0 - T))/((0.71d0 - T)**2 * (3.53d0 - T))
 
      RETURN
      END
C     ***

C     ****f* bca.F/F2NUCL
C  NAME
C       F2NUCL  --  Nucleon elastic form factor F2
C  DESCRIPTION
C        Given to me by DM
C  SYNOPSIS

      DOUBLE PRECISION FUNCTION F2NUCL (T)

      IMPLICIT NONE
      DOUBLE PRECISION T

C  PARENTS
C     TBH2, TINT
C  SOURCE
C

      F2NUCL = 3.2d0 / ((0.71d0 - T)**2 * (3.53d0 - T))
 
      RETURN
      END
C     ***
