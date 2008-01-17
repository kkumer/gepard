C     ****h* gepard/bca.f
C  FILE DESCRIPTION
C    calculation of DVCS beam charge asymmetry (BCA)
C    NB: Approximation of small XB is used!!!
C
C    $Id:$
C     *******


C     ****f* bca.f/BCA
C  NAME
C       BCA  --  DVCS beam charge asymmetry
C  SYNOPSIS

      SUBROUTINE BCA (BCA0)

      IMPLICIT NONE
      DOUBLE PRECISION BCA0

C  PARENTS
C     FCN
C  CHILDREN
C     TINT, TDVCS2, TBH2
C  SOURCE
C

      DOUBLE PRECISION SHERA, TAVG, Q2AVG, WAVG
      PARAMETER ( SHERA = 101568.0d0, TAVG = -0.228096d0, 
     &   Q2AVG = 8.0d0, WAVG = 82.0d0)
      DOUBLE PRECISION TINT, TDVCS2, TBH2
      INCLUDE 'header.f'

      PROCESS = 'DVCS'
      PHIAZ = 0.0d0

      Q2 = Q2AVG
*     NQS = 1
*     QS(1) = Q2
*     CALL EVOLC(1)
*  FIXME: this works only if Q2AVG is already in QS()

      DEL2 = TAVG
      W2 = WAVG**2
      XB = 2.0d0 * Q2AVG / (2.0d0*WAVG**2 + Q2AVG)
      YB = (W2 + Q2) / SHERA

      XI = XB / 2.0d0

      BCA0 = - TINT() / ( TDVCS2() + TBH2() )

      RETURN
      END
C     ***



C     ****f* bca.f/TBH2
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
      INCLUDE 'header.f'

      C0BH = 16.d0 * KKIN()**2 * (Q2/DEL2) *
     &  ( F1NUCL(DEL2)**2 - (DEL2/(4.0d0*MP2)) * F2NUCL(DEL2)**2 ) + 
     &        8.d0 * (2.0d0 - YB)**2 *
     &  ( F1NUCL(DEL2)**2 - (DEL2/(4.0d0*MP2)) * F2NUCL(DEL2)**2 )

      TBH2 =  C0BH / (XB**2 * YB**2 * DEL2 * P1() * P2()) 

      RETURN
      END
C     ***



C     ****f* bca.f/TDVCS2
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
      INCLUDE 'header.f'

      CALL CFFF

*     Multiplying with charge factor

      IF (NF .EQ. 3) THEN
         CFF(P) = CFF(P) * 2.0D0 / 9.0D0
         CFFE(P) = CFFE(P) * 2.0D0 / 9.0D0
      ELSE IF (NF .EQ. 4) THEN
         CFF(P) = CFF(P) * 5.0D0 / 18.0D0
         CFFE(P) = CFFE(P) * 5.0D0 / 18.0D0
      ELSE
         CALL ERROR ('GeParD', '  CFF',
     &   'NF is not integer equal to 3 or 4!                          ',
     &   5, 3)
      END IF

      CCALDVCS = ABS(CFF(P))**2 - (DEL2/(4.0d0*MP2)) * ABS(CFFE(P))**2
       
      C0DVCS = 2.0d0 * (2.0d0 - 2.0d0*YB + YB**2) * CCALDVCS

      TDVCS2 =  C0DVCS / ( YB**2 * Q2 )

      RETURN
      END
C     ***


C     ****f* bca.f/TINT
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
      INCLUDE 'header.f'

      CALL CFFF

*     Multiplying with charge factor

      IF (NF .EQ. 3) THEN
         CFF(P) = CFF(P) * 2.0D0 / 9.0D0
         CFFE(P) = CFFE(P) * 2.0D0 / 9.0D0
      ELSE IF (NF .EQ. 4) THEN
         CFF(P) = CFF(P) * 5.0D0 / 18.0D0
         CFFE(P) = CFFE(P) * 5.0D0 / 18.0D0
      ELSE
         CALL ERROR ('GeParD', '  CFF',
     &   'NF is not integer equal to 3 or 4!                          ',
     &   5, 3)
      END IF

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


C     ****f* bca.f/P1
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
      INCLUDE 'header.f'


      JJ = (1.0d0 - YB) * (1.0d0 + DEL2/Q2) - (2.0d0 - YB) * DEL2 / Q2 

      P1 = - ( JJ + 2.0d0 * KKIN() * COS(PHIAZ) ) / YB
 
      RETURN
      END
C     ***


C     ****f* bca.f/P2
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
      INCLUDE 'header.f'


      JJ = (1.0d0 - YB) * (1.0d0 + DEL2/Q2) - (2.0d0 - YB) * DEL2 / Q2 

      P2 = 1.0d0 + DEL2/Q2 + ( JJ + 2.0d0 * KKIN() * COS(PHIAZ) ) / YB
 
      RETURN
      END
C     ***


C     ****f* bca.f/KKIN
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
      INCLUDE 'header.f'

      KKIN = SQRT( -DEL2 * (1.d0 - YB) / Q2 )
 
      RETURN
      END
C     ***


C     ****f* bca.f/F1NUCL
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

C     ****f* bca.f/F2NUCL
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
