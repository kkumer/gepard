C     ****h* gepard/projectors.f
C  FILE DESCRIPTION
C    projectors to eigenstates of the LO singlet 
C    anomalous dimensions matrix.
C
C    $Id$
C     *******


C     ****s* projectors.f/PROJECTORSF
C  NAME
C     PROJECTORSF  --   projectors P^+ and P^-
C  DESCRIPTION
C    projectors to eigenstates of the LO singlet 
C    anomalous dimensions matrix.
C  SYNOPSIS
C     SUBROUTINE PROJECTORSF (K, LAM, PR)
C
C     INTEGER K
C     DOUBLE COMPLEX LAM(2), PR(2,2,2)
C  INPUTS
C           K -- Mellin-Barnes integration point index
C         LAM -- eigenvalues of LO evolution operator
C  OUTPUT
C          PR -- Two projector matrices PR(1,a,b)=P^+
C                and PR(2,a,b)=P^-
C  IDENTIFIERS
C        NGAM -- common block with anomalous dimensions of DIS
C  PARENTS
C      RNNLOF
C  CHILDREN
C      LAMBDAF
C  SOURCE
C

      SUBROUTINE PROJECTORSF (SEC, K, LAM, PR)

      IMPLICIT NONE
      INTEGER SEC, K
      DOUBLE COMPLEX LAM(2), PR(2,2,2)
      DOUBLE COMPLEX DEN
      INCLUDE 'header.f'


      DEN = 1.0d0 / (LAM(1) - LAM(2)) 

*     P(+)

      PR(1,1,1) = DEN * (GAM(SEC,K,0,1,1) - LAM(2)) 
      PR(1,1,2) = DEN *  GAM(SEC,K,0,1,2)
      PR(1,2,1) = DEN *  GAM(SEC,K,0,2,1)
      PR(1,2,2) = DEN * (GAM(SEC,K,0,2,2) - LAM(2)) 

*     Using P(-) = 1 - P(+) 

      PR(2,1,1) = 1.0d0 - PR(1,1,1) 
      PR(2,1,2) = - PR(1,1,2) 
      PR(2,2,1) = - PR(1,2,1) 
      PR(2,2,2) = 1.0d0 - PR(1,2,2) 

      RETURN
      END
C     ****

C     ****s* projectors.f/PROJECTORSNDF
C  NAME
C     PROJECTORSNDF  --   projectors P^+ and P^-
C  DESCRIPTION
C    projectors to eigenstates of the LO singlet 
C    anomalous dimensions matrix. - version needed for non-diagonal
C    evolution
C  SYNOPSIS
C     SUBROUTINE PROJECTORSF (K, PR)
C
C     INTEGER K
C     DOUBLE COMPLEX PR(2,2,2)
C  INPUTS
C           K -- Mellin-Barnes integration point index
C  OUTPUT
C          PR -- Two projector matrices PR(1,a,b)=P^+
C                and PR(2,a,b)=P^-
C  IDENTIFIERS
C        NGAM -- common block with anomalous dimensions of DIS
C  PARENTS
C      RNNLOF
C  CHILDREN
C      LAMBDAF
C  SOURCE
C

      SUBROUTINE PROJECTORSNDF (GAMN, GAMK, LAMN, LAMK, PN, PK)

      IMPLICIT NONE
      INTEGER K
      DOUBLE COMPLEX GAMN(2,2), GAMK(2,2), LAMN(2), LAMK(2)
      DOUBLE COMPLEX PN(2,2,2), PK(2,2,2)
      DOUBLE COMPLEX DEN
      INCLUDE 'header.f'


      DEN = 1.0d0 / (LAMN(1) - LAMN(2)) 

*     P(+)

      PN(1,1,1) = DEN * (GAMN(1,1) - LAMN(2)) 
      PN(1,1,2) = DEN *  GAMN(1,2)
      PN(1,2,1) = DEN *  GAMN(2,1)
      PN(1,2,2) = DEN * (GAMN(2,2) - LAMN(2)) 

*     Using P(-) = 1 - P(+) 

      PN(2,1,1) = 1.0d0 - PN(1,1,1) 
      PN(2,1,2) = - PN(1,1,2) 
      PN(2,2,1) = - PN(1,2,1) 
      PN(2,2,2) = 1.0d0 - PN(1,2,2) 

      DEN = 1.0d0 / (LAMK(1) - LAMK(2)) 

*     P(+)

      PK(1,1,1) = DEN * (GAMK(1,1) - LAMK(2)) 
      PK(1,1,2) = DEN *  GAMK(1,2)
      PK(1,2,1) = DEN *  GAMK(2,1)
      PK(1,2,2) = DEN * (GAMK(2,2) - LAMK(2)) 

*     Using P(-) = 1 - P(+) 

      PK(2,1,1) = 1.0d0 - PK(1,1,1) 
      PK(2,1,2) = - PK(1,1,2) 
      PK(2,2,1) = - PK(1,2,1) 
      PK(2,2,2) = 1.0d0 - PK(1,2,2) 

      RETURN
      END
C     ****
