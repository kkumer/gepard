C     ****h* gepard/projectors.f
C  FILE DESCRIPTION
C    projectors to eigenstates of the LO singlet 
C    anomalous dimensions matrix.
C
C    $Id: projectors.f,v 1.2 2006-05-12 19:55:34+02 kuk05260 Exp kuk05260 $
C     *******


C     ****s* projectors.f/PROJECTORSF
C  NAME
C     PROJECTORSF  --   projectors P^+ and P^-
C  DESCRIPTION
C    projectors to eigenstates of the LO singlet 
C    anomalous dimensions matrix.
C  SYNOPSIS
C     SUBROUTINE PROJECTORSF (NF, PR)
C
C     INTEGER NF
C     DOUBLE COMPLEX PR(2,2,2)
C  INPUTS
C          NF -- number of active flavours
C  OUTPUT
C          PR -- Two projector matrices PR(1,a,b)=P^+
C                and PR(2,a,b)=P^-
C  IDENTIFIERS
C      WGAMMA -- common block with anomalous dimensions of DIS
C  PARENTS
C      RNNLOF
C  CHILDREN
C      LAMBDAF
C  SOURCE
C

      SUBROUTINE PROJECTORSF (NF, PR)

      IMPLICIT NONE
      INTEGER NF
      DOUBLE COMPLEX PR(2,2,2)
      DOUBLE COMPLEX LAM(2), GAM0(2,2), GAM1(2,2), GAM2(2,2)
      DOUBLE COMPLEX DEN

*     Input common-block

      COMMON / WGAMMA  /  GAM0, GAM1, GAM2

      CALL LAMBDAF (NF, LAM)

      DEN = 1.0d0 / (LAM(1) - LAM(2)) 

*     P(+)

      PR(1,1,1) = DEN * (GAM0(1,1) - LAM(2)) 
      PR(1,1,2) = DEN *  GAM0(1,2)
      PR(1,2,1) = DEN *  GAM0(2,1)
      PR(1,2,2) = DEN * (GAM0(2,2) - LAM(2)) 

*     Using P(-) = 1 - P(+) 

      PR(2,1,1) = 1.0d0 - PR(1,1,1) 
      PR(2,1,2) = - PR(1,1,2) 
      PR(2,2,1) = - PR(1,2,1) 
      PR(2,2,2) = 1.0d0 - PR(1,2,2) 

      RETURN
      END
C     ****
