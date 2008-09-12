C     ****h* gepard/splice.F
C  FILE DESCRIPTION
C    Conformal moment obtained by Mathematica splicing into
C    splice_template.F
C
C    $Id$
C     *******


C     ****s* splice.F/SPLICE
C  NAME
C     SPLICE  --  conformal moment of input-scale singlet GPD  H_{J}
C  DESCRIPTION
C     returns H_{J} for ansatz spliced from Mathematica
C  SYNOPSIS

      SUBROUTINE SPLICE(J, FCM)

      IMPLICIT NONE
      DOUBLE COMPLEX J, FCM(2) 

C  INPUTS
C           J -- conformal moment
C  OUTPUT
C         FCM -- input scale singlet GPD H_{J}
C  PARENTS
C     HJ
C  CHILDREN
C     CLNGAMMA
C  BUGS
C     This is very broken!!
C  SOURCE
C
#include "header.F"
C     DOUBLE COMPLEX CLNGAMMA

      IF ( PROCESS .EQ. 'DIS' ) THEN 
                FCM(1)=(PAR(11)*(2-PAR(12))*(3-PAR(12))*(4
     &  -PAR(12))*(5-PAR(12))*(6-PAR(12))*(7-PAR(1
     &  2))*(8-PAR(12))*(9-PAR(12))*(1/(1-DEL2/(PA
     &  R(14)+j*PAR(15))))**PAR(16))/((2+j-PAR(12)
     &  )*(3+j-PAR(12))*(4+j-PAR(12))*(5+j-PAR(12)
     &  )*(6+j-PAR(12))*(7+j-PAR(12))*(8+j-PAR(12)
     &  )*(1+j-PAR(12)-DEL2*PAR(13))) 
                FCM(2)=(PAR(21)*(2-PAR(22))*(3-PAR(22))*(4
     &  -PAR(22))*(5-PAR(22))*(6-PAR(22))*(7-PAR(2
     &  2))*(1/(1-DEL2/(PAR(24)+j*PAR(25))))**PAR(
     &  26))/((2+j-PAR(22))*(3+j-PAR(22))*(4+j-PAR
     &  (22))*(5+j-PAR(22))*(6+j-PAR(22))*(1+j-PAR
     &  (22)-DEL2*PAR(23))) 
      ELSE  
                FCM(1)=(PAR(11)*(2-PAR(12))*(3-PAR(12))*(4
     &  -PAR(12))*(5-PAR(12))*(6-PAR(12))*(7-PAR(1
     &  2))*(8-PAR(12))*(9-PAR(12))*(1/(1-DEL2/(PA
     &  R(14)+j*PAR(15))))**PAR(16))/((2+j-PAR(12)
     &  )*(3+j-PAR(12))*(4+j-PAR(12))*(5+j-PAR(12)
     &  )*(6+j-PAR(12))*(7+j-PAR(12))*(8+j-PAR(12)
     &  )*(1+j-PAR(12)-DEL2*PAR(13))) 
                FCM(2)=(PAR(21)*(2-PAR(22))*(3-PAR(22))*(4
     &  -PAR(22))*(5-PAR(22))*(6-PAR(22))*(7-PAR(2
     &  2))*(1/(1-DEL2/(PAR(24)+j*PAR(25))))**PAR(
     &  26))/((2+j-PAR(22))*(3+j-PAR(22))*(4+j-PAR
     &  (22))*(5+j-PAR(22))*(6+j-PAR(22))*(1+j-PAR
     &  (22)-DEL2*PAR(23))) 
      ENDIF

      RETURN
      END
C     ****
