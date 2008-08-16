C     ****h* gepard/splice.f
C  FILE DESCRIPTION
C    Conformal moment obtained by Mathematica splicing into
C    splice_template.f
C
C    $Id:$
C     *******


C     ****s* splice.f/SPLICE
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
C     DCTAN, CLNGAMMA
C  SOURCE
C
#include "header.F"
      DOUBLE COMPLEX DCTAN, CLNGAMMA

      IF ( PROCESS .EQ. 'DIS' ) THEN 
        <* PDFQsplice *> 
        <* PDFGsplice *> 
      ELSE  
        <* GPDQsplice *> 
        <* GPDGsplice *> 
      ENDIF

      RETURN
      END
C     ****
