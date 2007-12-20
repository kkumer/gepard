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
      INCLUDE 'header.f'
      DOUBLE COMPLEX DCTAN, CLNGAMMA

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
     &  R(14)+j*PAR(15))))**PAR(16)*(1-(2**(-1-j+P
     &  AR(12)+DEL2*PAR(13)-2*(PAR(12)+DEL2*PAR(13
     &  )))*Pi**1.5d0*xi**(1+j-PAR(12)-DEL2*PAR(13
     &  ))*EXP(CLNGAMMA(1+j))*EXP(CLNGAMMA(1+j+PAR
     &  (12)+DEL2*PAR(13)))*PAR(19)*(1-(-sin(j*Pi)
     &  +sin(Pi*(PAR(12)+DEL2*PAR(13))))/sin(Pi*(j
     &  -PAR(12)-DEL2*PAR(13)))))/(EXP(CLNGAMMA(1.
     &  5d0+j))*EXP(CLNGAMMA(1+j-PAR(12)-DEL2*PAR(
     &  13)))*EXP(CLNGAMMA(PAR(12)+DEL2*PAR(13)))*
     &  *2*sin(Pi*(PAR(12)+DEL2*PAR(13))))))/((2+j
     &  -PAR(12))*(3+j-PAR(12))*(4+j-PAR(12))*(5+j
     &  -PAR(12))*(6+j-PAR(12))*(7+j-PAR(12))*(8+j
     &  -PAR(12))*(1+j-PAR(12)-DEL2*PAR(13))) 
                FCM(2)=(PAR(21)*(2-PAR(22))*(3-PAR(22))*(4
     &  -PAR(22))*(5-PAR(22))*(6-PAR(22))*(7-PAR(2
     &  2))*(1/(1-DEL2/(PAR(24)+j*PAR(25))))**PAR(
     &  26)*(1-(2**(-1-j+PAR(22)+DEL2*PAR(23)-2*(P
     &  AR(22)+DEL2*PAR(23)))*Pi**1.5d0*xi**(1+j-P
     &  AR(22)-DEL2*PAR(23))*EXP(CLNGAMMA(1+j))*EX
     &  P(CLNGAMMA(1+j+PAR(22)+DEL2*PAR(23)))*PAR(
     &  29)*(1-(-sin(j*Pi)+sin(Pi*(PAR(22)+DEL2*PA
     &  R(23))))/sin(Pi*(j-PAR(22)-DEL2*PAR(23))))
     &  )/(EXP(CLNGAMMA(1.5d0+j))*EXP(CLNGAMMA(1+j
     &  -PAR(22)-DEL2*PAR(23)))*EXP(CLNGAMMA(PAR(2
     &  2)+DEL2*PAR(23)))**2*sin(Pi*(PAR(22)+DEL2*
     &  PAR(23))))))/((2+j-PAR(22))*(3+j-PAR(22))*
     &  (4+j-PAR(22))*(5+j-PAR(22))*(6+j-PAR(22))*
     &  (1+j-PAR(22)-DEL2*PAR(23))) 
      ENDIF

      RETURN
      END
C     ****
