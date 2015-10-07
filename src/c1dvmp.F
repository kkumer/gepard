C     ****h* gepard/c1dvmp.F
C  FILE DESCRIPTION
C    calculation of NLO Wilson coefficients for DVEM 
C    according to DM's c1 notebook
C
C    $Id:$
C     *******


C     ****s* c1dvmp.F/C1DVMP
C  NAME
C     C1DVMP  --  NLO Wilson coefficients for DVEM
C  DESCRIPTION
C    calculates  Wilson coefficients for DVEM
C    normalized w.r.t. LO C=1.
C    OUTPUT MCU(L) L=1,2,3 for quark, pure singlet quark and gluon contrib
C  SYNOPSIS

      SUBROUTINE C1DVMP (SGNTR, J, K, MCU)

!f2py integer intent(in) :: SGNTR
!f2py complex(8) intent(in) :: J, K
!f2py complex(8) intent(out) :: MCU

      IMPLICIT NONE
      INTEGER SGNTR
      DOUBLE COMPLEX J, K, MCU(3)

C  SOURCE
C

      DOUBLE COMPLEX SUMA, SUMB, SUMANDB
      DOUBLE COMPLEX DELc1aGJK, DELc1bGKJ, LI
      DOUBLE COMPLEX MCQ1, MCQ1CF, MCQ1CG, MCQ1BET0
      DOUBLE COMPLEX MCG1, MCG1CF, MCG1CA, MCG1BET0
      DOUBLE COMPLEX MCPS1
      DOUBLE PRECISION LRGPDF2, LRDAF2, LRR2
      DOUBLE COMPLEX HS1, HS2, HS3, HSM1, HSM2, HSM3, HS1M2, SB3
      EXTERNAL HS1, HS2, HS3, HSM1, HSM2, HSM3, HS1M2, SB3
      INTEGER PTY, L, KINT
      EXTERNAL PTY
#include "header.F"


      CALL BETAF

      LRGPDF2 = LOG(RGPDF2)
      LRDAF2 = LOG(RDAF2)
      LRR2 = LOG(RR2)

      KINT = INT(K)

*  Spliced from Mathematica:

*  ... quark part

        MCQ1CF=-7.666666666666667+(0.5*(1.+3.*(1.+J)*(2.+J)))/((1
     &  .+J)**2*(2.+J)**2)+(0.5*(1.+3.*(1.+K)*(2.+K)))/((1.+K)**2
     &  *(2.+K)**2)+0.5*(-3.-2./((1.+J)*(2.+J))+4.*HS1(1.+J))*((-
     &  0.5*(1.+(1.+J)*(2.+J)))/((1.+J)*(2.+J))-(0.5*(1.+(1.+K)*(
     &  2.+K)))/((1.+K)*(2.+K))-LRGPDF2+HS1(1.+J)+HS1(1.+K))+0.5*
     &  ((-0.5*(1.+(1.+J)*(2.+J)))/((1.+J)*(2.+J))-(0.5*(1.+(1.+K
     &  )*(2.+K)))/((1.+K)*(2.+K))-LRDAF2+HS1(1.+J)+HS1(1.+K))*(-
     &  3.-2./((1.+K)*(2.+K))+4.*HS1(1.+K))

        MCQ1BET0=-0.8333333333333334+0.5/((1.+J)*(2.+J))+0.5/((1.
     &  +K)*(2.+K))+0.5*LRR2-HS1(1.+J)-HS1(1.+K)


      SUMA = (0.d0, 0.d0)
      SUMB = (0.d0, 0.d0)

      DO 10 L = 0, KINT

        LI = DBLCMPLX(L, 0)

        SUMANDB=(0.125*(1.+2.*LI)*(HS2(0.5*(1.+J))-HS2(-0.5+0.5*(
     &  1.+J))+HS2(-0.5+0.5*LI)-HS2(0.5*LI)))/((0.5*(1.+J)-0.5*LI
     &  )*(2.+J+LI))

        SUMB = SUMB + SUMANDB
        SUMA = SUMA + PTY(L)*SUMANDB

 10   CONTINUE

c     WRITE (*,*) 'SUMA = ', SUMA
c     WRITE (*,*) 'SUMB = ', SUMB

        DELc1aGJK=(-2.*PTY(KINT))/((1.+J)*(2.+J)*(1.+K)*(2.+K))+(0.5
     &  *(-1.+(1.+J)*(2.+J))*(-2.+(1.+K)*(2.+K)*(HS2(0.5*(1.+K))-
     &  HS2(-0.5+0.5*(1.+K))))*PTY(KINT))/((1.+J)*(2.+J))+(1.+K)*(2.
     &  +K)*(2.+0.5*(1.+K)*(2.+K))*((0.25*K*(2.+(1.+K)**2)*(HS2(0
     &  .5*(1.+J))-HS2(-0.5+0.5*(1.+J))+HS2(-0.5+0.5*K)-HS2(0.5*K
     &  )))/((0.5*(1.+J)-0.5*K)*(2.+J+K)*(3.+2.*K)*(4.+(1.+K)*(2.
     &  +K)))-(0.25*(HS2(0.5*(1.+J))-HS2(-0.5+0.5*(1.+J))-HS2(0.5
     &  *(1.+K))+HS2(-0.5+0.5*(1.+K))))/((0.5*(1.+J)+0.5*(-1.-K))
     &  *(3.+J+K))+(0.25*(3.+K)*(2.+(2.+K)**2)*(HS2(0.5*(1.+J))-H
     &  S2(-0.5+0.5*(1.+J))-HS2(0.5*(2.+K))+HS2(-0.5+0.5*(2.+K)))
     &  )/((0.5*(1.+J)+0.5*(-2.-K))*(4.+J+K)*(3.+2.*K)*(4.+(1.+K)
     &  *(2.+K))))*PTY(KINT)+2.*(J-K)*(3.+J+K)*(-SUMA-ZETA3+HS3(1.+J
     &  )+(0.125*(1.+K)*(HS2(0.5*(1.+J))-HS2(-0.5+0.5*(1.+J))-HS2
     &  (0.5*(1.+K))+HS2(-0.5+0.5*(1.+K)))*PTY(KINT))/((0.5*(1.+J)+0
     &  .5*(-1.-K))*(3.+J+K)))

        DELc1bGKJ=1/((1.+K)*(2.+K))+0.5*(-2.-(1.+K)**2-((1.+J)*(2
     &  .+J))/((1.+K)*(2.+K)))*(HS2(0.5*(1.+J))-HS2(-0.5+0.5*(1.+
     &  J)))-0.5*(1.+K)*(HS2(0.5*(1.+K))-HS2(-0.5+0.5*(1.+K)))-(0
     &  .125*(1.+K)*(2.+K)*(4.+(1.+K)*(2.+K))*(HS2(0.5*(1.+J))-HS
     &  2(-0.5+0.5*(1.+J))-HS2(0.5*(1.+K))+HS2(-0.5+0.5*(1.+K))))
     &  /((0.5*(1.+J)+0.5*(-1.-K))*(3.+J+K))-(0.5*(1.+K)*(2.+K)*(
     &  (0.25*K*(2.+(1.+K)**2)*(HS2(0.5*(1.+J))-HS2(-0.5+0.5*(1.+
     &  J))+HS2(-0.5+0.5*K)-HS2(0.5*K)))/((0.5*(1.+J)-0.5*K)*(2.+
     &  J+K))+(0.25*(3.+K)*(2.+(2.+K)**2)*(HS2(0.5*(1.+J))-HS2(-0
     &  .5+0.5*(1.+J))-HS2(0.5*(2.+K))+HS2(-0.5+0.5*(2.+K))))/((0
     &  .5*(1.+J)+0.5*(-2.-K))*(4.+J+K))))/(3.+2.*K)+2.*(-J+K)*(3
     &  .+J+K)*(-SUMB-0.5*HS1(1.+K)*(HS2(0.5*(1.+J))-HS2(-0.5+0.5
     &  *(1.+J)))+Sb3(ONE+J))

        MCQ1CG=0.9565348003631189+DELc1aGJK-(2.*(1.+(1.+J)*(2.+J)
     &  )*(1.-SGNTR))/((1.+J)**2*(2.+J)**2)-DELc1bGKJ*SGNTR+(-(1/
     &  ((1.+K)*(2.+K)))+2.*HS1(1.+K))*(1.-SGNTR+0.5*(1.+J)*(2.+J
     &  )*SGNTR*(-HS2(0.5*J)+HS2(0.5*(1.+J))))+(2.*SGNTR*PTY(KINT))/
     &  ((1.+J)*(2.+J)*(1.+K)*(2.+K))-(2.*(1.+(1.+K)*(2.+K))*(1.+
     &  PTY(KINT)))/((1.+K)**2*(2.+K)**2)+(-(1/((1.+J)*(2.+J)))+2.*H
     &  S1(1.+J))*(1.+PTY(KINT)-0.5*(1.+K)*(2.+K)*(-HS2(0.5*K)+HS2(0
     &  .5*(1.+K)))*PTY(KINT))+2.*(1.+J)*(2.+J)*((-0.5*(-1.+(1.+J)*(
     &  2.+J)))/((1.+J)**2*(2.+J)**2)+ZETA3-(0.5*SGNTR*(-HS2(0.5*
     &  J)+HS2(0.5*(1.+J))))/((1.+J)*(2.+J))-HS3(1.+J)-SGNTR*Sb3(
     &  ONE+J))+2.*(1.+K)*(2.+K)*((-0.5*(-1.+(1.+K)*(2.+K)))/((1.+
     &  K)**2*(2.+K)**2)+ZETA3-HS3(1.+K)+(0.5*(-HS2(0.5*K)+HS2(0.
     &  5*(1.+K)))*PTY(KINT))/((1.+K)*(2.+K))+PTY(KINT)*Sb3(ONE+K))

c     WRITE (*,*) 'DELc1aGJK = ', DELc1aGJK
c     WRITE (*,*) 'DELc1bGKJ = ', DELc1bGKJ

c     WRITE (*,*) 'CF*MCQ1CF = ', CF*MCQ1CF
c     WRITE (*,*) 'BET1*MCQ1BET0 = ',BETA0(NF) * MCQ1BET0
c     WRITE (*,*) 'CG*MCQ1CG = ', CG*MCQ1CG

      MCQ1 = CF * MCQ1CF  +  CG * MCQ1CG  +  BETA0(NF) * MCQ1BET0



C  ... pure singlet quark part starting at NLO

        MCPS1=(-2.*(0.5+1/((1.+j)*(2.+j))+1/((1.+k)*(2.+k))))/((1
     &  .+j)*(2.+j))-(2.*(2.+(1.+j)*(2.+j))*(-1.-LRGPDF2+2.*HS1(1
     &  .+j)+2.*HS1(1.+k)))/(j*(1.+j)*(2.+j)*(3.+j))+(0.5*k*(1.+k
     &  )*(2.+k)*(3.+k)*((0.25*(HS2(0.5*(1.+j))-HS2(-0.5+0.5*(1.+
     &  j))+HS2(-0.5+0.5*k)-HS2(0.5*k)))/((0.5*(1.+j)-0.5*k)*(2.+
     &  j+k))-(0.25*(HS2(0.5*(1.+j))-HS2(-0.5+0.5*(1.+j))-HS2(0.5
     &  *(2.+k))+HS2(-0.5+0.5*(2.+k))))/((0.5*(1.+j)+0.5*(-2.-k))
     &  *(4.+j+k))))/(3.+2.*k)


C   ... gluon part

        MCG1CF=(-0.5*(2.+(1.+k)*(2.+k)))/((1.+j)*(2.+j)*(1.+k)*(2
     &  .+k))+HS1(1.+j)/((1.+j)*(2.+j))-(-1.5-3./((1.+j)*(2.+j))+
     &  0.5/((1.+k)*(2.+k))+HS1(1.+j))/((1.+k)*(2.+k))-(0.5*(2.+(
     &  1.+j)*(2.+j))*(-1.5-1/((1.+j)*(2.+j))+2./((1.+k)*(2.+k))-
     &  LRGPDF2+3.*HS1(1.+j)))/((1.+j)*(2.+j))+0.5*(-0.75-1/((1.+
     &  j)*(2.+j))-0.5/((1.+k)*(2.+k))-LRDAF2+HS1(1.+j)+HS1(1.+k)
     &  )*(-3.-2./((1.+k)*(2.+k))+4.*HS1(1.+k))+0.125*(-39.+(2.+(
     &  1.+k)*(2.+k))*(-HS2(0.5*k)+HS2(0.5*(1.+k))))+0.25*(1.+k)*
     &  (2.+k)*(2.+(1.+k)*(2.+k))*((0.5*(-HS2(0.5*j)+HS2(0.5*(1.+
     &  j))))/((1.+k)*(2.+k))-(0.5*((0.25*(-1.+k)*k*(HS2(0.5*(1.+
     &  j))-HS2(-0.5+0.5*(1.+j))+HS2(-0.5+0.5*k)-HS2(0.5*k)))/((0
     &  .5*(1.+j)-0.5*k)*(2.+j+k))-(0.25*(3.+k)*(4.+k)*(HS2(0.5*(
     &  1.+j))-HS2(-0.5+0.5*(1.+j))-HS2(0.5*(2.+k))+HS2(-0.5+0.5*
     &  (2.+k))))/((0.5*(1.+j)+0.5*(-2.-k))*(4.+j+k))))/(3.+2.*k)
     &  )

        MCG1CA=1.572467033424113-(4.+10.*(1.+j)*(2.+j))/((1.+j)**
     &  2*(2.+j)**2)-(3.*(-6.+2.*HS1(1.+j)+HS1(1.+k)))/(j*(3.+j))
     &  +0.5*(4./((1.+j)*(2.+j))-12./(j*(3.+j))+4.*HS1(1.+j))*((0
     &  .5*(2.+(1.+j)*(2.+j)))/((1.+j)*(2.+j))-LRGPDF2+HS1(1.+j)+
     &  1.5*HS1(1.+k))+(-2.+(1.+k)*(2.+k)*HS1(1.+k))/((1.+j)*(2.+
     &  j)*(1.+k)*(2.+k))+0.5*(HS2(0.5*j)-HS2(0.5*(1.+j)))+0.125*
     &  (2.-(1.+k)*(2.+k)*(-HS2(0.5*k)+HS2(0.5*(1.+k))))+0.25*k*(
     &  1.+k)*(2.+k)*(3.+k)*((-0.5*(-HS2(0.5*j)+HS2(0.5*(1.+j))))
     &  /((1.+k)*(2.+k))-((0.25*(-1.+k)*(HS2(0.5*(1.+j))-HS2(-0.5
     &  +0.5*(1.+j))+HS2(-0.5+0.5*k)-HS2(0.5*k)))/((0.5*(1.+j)-0.
     &  5*k)*(2.+j+k))+(0.25*(4.+k)*(HS2(0.5*(1.+j))-HS2(-0.5+0.5
     &  *(1.+j))-HS2(0.5*(2.+k))+HS2(-0.5+0.5*(2.+k))))/((0.5*(1.
     &  +j)+0.5*(-2.-k))*(4.+j+k)))/(3.+2.*k))

        MCG1BET0=-0.5*LRGPDF2+0.5*LRR2

      MCG1 = CF * MCG1CF  +  CA * MCG1CA  +  BETA0(NF) * MCG1BET0


      MCU(1) = MCQ1
      MCU(2) = MCPS1
      MCU(3) = MCG1


      RETURN
      END
C     *******



      DOUBLE COMPLEX FUNCTION SB3(J)
       
      IMPLICIT NONE
      DOUBLE COMPLEX J

      DOUBLE COMPLEX HS1, HS2, HS3, MELLINF2
      EXTERNAL HS1, HS2, HS3, MELLINF2

        SB3=0.5*HS1(j)*(-HS2(-0.5+0.5*j)+HS2(0.5*j))+0.125*(-HS3(
     &  -0.5+0.5*j)+HS3(0.5*j))-2.*(0.8224670334241131*(-HS1(0.5*
     &  (-1.+j))+HS1(0.5*j))-MELLINF2(1.+J))

      RETURN
      END


      INTEGER FUNCTION PTY(K)
       
      IMPLICIT NONE
      INTEGER K

C     parity of K i.e. (-1)^K
      PTY = (1 - 2*MOD(K, 2))
      RETURN
      END
