C     ****h* gepard/c1dvmp.F
C  FILE DESCRIPTION
C    calculation of NLO Wilson coefficients for DVEM 
C    according to DM's c1 notebook
C
C    $Id:$
C     *******


C     ****s* dvem.F/C1DVMP
C  NAME
C     C1DVMPF  --  NLO Wilson coefficients for DVEM
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

      DOUBLE COMPLEX SUMA, SUMB, SUMANDA, SUMANDB
      DOUBLE COMPLEX DELc1aGJK, DELc1bGKJ, LI
      DOUBLE COMPLEX MCQ1, MCQ1CF, MCQ1CG, MCQ1BET0
      DOUBLE COMPLEX MCG1, MCG1CF, MCG1CA, MCG1BET0
      DOUBLE COMPLEX MCPS1
      DOUBLE PRECISION LRGPDF2, LRDAF2, LRR2
      DOUBLE COMPLEX HS1, HS2, HS3, HSM1, HSM2, HSM3, HS1M2, SB3
      EXTERNAL HS1, HS2, HS3, HSM1, HSM2, HSM3, HS1M2, SB3
      INTEGER PTY, L
      EXTERNAL PTY
#include "header.F"


      CALL BETAF

      LRGPDF2 = LOG(RGPDF2)
      LRDAF2 = LOG(RDAF2)
      LRR2 = LOG(RR2)


*  Spliced from Mathematica:

*  ... quark part

<* MCQ1CFexpr *>

<* MCQ1BET0expr *>


      SUMA = (0.d0, 0.d0)
      SUMB = (0.d0, 0.d0)

      DO 10 L = 0, INT(K)

        LI = DBLCMPLX(L, 0)

<* SUMANDBexpr *>

        SUMB = SUMB + SUMANDB
        SUMA = SUMA + PTY(L)*SUMANDB

 10   CONTINUE

c     WRITE (*,*) 'SUMA = ', SUMA
c     WRITE (*,*) 'SUMB = ', SUMB

<* DELc1aGJKexpr *>

<* DELc1bGKJexpr *>

<* MCQ1CGexpr *>

c     WRITE (*,*) 'DELc1aGJK = ', DELc1aGJK
c     WRITE (*,*) 'DELc1bGKJ = ', DELc1bGKJ

c     WRITE (*,*) 'CF*MCQ1CF = ', CF*MCQ1CF
c     WRITE (*,*) 'BET1*MCQ1BET0 = ',BETA0(NF) * MCQ1BET0
c     WRITE (*,*) 'CG*MCQ1CG = ', CG*MCQ1CG

      MCQ1 = CF * MCQ1CF  +  CG * MCQ1CG  +  BETA0(NF) * MCQ1BET0



C  ... pure singlet quark part starting at NLO

<* MCPS1expr *>


C   ... gluon part

<* MCG1CFexpr *>

<* MCG1CAexpr *>

<* MCG1BET0expr *>

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

<* SB3expr *>

      RETURN
      END

