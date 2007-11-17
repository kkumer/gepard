C     ****h* gepard/astrong.f
C  FILE DESCRIPTION
C    running of QCD coupling constant
C
C    $Id$
C     *******

C     ****s* astrong.f/AS2PF
C  NAME
C     AS2PF  --  alpha strong divided by 2 Pi
C  DESCRIPTION
C    QCD coupling alpha strong is obtained by integrating the evolution
C    equation for a fixed number of massless flavours  NF.  Except at
C    leading order (LO), result is obtained using a fourth-order
C    Runge-Kutta integration. 
C  SYNOPSIS

      SUBROUTINE AS2PF (AS, R2, AS0, R20)
      
      IMPLICIT NONE
      DOUBLE PRECISION AS, R2, AS0, R20

C  INPUTS
C          R2 -- final momentum scale squared
C         AS0 -- initial value for a_strong/(2 Pi)
C         R20 -- initial momentum scale squared
C  OUTPUT
C          AS -- final value for a_strong/(2 Pi)
C  PARENTS
C     PARWAVF
C  NOTES
C     Requires a previous call of BETAF to initialize common block !
C     This is slightly modified routine from file asrgkt.f from package
C     'QCD-Pegasus' by Andreas Vogt, hep-ph/0408244
C  SOURCE
C

      DOUBLE PRECISION LRRAT, DLR, SXTH, A
      INTEGER NASTPS, K1
      DOUBLE PRECISION FBETA1, FBETA2, FBETA3
      DOUBLE PRECISION XK0, XK1, XK2, XK3
      PARAMETER (NASTPS = 20)
      PARAMETER ( SXTH = 0.16666 66666 66666 D0 )
      INCLUDE 'header.f'

*   The beta functions FBETAn at N^nLO for n = 1, 2, and 3

      FBETA1(A) = A**2 * ( BETA0(NF) + A *   BETA1(NF) )
      FBETA2(A) = A**2 * ( BETA0(NF) + A * ( BETA1(NF)
     &                       + A * BETA2(NF) ) )
      FBETA3(A) = A**2 * ( BETA0(NF) + A * ( BETA1(NF)
     &                        + A * (BETA2(NF) + A * BETA3(NF)) ) )

*   Initial value, evolution distance and step size

*   NB: AS below is from 4pi expansion and is returned to 2pi 
*          just before RETURN

      AS = 0.5D0 * AS0
      LRRAT = LOG (R2/R20)
      DLR = LRRAT / NASTPS

*   Solution of the evolution equation depending on  P
*   (fourth-order Runge-Kutta beyond the leading order)

      IF (P .EQ. 0) THEN

         AS = 0.5D0 * AS0 / (1.0D0 - 0.5D0 * BETA0(NF) * AS0 * LRRAT)

      ELSE IF (P .EQ. 1) THEN

      DO 2 K1 = 1, NASTPS
         XK0 = DLR * FBETA1 (AS)
         XK1 = DLR * FBETA1 (AS + 0.5D0 * XK0)
         XK2 = DLR * FBETA1 (AS + 0.5D0 * XK1)
         XK3 = DLR * FBETA1 (AS + XK2)
         AS = AS + SXTH * (XK0 + 2.0d0 * XK1 + 2.0D0 * XK2 + XK3)
  2   CONTINUE

      ELSE IF (P .EQ. 2) THEN

      DO 3 K1 = 1, NASTPS
         XK0 = DLR * FBETA2 (AS)
         XK1 = DLR * FBETA2 (AS + 0.5D0 * XK0)
         XK2 = DLR * FBETA2 (AS + 0.5D0 * XK1)
         XK3 = DLR * FBETA2 (AS + XK2)
         AS = AS + SXTH * (XK0 + 2.0D0 * XK1 + 2.0D0 * XK2 + XK3)
  3   CONTINUE

      ELSE IF (P .EQ. 3) THEN
*
      DO 4 K1 = 1, NASTPS
         XK0 = DLR * FBETA3 (AS)
         XK1 = DLR * FBETA3 (AS + 0.5D0 * XK0)
         XK2 = DLR * FBETA3 (AS + 0.5D0 * XK1)
         XK3 = DLR * FBETA3 (AS + XK2)
         AS = AS + SXTH * (XK0 + 2.0D0 * XK1 + 2.0D0 * XK2 + XK3)
  4   CONTINUE
      END IF

*  Return to .../(2pi)  expansion

      AS = 2.0D0 * AS

      RETURN
      END
C     ****
