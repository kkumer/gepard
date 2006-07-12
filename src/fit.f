C     ****h* gepard/fit.f
C  FILE DESCRIPTION
C    Fitting of DVCS (and DIS in future) experimental data.
C
C    $Id$
C     *******

C     ****p* fit.f/FIT
C  NAME
C    FIT   --  Determines parameters of GPDs (and low-energy
C              input point Q0), by fitting to DVCS and DIS
C              experimental data.
C            
C  DESCRIPTION
C             Just a driver. Calls minuit subroutine
C             which does everything else.
C
C  INPUTS
C          FIT.CMD -- file with specification of fit parameters
C                     and with minuit commands.
C
C  OUTPUT
C          Goes to standard output. Good way to both observe
C          progress and save results is to use:
C            time ./fit | tee -a savefile
C
C  CHILDREN
C      MINUIT, FCN (via MINUIT)
C
C  SOURCE
C


      PROGRAM FIT

      CHARACTER SCHEME*5, ANSATZ*6
      EXTERNAL FCN

      COMMON / LABELS   /  SCHEME, ANSATZ

      OPEN (UNIT=5,FILE='FIT.CMD',STATUS='OLD')

      SCHEME = 'CSBAR'

      CALL READPAR
      CALL INIT

      CALL MINUIT(FCN,0)  

      STOP
      END
C     ****


C     ****s* fit.f/FCN
C  NAME
C     FCN  --  minimization subroutine called by minuit
C            
C  DESCRIPTION
C    Reads experimental datasets and calculates chi-square by
C    comparing to theoretical prediction functions.
C    Writes out theoretical line points for fit plotting.
C  SYNOPSIS
C     SUBROUTINE FCN(NPAR,GIN,F,A,IFLAG,FUTIL)
C
C     IMPLICIT NONE
C     INTEGER NPAR, IFLAG
C     DOUBLE PRECISION F, FUTIL
C     DOUBLE PRECISION A(NPAR),GIN(NPAR)
C    
C  INPUTS 
C            DATASET?  --  ? = 1, ..., 9, 0.  (name of 8 chars exactly!!)
C                          files with experimental data formatted
C                          like:   x  y  dy_stat  dy_syst
C  OUTPUT
C             FIT.PLT  --  Files with fit results formatted like:
C                          x  y(x)_theor
C                          (with empty lines separating graphs)
C  IDENTIFIERS
C                   N  --  Number of points in a given data set
C                   X  --  Array with x-values of experimental points
C                   Y  --  Array with y-values of experimental points
C                  DY  --  Array with errors of y-values
C                   A  --  Array with fit parameters
C              FITPAR  --  Equal to A
C
C  CHILDREN
C      PARSIGMA, SIGMA, F2
C
C  PARENTS
C      MINUIT, FCN (via MINUIT)
C
C  BUGS
C    - dimension of FITPAR is hardcoded 
C    - processing of datasets should be done
C      by subroutine
C
C  SOURCE
C


      SUBROUTINE FCN(NPAR,GIN,F,A,IFLAG,FUTIL)

      IMPLICIT NONE
      INTEGER NPAR, IFLAG
      DOUBLE PRECISION F, FUTIL
      DOUBLE PRECISION A(NPAR),GIN(NPAR)
      LOGICAL INCLUDEDVCS, INCLUDEDIS
      INTEGER N, NMAX, I
      PARAMETER (NMAX=40)
      DOUBLE PRECISION FITPAR(10)
      DOUBLE PRECISION X(NMAX), Y(NMAX), DY(NMAX)
      DOUBLE PRECISION CHISQ, FITFN
      DOUBLE PRECISION W2, XI, DEL2, Q2, Q02
      DOUBLE PRECISION PARSIGMA, SIGMA
      DOUBLE PRECISION F2(0:2)

      COMMON / PARLOG /  INCLUDEDVCS, INCLUDEDIS

      COMMON / KINEMATICS /  XI, DEL2, Q2, Q02
      COMMON / FITPARAMS  /  FITPAR
      COMMON / F2P        /  F2

*     Clear output buffer
      CALL FLUSHOUT()
      
*     Transfer parameters to common block
      DO 10 I = 1, NPAR
 10     FITPAR(I) = A(I)


      IF (IFLAG .EQ. 3)  THEN
*       file for writing out points of final fit
        OPEN (UNIT = 16, FILE = 'FIT.PLT', STATUS = 'UNKNOWN')
      END IF


      CHISQ = 0.0d0


      IF ( INCLUDEDVCS ) THEN

*     DATASET 1  [H1,  Eur.Phys.J.C44:1-11,2005, hep-ex/0505061]
*       X = -T,   Y = DSIG/DT   
*       W = 71 GeV, Q2 = 4 GeV

      W2 = 71.0d0**2
      Q2 = 4.0d0
      XI = Q2 / (2.0d0 * W2 + Q2)
      N = 4
      CALL READDATA ('DATASET1', N, X, Y, DY)
      DO I= 1, N
      FITFN = PARSIGMA ( X(I) )
      CHISQ = CHISQ + ( (FITFN - Y(I))**2 / DY(I)**2 )
      IF (IFLAG .EQ. 3)  THEN
        WRITE (16, *) X(I), FITFN
      END IF
      END DO

      IF (IFLAG .EQ. 3)  THEN
*       empty line for separation of different fit lines
        WRITE (16, *) 
      END IF

*     DATASET 2  [H1,  Eur.Phys.J.C44:1-11,2005, hep-ex/0505061]
*       X = -T,   Y = DSIG/DT   
*       W = 82 GeV, Q2 = 8 GeV

      W2 = 82.0d0**2
      Q2 = 8.2d0
      XI = Q2 / ( 2.0d0 * W2 + Q2)
      N = 4
      CALL READDATA ('DATASET2', N, X, Y, DY)
      DO I= 1, N
      FITFN = PARSIGMA ( X(I) )
      CHISQ = CHISQ + ( (FITFN - Y(I))**2 / DY(I)**2 )
      IF (IFLAG .EQ. 3)  THEN
        WRITE (16, *) X(I), FITFN
      END IF
      END DO

      IF (IFLAG .EQ. 3)  THEN
*       empty line for separation of different fit lines
        WRITE (16, *) 
      END IF

*     DATASET 3  [H1,  Eur.Phys.J.C44:1-11,2005, hep-ex/0505061]
*       X = Q2,   Y = SIG   (TCUT = - 1 GeV)
*       W = 82 GeV

      W2 = 82.0d0**2
      N = 6
      CALL READDATA ('DATASET3', N, X, Y, DY)
      DO I= 1, N
      Q2 = X(I)
      XI = Q2 / ( 2.0d0 * W2 + Q2)
      FITFN = SIGMA()
      CHISQ = CHISQ + ( (FITFN - Y(I))**2 / DY(I)**2 )
      IF (IFLAG .EQ. 3)  THEN
        WRITE (16, *) X(I), FITFN
      END IF
      END DO

      IF (IFLAG .EQ. 3)  THEN
*       empty line for separation of different fit lines
        WRITE (16, *) 
      END IF


*     DATASET 4  [ZEUS, Phys.Lett.B573:46-62,2003, hep-ex/0305028]
*       X = Q2,   Y = SIG   (TCUT = - \inf GeV ???)
*       W = 89 GeV

      W2 = 89.0d0**2
      N = 6
      CALL READDATA ('DATASET4', N, X, Y, DY)
      DO I= 1, N
      Q2 = X(I)
      XI = Q2 / ( 2.0d0 * W2 + Q2)
      FITFN = SIGMA()
      CHISQ = CHISQ + ( (FITFN - Y(I))**2 / DY(I)**2 )
      IF (IFLAG .EQ. 3)  THEN
        WRITE (16, *) X(I), FITFN
      END IF
      END DO

      IF (IFLAG .EQ. 3)  THEN
*       empty line for separation of different fit lines
        WRITE (16, *) 
      END IF


*     DATASET 5  [H1,  Eur.Phys.J.C44:1-11,2005, hep-ex/0505061]
*       X = W,   Y = SIG   (TCUT = - 1 GeV)
*       Q2 = 8 GeV^2

      Q2 = 8.0d0
      N = 5
      CALL READDATA ('DATASET5', N, X, Y, DY)
      DO I= 1, N
      W2 = X(I)**2
      XI = Q2 / ( 2.0d0 * W2 + Q2)
      FITFN = SIGMA()
      CHISQ = CHISQ + ( (FITFN - Y(I))**2 / DY(I)**2 )
      IF (IFLAG .EQ. 3)  THEN
        WRITE (16, *) X(I), FITFN
      END IF
      END DO

      IF (IFLAG .EQ. 3)  THEN
*       empty line for separation of different fit lines
        WRITE (16, *) 
      END IF


*     DATASET 6  [ZEUS, Phys.Lett.B573:46-62,2003, hep-ex/0305028]
*       X = W,   Y = SIG   (TCUT = - \inf GeV ???)
*       Q2 = 9.6 GeV

      Q2 = 9.6d0
      N = 10
      CALL READDATA ('DATASET6', N, X, Y, DY)
      DO I= 1, N
      W2 = X(I)**2
      XI = Q2 / ( 2.0d0 * W2 + Q2)
      FITFN = SIGMA()
      CHISQ = CHISQ + ( (FITFN - Y(I))**2 / DY(I)**2 )
      IF (IFLAG .EQ. 3)  THEN
        WRITE (16, *) X(I), FITFN
      END IF
      END DO

      END IF



      IF ( INCLUDEDIS ) THEN

*     DATASET 7  [H1,  Nucl.Phys.B470(96)3]
*       X = Q2,   Y = F2
*       X_BJ = 0.002

      XI = 0.002d0
      N = 9
      CALL READDATA ('DATASET7', N, X, Y, DY)
      DO I= 1, N
      Q2 = X(I)
      CALL F2F
      FITFN = F2(1)
      CHISQ = CHISQ + ( (FITFN - Y(I))**2 / DY(I)**2 )
      IF (IFLAG .EQ. 3)  THEN
        WRITE (16, *) X(I), FITFN
      END IF
      END DO

      IF (IFLAG .EQ. 3)  THEN
*       empty line for separation of different fit lines
        WRITE (16, *) 
      END IF

*     DATASET 8  [H1,  Nucl.Phys.B470(96)3]
*       X = Q2,   Y = F2
*       X_BJ = 0.0005

      XI = 0.0005d0
      N = 8
      CALL READDATA ('DATASET8', N, X, Y, DY)
      DO I= 1, N
      Q2 = X(I)
      CALL F2F
      FITFN = F2(1)
      CHISQ = CHISQ + ( (FITFN - Y(I))**2 / DY(I)**2 )
      IF (IFLAG .EQ. 3)  THEN
        WRITE (16, *) X(I), FITFN
      END IF
      END DO

      IF (IFLAG .EQ. 3)  THEN
*       empty line for separation of different fit lines
        WRITE (16, *) 
      END IF

*     DATASET 9  [H1,  Nucl.Phys.B470(96)3]
*       X = X_BJ,   Y = F2
*       Q2 = 8.5 GeV^2

      Q2 = 8.5d0
      N = 9
      CALL READDATA ('DATASET9', N, X, Y, DY)
      DO I= 1, N
      XI = X(I)
      CALL F2F
      FITFN = F2(1)
      CHISQ = CHISQ + ( (FITFN - Y(I))**2 / DY(I)**2 )
      IF (IFLAG .EQ. 3)  THEN
        WRITE (16, *) X(I), FITFN
      END IF
      END DO

      IF (IFLAG .EQ. 3)  THEN
*       empty line for separation of different fit lines
        WRITE (16, *) 
      END IF


*     DATASET 0  [H1,  Nucl.Phys.B470(96)3]
*       X = X_BJ,   Y = F2
*       Q2 = 15 GeV^2

      Q2 = 15.0d0
      N = 9
      CALL READDATA ('DATASET0', N, X, Y, DY)
      DO I= 1, N
      XI = X(I)
      CALL F2F
      FITFN = F2(1)
      CHISQ = CHISQ + ( (FITFN - Y(I))**2 / DY(I)**2 )
      IF (IFLAG .EQ. 3)  THEN
        WRITE (16, *) X(I), FITFN
      END IF
      END DO

      ENDIF

      F = CHISQ

*  Final commands executed after fit is done
      IF (IFLAG .EQ. 3)  THEN
         WRITE (*,*) 'CHISQ = ', CHISQ
         CLOSE (16)
      END IF

      RETURN
      END
C     ****


      SUBROUTINE READDATA (FILENAME, N, X, Y, DY)

      CHARACTER  FILENAME*8
      INTEGER I, N, NMAX
      PARAMETER (NMAX=40)
      DOUBLE PRECISION STAT, SYST
      DOUBLE PRECISION X(NMAX), Y(NMAX), DY(NMAX)

      OPEN (UNIT = 15, FILE = FILENAME, STATUS = 'OLD')
      DO 11 I=1,N
      READ (15, *) X(I), Y(I), STAT, SYS
*      Adding stat. i syst. errors in quadrature:
 11   DY(I) = SQRT( STAT*STAT + SYS*SYS )
      CLOSE (15)

      RETURN
      END
C     ****
