/****h* gepard/fit.c
*  FILE DESCRIPTION
*    Mathematica - Minuit interface (C)
*
*    $Id:$
******
*/

#include "mathlink.h"
#include "header.h"
#include <string.h>

/*  Commands that call Minuit subroutines via Fortran intermediaries */

/****f* fit.c/MinuitInit
*  NAME
*     MinuitInit  --   Initialization of fitting procedure
*  DESCRIPTION
*     Wrapper around FITINIT, which additionally returns to
*     Mathematica a list of complex
*     coordinates of points on Mellin-Barnes contour
*  SYNOPSIS
*/

void MinuitInit(int a) 

/*
*  INPUTS 
*              a  --  does nothing
*  PARENTS
*     MinuitSetParameter
*  CHILDREN
*     FITINIT, MLPut*
*  BUGS    
*     should be without arguments?
*  SOURCE
*/
{
        int i;
        struct dblcomplex nc;

        fitinit_();

          MLPutFunction(stdlink, "List", contour_.npts); /* passing MB points */
            for (i = 0; i < contour_.npts; i++){
                    nc = npoints_.n[i][0];
                    MLPutFunction(stdlink, "Complex", 2);
                      MLPutReal(stdlink, nc.dr - 1);
                      MLPutReal(stdlink, nc.di);
            }

        return;
};
/******/


/****f* fit.c/GepardInitInternal
*  NAME
*     GepardInitInternal  --   Optional overriding of GEPARD.INI
*  DESCRIPTION
*     Performs initialization of GeParD parameters. First it reads default values
*     from GEPARD.INI, which are then overriden by the function arguments with
*     the same name if they are positive numbers or strings different then
*     'DFLT'.
*  SYNOPSIS
*/

int *GepardInitInternal(int speed, int p, char *process,
        char *scheme, char *ansatz, char *datfile, char *outfile)
 
/*
*  INPUTS 
*           speed -- SPEED
*               p -- P
*          process -- PROCESS
*          scheme -- SCHEME
*          ansatz -- ANSATZ
*          datfile -- DATFILE
*          outfile -- OUTFILE
*  PARENTS
*     GepardInit
*  CHILDREN
*     PARINIT
*  SOURCE
*/
{
        long int processlen, schemelen, ansatzlen, datfilelen, outfilelen;

        processlen = strlen(process);
        schemelen = strlen(scheme);
        ansatzlen = strlen(ansatz);
        datfilelen = strlen(datfile);
        outfilelen = strlen(outfile);

        parinit_(&speed, &p, &processlen, process,
                        &schemelen, scheme, &ansatzlen, ansatz,
                        &datfilelen, datfile, &outfilelen, outfile, 
               processlen, schemelen, ansatzlen, datfilelen, outfilelen);

        return 0;
};
/******/


/****f* fit.c/MinuitSetParameter
*  NAME
*     MinuitSetParameter  --   C wrapper for MPAR
*  SYNOPSIS
*/

int MinuitSetParameter(long int id, char *pnam, double vstrt, double step, double lo, double hi)

/*
*  INPUTS 
*                  id  --  parameter number
*                pnam  --  parameter name
*               vstrt  --  starting value of parameter
*                step  --  starting step size or approx. parameter error
*                  lo  --  lower bound on parameter
*                  hi  --  upper bound on parameter
*                          (if lo=hi=0 parameter is considered unbounded)
*  CHILDREN
*            MPAR
*  SOURCE
*/
{        
        long int tlen;

        tlen = strlen(pnam);
        mpar_(&id, &tlen, pnam, &vstrt, &step, &lo, &hi, tlen);

        return 0;
};
/******/

/****f* fit.c/MinuitCommand
*  NAME
*     MinuitCommand  --  C wrapper for MCOM
*  SYNOPSIS
*/

int MinuitCommand(char *cmd)

/*
*  INPUTS 
*           cmd -- Minuit command as character string
*  CHILDREN
*            MCOM
*  SOURCE
*/
{
        long int tlen, ierflg;

        tlen = strlen(cmd);
        mcom_(&tlen, cmd, &ierflg,tlen);
        flushout_();

        return ierflg;
};
/******/

/****f* fit.c/MinuitGetParameter
*  NAME
*     MinuitGetParameter  --  C wrapper for GETPAR
*  DESCRIPTION
*     It returns current value of parameter and
*     its error, and internal parameter number to Mathematica
*  SYNOPSIS
*/

void MinuitGetParameter(long int id)

/*
*  INPUTS 
*           id -- MINUIT's parameter number
*  CHILDREN
*            GETPAR, MLPut*
*  SOURCE
*/
{
        long int ivarbl;
        double val, error;

        getpar_(&id, &val, &error, &ivarbl);

        MLPutFunction(stdlink, "List", 3);
          MLPutReal(stdlink, val);
          MLPutReal(stdlink, error);
          MLPutInteger(stdlink, ivarbl);

        return;
};
/******/

/****f* fit.c/MinuitContour
*  NAME
*     MinuitContour  --  C wrapper for MNCONT via CONTOUR
*  DESCRIPTION
*     It returns list of contour coordinates
*  SYNOPSIS
*/

void MinuitContour(long int num1, long int num2, long int npt)

/*
*  INPUTS 
*                num1  --  first parameter number (external)
*                num2  --  second parameter number (external)
*                 npt  --  number of points required on a contour
*  OUTPUT 
*                xpt  --  array of contour x-coordinates (x: NUM1)
*                ypt  --  array of contour y-coordinates (y: NUM2)
*  PARENTS
*      PlotMinuitContour
*  CHILDREN
*            MCONT
*  SOURCE
*/

{
        int nfound, i;
        double xpt[250], ypt[250];

        mcont_(&num1, &num2, &npt, xpt, ypt, &nfound);

        MLPutFunction(stdlink, "List", nfound);
            for (i = 0; i < nfound; i++){
                    MLPutFunction(stdlink, "List", 2);
                      MLPutReal(stdlink, xpt[i]);
                      MLPutReal(stdlink, ypt[i]);
            }

        return;
};
/******/

/*  Commands that directly call Minuit subroutines */

/****f* fit.c/MinuitStatus
*  NAME
*     MinuitStatus -- C Wrapper for MINUIT's MNSTAT
*  SYNOPSIS
*/

void MinuitStatus(void)

/*
*  OUTPUT 
*          rfmin  --  chi-square  
*          istat  --  quality of covariance matrix
*  CHILDREN
*     mnstat, MLPut*
*  SOURCE
*/
{
        long int npari, nparx, istat;
        double fmin, fedm, errdef;

        mnstat_(&fmin, &fedm, &errdef, &npari, &nparx, &istat);

        MLPutFunction(stdlink, "List", 2);
          MLPutReal(stdlink, fmin);
          MLPutInteger(stdlink, istat);

        return;
};
/******/

/****f* fit.c/MinuitCovarianceMatrix
*  NAME
*     MinuitCovarianceMatrix  --  C wrapper for MINUIT's MNEMAT
*  SYNOPSIS
*/

void MinuitCovarianceMatrix(int adim)

/*
*  INPUTS 
*            adim -- Number of variable parameters
*  CHILDREN
*           MNEMAT, MLPut*
*  SOURCE
*/
{
        int i,j;
        long int ndim=NPARMAX;
        double emat[NPARMAX][NPARMAX];

        mnemat_(emat, &ndim);

        MLPutFunction(stdlink, "List", adim);
          for (i = 0; i < adim; i++){
            MLPutFunction(stdlink, "List", adim);
            for (j = 0; j < adim; j++){
              MLPutReal(stdlink, emat[j][i]);
            }
          }
        return;
};
/******/

/*  Other commands (that don't communicate with Minuit) */


/****f* fit.c/cffHInternal
*  NAME
*     cffHInternal  --  calculates CFF H
*  SYNOPSIS
*/

void cffHInternal(double xi, double t, double q2, double q02, int speed, int p, 
        char *process, char *scheme, char *ansatz) 

/*
*  INPUTS 
*         xi -- XI;   t -- DEL2;  q2 -- Q2;  q02 -- Q02;
*      speed -- SPEED; p -- P; process -- PROCESS; scheme -- SCHEME; ansatz -- ANSATZ
*  CHILDREN
*     GepardInitInternal, SETPROC, INIT, EVOLC, CFFF, MLPut*
*  SOURCE
*/
{
        int i;
        double xreal;
        struct dblcomplex nc;
        int nargs;
        const char *fname;
        long int evoli=1, evolj=1, evols=2, evolt=3;
        char datfile[] = "DFLT";
        char outfile[] = "DFLT";


        GepardInitInternal(speed, p, process, scheme, ansatz, datfile, outfile);

        kinematics_.xi = xi;
        kinematics_.del2 = t;
        kinematics_.q2 = q2;
        nqs_.nqs = 2; /* Why not 1? */
        qs_.qs[0] = kinematics_.q2;
        parflt_.q02 = q02;

        init_();

        MLPutFunction(stdlink, "EvaluatePacket", 1);
          MLPutFunction(stdlink, "Set", 2);
            MLPutSymbol(stdlink, "jValues");
            MLPutFunction(stdlink, "List", contour_.npts); /* passing MB points */
              for (i = 0; i < contour_.npts; i++){
                      nc = npoints_.n[i][0];
                      MLPutFunction(stdlink, "Complex", 2);
                        MLPutReal(stdlink, nc.dr - 1);
                        MLPutReal(stdlink, nc.di);
              }
        MLEndPacket(stdlink);
        MLNextPacket(stdlink);
        MLNewPacket(stdlink);

        MLPutFunction(stdlink, "EvaluatePacket", 1);
          MLPutFunction(stdlink, "AllParameterValues", 0);
        MLEndPacket(stdlink);

        MLNextPacket(stdlink);
        MLGetFunction(stdlink, &fname, &nargs); /* fname = List */
            for (i = 0; i < NPARMAX; i++){
                    MLGetReal(stdlink, &xreal);
                    par_.par[i] = xreal;
            }
        MLEndPacket(stdlink);

evolc_(&evoli, &evolj);
evolc_(&evols, &evolj);
evolc_(&evolt, &evolj);

cfff_();

/* returning result via MathLink  */
        MLPutFunction(stdlink, "Complex", 2);
        MLPutReal(stdlink, cff_.cff[p].dr);
        MLPutReal(stdlink, cff_.cff[p].di);

        return;
};
/******/

/****f* fit.c/cffEInternal
*  NAME
*     cffEInternal  --  calculates CFF E
*  SYNOPSIS
*/

void cffEInternal(double xi, double t, double q2, double q02, int speed, int p, 
        char *process, char *scheme, char *ansatz) 

/*
*  INPUTS 
*         xi -- XI;   t -- DEL2;  q2 -- Q2;  q02 -- Q02;
*      speed -- SPEED; p -- P; process -- PROCESS; scheme -- SCHEME; ansatz -- ANSATZ
*  CHILDREN
*     GepardInitInternal, SETPROC, INIT, EVOLC, CFFFE, MLPut*
*  SOURCE
*/
{
        int i;
        double xreal;
        struct dblcomplex nc;
        int nargs;
        const char *fname;
        long int evoli=1, evolj=1, evols=2, evolt=3;
        char datfile[] = "DFLT";
        char outfile[] = "DFLT";


        GepardInitInternal(speed, p, process, scheme, ansatz, datfile, outfile);

        kinematics_.xi = xi;
        kinematics_.del2 = t;
        kinematics_.q2 = q2;
        nqs_.nqs = 2; /* Why not 1? */
        qs_.qs[0] = kinematics_.q2;
        parflt_.q02 = q02;

        init_();

        MLPutFunction(stdlink, "EvaluatePacket", 1);
          MLPutFunction(stdlink, "Set", 2);
            MLPutSymbol(stdlink, "jValues");
            MLPutFunction(stdlink, "List", contour_.npts); /* passing MB points */
              for (i = 0; i < contour_.npts; i++){
                      nc = npoints_.n[i][0];
                      MLPutFunction(stdlink, "Complex", 2);
                        MLPutReal(stdlink, nc.dr - 1);
                        MLPutReal(stdlink, nc.di);
              }
        MLEndPacket(stdlink);
        MLNextPacket(stdlink);
        MLNewPacket(stdlink);

        MLPutFunction(stdlink, "EvaluatePacket", 1);
          MLPutFunction(stdlink, "AllParameterValues", 0);
        MLEndPacket(stdlink);

        MLNextPacket(stdlink);
        MLGetFunction(stdlink, &fname, &nargs); /* fname = List */
            for (i = 0; i < NPARMAX; i++){
                    MLGetReal(stdlink, &xreal);
                    par_.par[i] = xreal;
            }
        MLEndPacket(stdlink);

evolc_(&evoli, &evolj);
evolc_(&evols, &evolj);
evolc_(&evolt, &evolj);

cfff_();

/* returning result via MathLink  */
        MLPutFunction(stdlink, "Complex", 2);
        MLPutReal(stdlink, cff_.cffe[p].dr);
        MLPutReal(stdlink, cff_.cffe[p].di);

        return;
};
/******/

/****f* fit.c/F2Internal
*  NAME
*     F2Internal  --  calculates singlet F2(x, Q2)
*  SYNOPSIS
*/

void F2Internal(double xbj, double q2, double q02, int speed, int p, 
        char *process, char *scheme, char *ansatz) 

/*
*  INPUTS 
*         xbj -- X_BJ;   q2 -- Q2;  q02 -- Q02;
*      speed -- SPEED; p -- P; process -- PROCESS; scheme -- SCHEME; ansatz -- ANSATZ
*  CHILDREN
*     GepardInitInternal, SETPROC, INIT, EVOLC, F2, MLPut*
*  SOURCE
*/
{
        int i;
        double xreal;
        struct dblcomplex nc;
        int nargs;
        const char *fname;
        long int evoli=1, evolj=1;
        char datfile[] = "DFLT";
        char outfile[] = "DFLT";


        GepardInitInternal(speed, p, process, scheme, ansatz, datfile, outfile);

        kinematics_.xi = xbj;
        kinematics_.del2 = 0;
        kinematics_.q2 = q2;
        nqs_.nqsdis = 2; /* Why not 1? */
        qs_.qsdis[0] = kinematics_.q2;
        parflt_.q02 = q02;

        init_();

        MLPutFunction(stdlink, "EvaluatePacket", 1);
          MLPutFunction(stdlink, "Set", 2);
            MLPutSymbol(stdlink, "jValues");
            MLPutFunction(stdlink, "List", contour_.npts); /* passing MB points */
              for (i = 0; i < contour_.npts; i++){
                      nc = npoints_.n[i][0];
                      MLPutFunction(stdlink, "Complex", 2);
                        MLPutReal(stdlink, nc.dr - 1);
                        MLPutReal(stdlink, nc.di);
              }
        MLEndPacket(stdlink);
        MLNextPacket(stdlink);
        MLNewPacket(stdlink);

        MLPutFunction(stdlink, "EvaluatePacket", 1);
          MLPutFunction(stdlink, "AllParameterValues", 0);
        MLEndPacket(stdlink);

        MLNextPacket(stdlink);
        MLGetFunction(stdlink, &fname, &nargs); /* fname = List */
            for (i = 0; i < NPARMAX; i++){
                    MLGetReal(stdlink, &xreal);
                    par_.par[i] = xreal;
            }
        MLEndPacket(stdlink);

        evolc_(&evoli, &evolj);

        if (parchr_[5] == 'M' && parchr_[6] == 'M' && parchr_[7] == 'A')  /* MMA */
          getmbgpdmma_();
        else                    /* FIT or other Fortran */
          getmbgpd_();

        for (i = 0; i < contour_.npts; i++){
                hgrid_.hgrid[0][i][0] = mbgpd_.mbgpd[0][i];
                hgrid_.hgrid[1][i][0] = mbgpd_.mbgpd[1][i];
        }


        f2f_();

/* returning result via MathLink. */
        MLPutReal(stdlink, f2_.f2[p]);

        return;
};
/******/



/****f* fit.c/getmbgpdmma_
*  NAME
*     getmbgpdmma_  --  get GPD values from Mathematica
*  DESCRIPTION
*     This Fortran-callable function uses MathLink to
*     evaluate functions GPD[{fitting parameters}, DEL2, XI] and 
*     PDF[...the same...] and thus uses GPD ansatz as defined
*     within Mathematica session. Mathematica returns values
*     on MB contour which are written into MBGPD common block.
*  SYNOPSIS
*/

void getmbgpdmma_(void) 

/*
*  PARENTS
*     FCN, CFFF
*  SOURCE
*/
{

        int i;
        double xr, xi;
        struct dblcomplex xc;
        int nargs;
        const char *fname;

/* calling Mathematica function GPD[{params}, DEL2, XI]  */

        MLPutFunction(stdlink, "EvaluatePacket", 1);
        if (parchr_[12] == 'I') /* DIS */
          MLPutFunction(stdlink, "PDF", 3);
        else                    /* DVCS */
          MLPutFunction(stdlink, "GPD", 3);
        
          MLPutFunction(stdlink, "List", NPARMAX); /* passing parameters */
            for (i = 0; i < NPARMAX; i++){
                    MLPutReal(stdlink, par_.par[i]);
            }
          MLPutReal(stdlink, kinematics_.del2);
          MLPutReal(stdlink, kinematics_.xi);
        MLEndPacket(stdlink);

/* getting returned result and writing it to COMMON block  */

        MLNextPacket(stdlink);
        MLGetFunction(stdlink, &fname, &nargs); /* fname = List */
            for (i = 0; i < contour_.npts; i++){
                    /* quark GPD */
                    MLGetFunction(stdlink, &fname, &nargs); /* fname = Complex */
                      MLGetReal(stdlink, &xr);
                      MLGetReal(stdlink, &xi);
                    xc.dr = xr;
                    xc.di = xi;
                    mbgpd_.mbgpd[0][i] = xc;
                    /* gluon GPD */
                    MLGetFunction(stdlink, &fname, &nargs); /* fname = Complex */
                      MLGetReal(stdlink, &xr);
                      MLGetReal(stdlink, &xi);
                    xc.dr = xr;
                    xc.di = xi;
                    mbgpd_.mbgpd[1][i] = xc;
            }
        MLEndPacket(stdlink);

        return;
};
/******/


/****f* fit.c/GetChiSquares
*  NAME
*     GetChiSquares --  Read chi-squares from common block
*  SYNOPSIS
*/

void GetChiSquares(void) 

/*
*  PARENTS
*     JustCali3, GepardFit
*  SOURCE
*/
{

        int i, nmax;

        nmax = chinblk_.chin + 1;

        MLPutFunction(stdlink, "List", nmax); /* passing parameters */
            for (i = 0; i < nmax; i++){
                    MLPutFunction(stdlink, "List", 2);
                    MLPutReal(stdlink, chisqblk_.chisq[i]);
                    MLPutInteger(stdlink, ndataptsblk_.ndatapts[i]);
            }
        MLEndPacket(stdlink);

        return;
};
/******/

/****f* fit.c/BCAInternal
*  NAME
*     BCAInternal  --  calculates beam charge asymmetry
*  SYNOPSIS
*/

void BCAInternal(double wavg, double q2avg, double phiin, 
        int speed, int p, char *process, char *scheme, char *ansatz) 

/*
*  INPUTS 
*         wavg -- average  W;  t -- DEL2;  q2avg -- average Q2;
*      speed -- SPEED; p -- P; process -- PROCESS; scheme -- SCHEME; ansatz -- ANSATZ
*  CHILDREN
*     GepardInitInternal, INIT, EVOLC, BCA, MLPut*
*  SOURCE
*/
{
        int i;
        double xreal, herabca;
        struct dblcomplex nc;
        int nargs; 
        long int fitflag=0;
        const char *fname;
        char datfile[] = "DFLT";
        char outfile[] = "DFLT";


        GepardInitInternal(speed, p, process, scheme, ansatz, datfile, outfile);

        init_();

        MLPutFunction(stdlink, "EvaluatePacket", 1);
          MLPutFunction(stdlink, "Set", 2);
            MLPutSymbol(stdlink, "jValues");
            MLPutFunction(stdlink, "List", contour_.npts); /* passing MB points */
              for (i = 0; i < contour_.npts; i++){
                      nc = npoints_.n[i][0];
                      MLPutFunction(stdlink, "Complex", 2);
                        MLPutReal(stdlink, nc.dr - 1);
                        MLPutReal(stdlink, nc.di);
              }
        MLEndPacket(stdlink);
        MLNextPacket(stdlink);
        MLNewPacket(stdlink);

        MLPutFunction(stdlink, "EvaluatePacket", 1);
          MLPutFunction(stdlink, "AllParameterValues", 0);
        MLEndPacket(stdlink);

        MLNextPacket(stdlink);
        MLGetFunction(stdlink, &fname, &nargs); /* fname = List */
            for (i = 0; i < NPARMAX; i++){
                    MLGetReal(stdlink, &xreal);
                    par_.par[i] = xreal;
            }
        MLEndPacket(stdlink);


        bca_(&fitflag, &wavg, &q2avg, &phiin, &herabca);

/* returning result via MathLink  */
        MLPutReal(stdlink, herabca);

        return;
};
/******/

/****p* fit.c/main
*  NAME
*     main  --  main MathLink program
*  SOURCE
*/

int main(int argc, char *argv[]) {
           return MLMain(argc, argv);
}
/******/
