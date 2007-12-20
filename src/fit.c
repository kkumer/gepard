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
                    nc = npoints_.n[i];
                    MLPutFunction(stdlink, "Complex", 2);
                      MLPutReal(stdlink, nc.dr - 1);
                      MLPutReal(stdlink, nc.di);
            }

        return;
};
/******/

/****f* fit.c/spliceparchr
*  NAME
*     spliceparchr  --   Takes a substring of PARCHR block
*  DESCRIPTION
*     Since Fortran CHARACTER common blocks are recognized from
*     C as just a char array, to access members of this block
*     we extract part of the array.
*  SYNOPSIS
*/

void spliceparchr(char *out, int start, int end)

/*
*  INPUTS 
*           start -- (Fortran) index of first character we want - 1
*             end -- (Fortran) index of last character we want
*  PARENTS
*     GepardInitInternal
*  SOURCE
*/
{
        int i, k;

        for (i = start, k=0; i < end; i++, k++)
                out[k] = parchr_[i];

        return;
}
/******/

/****f* fit.c/GepardInitInternal
*  NAME
*     GepardInitInternal  --   Optional overriding of GEPARD.INI
*  DESCRIPTION
*     Performs initialization of GeParD parameters. First reads default values
*     from GEPARD.INI, which are then overriden by the function arguments with
*     the same name if they are positive numbers or strings different then
*     'DFLT'.
*  SYNOPSIS
*/

char *GepardInitInternal(int speed, int p, char *scheme, char *ansatz)
 
/*
*  INPUTS 
*           speed -- SPEED
*               p -- P
*          scheme -- SCHEME
*          ansatz -- ANSATZ
*  PARENTS
*     GepardInit
*  CHILDREN
*     READPAR, spliceparchr
*  SOURCE
*/
{
        const char dflt[4] = "DFLT"; 
        char inischeme[6], iniansatz[7];

        readpar_();

        if (speed >= 0) /* override GEPARD.INI */
          parint_.speed = speed;
        if (p >= 0) /* override GEPARD.INI */
          parint_.p = p;


        spliceparchr(inischeme, 0, 5);
        spliceparchr(iniansatz, 5, 11);

        strcat(scheme, "     "); /* want len(scheme) at least 5 */ 
        strcat(ansatz, "      "); /* want len(ansatz) at least 6 */ 

        if (strncmp(scheme, dflt, 4) != 0)   /* override GEPARD.INI */  
          strncpy(inischeme, scheme, 5);

        if (strncmp(ansatz, dflt, 4) != 0)   /* override GEPARD.INI */  
          strncpy(iniansatz, ansatz, 6);

        strcpy(parchr_, "");
        strncat(parchr_, inischeme, 5);
        strncat(parchr_, iniansatz, 6);
        strcat(parchr_, "DVCS  SINGLET   ");

        /*FALLBACK: strcpy(parchr_, "MSBARMMA   DVCS  SINGLET");*/
        return parchr_;
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

void cffHInternal(double xi, double t, double q2, double q02, int speed, int p, char *scheme, char *ansatz) 

/*
*  INPUTS 
*         xi -- XI;   t -- DEL2;  q2 -- Q2;  q02 -- Q02;
*      speed -- SPEED; p -- P; scheme -- SCHEME; ansatz -- ANSATZ
*  CHILDREN
*     GepardInitInternal, INIT, EVOLC, CFFF, MLPut*
*  SOURCE
*/
{
        int i, j, xint, dttype;
        double xreal, xim, xre;
        struct dblcomplex xc, nc;
        double args[10], mt;
        long int nargs;
        const char *fname, *sname;
        long int evoli=1, evolj=1;


        GepardInitInternal(speed, p, scheme, ansatz);

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
                      nc = npoints_.n[i];
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

cfff_();

/* returning result via MathLink  */
        MLPutFunction(stdlink, "Complex", 2);
        MLPutReal(stdlink, cff_.cff[parint_.p].dr);
        MLPutReal(stdlink, cff_.cff[parint_.p].di);

        return;
};
/******/

/****f* fit.c/cffEInternal
*  NAME
*     cffEInternal  --  calculates CFF E
*  SYNOPSIS
*/

void cffEInternal(double xi, double t, double q2, double q02, int speed, int p, char *scheme, char *ansatz) 

/*
*  INPUTS 
*         xi -- XI;   t -- DEL2;  q2 -- Q2;  q02 -- Q02;
*      speed -- SPEED; p -- P; scheme -- SCHEME; ansatz -- ANSATZ
*  CHILDREN
*     GepardInitInternal, INIT, EVOLC, CFFFE, MLPut*
*  SOURCE
*/
{
        int i, j, xint, dttype;
        double xreal, xim, xre;
        struct dblcomplex xc, nc;
        double args[10], mt;
        long int nargs;
        const char *fname, *sname;
        long int evoli=1, evolj=1;


        GepardInitInternal(speed, p, scheme, ansatz);

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
                      nc = npoints_.n[i];
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

cfff_();

/* returning result via MathLink  */
        MLPutFunction(stdlink, "Complex", 2);
        MLPutReal(stdlink, cff_.cffe[parint_.p].dr);
        MLPutReal(stdlink, cff_.cffe[parint_.p].di);

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

        int i, j;
        double xr, xi;
        struct dblcomplex nc, xc;
        long int effacc, nargs;
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

/****p* fit.c/main
*  NAME
*     main  --  main MathLink program
*  SOURCE
*/

int main(int argc, char *argv[]) {
           return MLMain(argc, argv);
}
/******/
