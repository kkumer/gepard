
#include "mathlink.h"
#include "header.h"
#include <string.h>

/*  Commands that call Minuit subroutines via Fortran intermediaries */

void MinuitInit(int a) {
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

void spliceparchr(char *out, int start, int end) {
        int i, k;

        for (i = start, k=0; i < end; i++, k++)
                out[k] = parchr_[i];

        return;
}

char *GepardInitInternal(int speed, int p, char *scheme, char *ansatz) {
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


int MinuitSetParameter(long int id, char *pnam, double vstrt, double step, double lo, double hi){
        long int tlen;

        tlen = strlen(pnam);
        mpar_(&id, &tlen, pnam, &vstrt, &step, &lo, &hi, tlen);

        return 0;
};

int MinuitCommand(char *cmd){
        long int tlen, ierflg;

        tlen = strlen(cmd);
        mcom_(&tlen, cmd, &ierflg,tlen);

        return ierflg;
};

void MinuitGetParameter(long int id){
        long int ivarbl;
        double val, error;

        getpar_(&id, &val, &error, &ivarbl);

        MLPutFunction(stdlink, "List", 3);
          MLPutReal(stdlink, val);
          MLPutReal(stdlink, error);
          MLPutInteger(stdlink, ivarbl);

        return;
};

/*  Commands that directly call Minuit subroutines */

void MinuitStatus(void){
        long int npari, nparx, istat;
        double fmin, fedm, errdef;

        mnstat_(&fmin, &fedm, &errdef, &npari, &nparx, &istat);

        MLPutFunction(stdlink, "List", 2);
          MLPutReal(stdlink, fmin);
          MLPutInteger(stdlink, istat);

        return;
};

void MinuitCovarianceMatrix(int adim){
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

/*  Other commands (that don't communicate with Minuit) */

void getRealList(int n, double list[10]){
        int j, xint, dttype;
        double xreal;
        long int nargs;
        const char *fname, *sname;

/* getting input via MathLink.  */

        /*MLGetFunction(stdlink, &fname, &nargs); [> fname = List <]*/
            for (j = 0; j < n; j++){

                dttype=MLGetType(stdlink); 
                
                switch (MLGetType(stdlink)) {
                        case MLTKINT:   /* Mma sent integer number*/
                                MLGetInteger(stdlink, &xint);
                                xreal = (double) xint;
                                break;
                        case MLTKREAL:  /* Mma sent real number */
                                MLGetReal(stdlink, &xreal);
                                break;
                       /*  case MLTKSYM:
                                MLGetSymbol(stdlink, &sname);
                                MLPutSymbol(stdlink, "err-sym");
                                return;
                       */
                        default:  /* Return an error */
                                MLPutSymbol(stdlink, "errnan");
                                return;
                }
                list[j] = xreal;
        }
        return;
}

void cffHInternal(double xi, double t, double q2, double q02, int speed, int p, char *scheme, char *ansatz) {

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


void getmbgpdmma_(void) {

        int i, j;
        double xr, xi;
        struct dblcomplex nc, xc;
        long int effacc, nargs;
        const char *fname;

/* calling Mathematica function GPD[{params}]  */

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

int main(int argc, char *argv[]) {
           return MLMain(argc, argv);
}

