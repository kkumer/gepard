
#include "mathlink.h"
#include "header.h"
#include <string.h>


void MinuitInit(int a) {
        int i;
        struct dblcomplex nc;

        fitinit_();

          MLPutFunction(stdlink, "List", 2);
          /* FIXME: 4 below should be mt_.nmtsexp ! */
          MLPutFunction(stdlink, "List", mt_.nmts + 4 + 1); /* passing -t points */
            for (i = 0; i < (mt_.nmts + 4 + 1); i++){
                    MLPutReal(stdlink, mts_.mts[i]);
            }
          MLPutFunction(stdlink, "List", contour_.npts); /* passing MB points */
            for (i = 0; i < contour_.npts; i++){
                    nc = npoints_.n[i][0];
                    MLPutFunction(stdlink, "Complex", 2);
                      MLPutReal(stdlink, nc.dr - 1);
                      MLPutReal(stdlink, nc.di);
            }

        return;
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

void MinuitStatus(void){
        long int npari, nparx, istat;
        double fmin, fedm, errdef;

        mnstat_(&fmin, &fedm, &errdef, &npari, &nparx, &istat);

        MLPutFunction(stdlink, "List", 2);
          MLPutReal(stdlink, fmin);
          MLPutInteger(stdlink, istat);

        return;
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

void cffH(void) {

        int i, j, xint, dttype;
        double xreal, xi, xr;
        struct dblcomplex xc, nc;
        double args[10], mt;
        long int nargs;
        const char *fname, *sname;
        long int evoli=1, evolj=1;

                       

strcpy(parchr_, "12345678901DVCS  SINGLET");

readpar_();
par_.par[49] = 0.0;
par_.par[50] = 0.0;

getRealList(6, args);

kinematics_.xi = args[0];
mt = - args[1];
kinematics_.q2 = args[2];
nqs_.nqs = 2;
qs_.qs[0] = kinematics_.q2;

parflt_.q02 = args[3];

parint_.nf = (int) args[4];
parint_.p = (int) args[5];

astrong_.asp[parint_.p] = 0.05;

init_();

/*initgpd_(&mt);*/


        MLPutFunction(stdlink, "EvaluatePacket", 1);
          MLPutFunction(stdlink, "Set", 2);
            MLPutSymbol(stdlink, "tValues");
            MLPutFunction(stdlink, "List", 1);
              MLPutReal(stdlink, mt);
        MLEndPacket(stdlink);
        MLNextPacket(stdlink);
        MLNewPacket(stdlink);
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

        MLEvaluate(stdlink, "GPD[InitialValues]");

        MLNextPacket(stdlink);
        MLGetFunction(stdlink, &fname, &nargs); /* fname = List */
            for (j = 0; j < 1; j++){
                    for (i = 0; i < contour_.npts; i++){
                            /* quark GPD */
                            MLGetFunction(stdlink, &fname, &nargs); /* fname = Complex */
                              MLGetReal(stdlink, &xr);
                              MLGetReal(stdlink, &xi);
                            xc.dr = xr;
                            xc.di = xi;
                            hgrid_.hgrid[0][i][j] = xc;
                            /* gluon GPD */
                            MLGetFunction(stdlink, &fname, &nargs); /* fname = Complex */
                              MLGetReal(stdlink, &xr);
                              MLGetReal(stdlink, &xi);
                            xc.dr = xr;
                            xc.di = xi;
                            hgrid_.hgrid[1][i][j] = xc;
                    }
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



void initgpdmma_(double mt) {

        int i, j;
        double xr, xi;
        struct dblcomplex nc, xc;
        long int effacc, nargs;
        const char *fname;

/* calling Mathematica function GPD[{params}]  */

        MLPutFunction(stdlink, "EvaluatePacket", 1);
          MLPutFunction(stdlink, "GPD", 1);
          MLPutFunction(stdlink, "List", NPARMAX); /* passing parameters */
            for (i = 0; i < NPARMAX; i++){
                    MLPutReal(stdlink, par_.par[i]);
            }
        MLEndPacket(stdlink);

/* getting returned result and writing it to COMMON block  */

        MLNextPacket(stdlink);
        MLGetFunction(stdlink, &fname, &nargs); /* fname = List */
            /* FIXME: 4 below should be mt_.nmtsexp ! */
            for (j = 0; j < (mt_.nmts + 4 + 1); j++){
                    for (i = 0; i < contour_.npts; i++){
                            /* quark GPD */
                            MLGetFunction(stdlink, &fname, &nargs); /* fname = Complex */
                              MLGetReal(stdlink, &xr);
                              MLGetReal(stdlink, &xi);
                            xc.dr = xr;
                            xc.di = xi;
                            hgrid_.hgrid[0][i][j] = xc;
                            /* gluon GPD */
                            MLGetFunction(stdlink, &fname, &nargs); /* fname = Complex */
                              MLGetReal(stdlink, &xr);
                              MLGetReal(stdlink, &xi);
                            xc.dr = xr;
                            xc.di = xi;
                            hgrid_.hgrid[1][i][j] = xc;
                    }
            }
        MLEndPacket(stdlink);

        return;
};

int main(int argc, char *argv[]) {
           return MLMain(argc, argv);
}

