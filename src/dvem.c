/****h* gepard/int2F1.c
*  FILE DESCRIPTION
*    Mathematica - DVEM interface (C)
*
*    $Id:$
******
*/

#include "mathlink.h"
#include "header.h"
/* #include <string.h>
*/



void int2f1(void) {

        int xint, dttype;
        double xreal, xi, xr;
        float mu, nu;
        float ind, acc;
        struct dblcomplex j, k;
        struct dblcomplex result;
        int nargs;
        const char *fname;


/* getting input mu via MathLink.  */

        dttype=MLGetType(stdlink); 
        switch (MLGetType(stdlink)) {
                case MLTKREAL:  /* Mma sent real number */
                        MLGetReal(stdlink, &xreal);
                        mu=xreal;
                        break;
                default:  /* Return an error */
                        MLPutSymbol(stdlink, "errnan");
                        return;
        }

/* getting input j via MathLink.  */

        dttype=MLGetType(stdlink); 
        switch (MLGetType(stdlink)) {
                case MLTKINT:   /* Mma sent integer number*/
                        MLGetInteger(stdlink, &xint);
                        j.dr=xint;
                        j.di=0.0;
                        break;
                case MLTKREAL:  /* Mma sent real number */
                        MLGetReal(stdlink, &xreal);
                        j.dr=xreal;
                        j.di=0.0;
                        break;
                case MLTKFUNC:  /* Mma sent Complex[] */
                        /* FIXME: Check if it realy is Complex[] */
                        MLGetFunction(stdlink, &fname, &nargs);
                        MLGetReal(stdlink, &xr);
                        MLGetReal(stdlink, &xi);
                        j.dr=xr;
                        j.di=xi;
                        break;
                default:  /* Return an error */
                        MLPutSymbol(stdlink, "errnan");
                        return;
        }
                       
/* getting input nu via MathLink.  */

        dttype=MLGetType(stdlink); 
        switch (MLGetType(stdlink)) {
                case MLTKREAL:  /* Mma sent real number */
                        MLGetReal(stdlink, &xreal);
                        nu=xreal;
                        break;
                default:  /* Return an error */
                        MLPutSymbol(stdlink, "errnan");
                        return;
        }

/* getting input k via MathLink.  */

        dttype=MLGetType(stdlink); 
        switch (MLGetType(stdlink)) {
                case MLTKINT:   /* Mma sent integer number*/
                        MLGetInteger(stdlink, &xint);
                        k.dr=xint;
                        k.di=0.0;
                        break;
                case MLTKREAL:  /* Mma sent real number */
                        MLGetReal(stdlink, &xreal);
                        k.dr=xreal;
                        k.di=0.0;
                        break;
                case MLTKFUNC:  /* Mma sent Complex[] */
                        /* FIXME: Check if it realy is Complex[] */
                        MLGetFunction(stdlink, &fname, &nargs);
                        MLGetReal(stdlink, &xr);
                        MLGetReal(stdlink, &xi);
                        k.dr=xr;
                        k.di=xi;
                        break;
                default:  /* Return an error */
                        MLPutSymbol(stdlink, "errnan");
                        return;
        }

/* getting input ind via MathLink.  */

        dttype=MLGetType(stdlink); 
        switch (MLGetType(stdlink)) {
                case MLTKINT:   /* Mma sent integer number*/
                        MLGetInteger(stdlink, &xint);
                        ind=xint;
                        break;
                default:  /* Return an error */
                        MLPutSymbol(stdlink, "errnan");
                        return;
        }
                       
/* getting input acc via MathLink.  */

        dttype=MLGetType(stdlink); 
        switch (MLGetType(stdlink)) {
                case MLTKINT:   /* Mma sent integer number*/
                        MLGetInteger(stdlink, &xint);
                        acc=xint;
                        break;
                default:  /* Return an error */
                        MLPutSymbol(stdlink, "errnan");
                        return;
        }

/* calling FORTRAN double complex function  */
        int2f1f_(&mu, &j, &nu, &k, &ind, &acc, &result);

/* returning result via MathLink  */
        MLPutFunction(stdlink, "Complex", 2);
        MLPutReal(stdlink, result.dr);
        MLPutReal(stdlink, result.di);

        return;
};

void cdvemf_(int sgntr, struct dblcomplex *j, struct dblcomplex *k, struct dblcomplex (*mcu)[]);

void cdvem(void) {

        int xint, dttype, nf, sgntr;
        double xreal, xi, xr;
        double rgpdf2, rdaf2, rr2;
        struct dblcomplex j, k, z;
        struct dblcomplex mcu[3];
        int i, nargs;
        const char *fname;

/* getting input sgntr via MathLink.  */

        dttype=MLGetType(stdlink); 
        switch (MLGetType(stdlink)) {
                case MLTKINT:   /* Mma sent integer number*/
                        MLGetInteger(stdlink, &xint);
                        sgntr=xint;
                        break;
                case MLTKREAL:  /* Mma sent real number */
                        MLGetReal(stdlink, &xreal);
                        sgntr = (int) xreal ;
                        break;
                default:  /* Return an error */
                        MLPutSymbol(stdlink, "errnan");
                        return;
        }



/* getting input j via MathLink.  */

        dttype=MLGetType(stdlink); 
        switch (MLGetType(stdlink)) {
                case MLTKINT:   /* Mma sent integer number*/
                        MLGetInteger(stdlink, &xint);
                        j.dr=xint;
                        j.di=0.0;
                        break;
                case MLTKREAL:  /* Mma sent real number */
                        MLGetReal(stdlink, &xreal);
                        j.dr=xreal;
                        j.di=0.0;
                        break;
                case MLTKFUNC:  /* Mma sent Complex[] */
                        /* FIXME: Check if it realy is Complex[] */
                        MLGetFunction(stdlink, &fname, &nargs);
                        MLGetReal(stdlink, &xr);
                        MLGetReal(stdlink, &xi);
                        j.dr=xr;
                        j.di=xi;
                        break;
                default:  /* Return an error */
                        MLPutSymbol(stdlink, "errnan");
                        return;
        }
                       
/* getting input k via MathLink.  */

        dttype=MLGetType(stdlink); 
        switch (MLGetType(stdlink)) {
                case MLTKINT:   /* Mma sent integer number*/
                        MLGetInteger(stdlink, &xint);
                        k.dr=xint;
                        k.di=0.0;
                        break;
                case MLTKREAL:  /* Mma sent real number */
                        MLGetReal(stdlink, &xreal);
                        k.dr=xreal;
                        k.di=0.0;
                        break;
                case MLTKFUNC:  /* Mma sent Complex[] */
                        /* FIXME: Check if it realy is Complex[] */
                        MLGetFunction(stdlink, &fname, &nargs);
                        MLGetReal(stdlink, &xr);
                        MLGetReal(stdlink, &xi);
                        k.dr=xr;
                        k.di=xi;
                        break;
                default:  /* Return an error */
                        MLPutSymbol(stdlink, "errnan");
                        return;
        }

/* getting input nf via MathLink.  */

        dttype=MLGetType(stdlink); 
        switch (MLGetType(stdlink)) {
                case MLTKINT:   /* Mma sent integer number*/
                        MLGetInteger(stdlink, &xint);
                        nf=xint;
                        break;
                case MLTKREAL:  /* Mma sent real number */
                        MLGetReal(stdlink, &xreal);
                        nf = (int) xreal ;
                        break;
                default:  /* Return an error */
                        MLPutSymbol(stdlink, "errnan");
                        return;
        }


/* getting input rgpdf2 via MathLink.  */

        dttype=MLGetType(stdlink); 
        switch (MLGetType(stdlink)) {
                case MLTKINT:   /* Mma sent integer number*/
                        MLGetInteger(stdlink, &xint);
                        rgpdf2=xint;
                        break;
                case MLTKREAL:  /* Mma sent real number */
                        MLGetReal(stdlink, &xreal);
                        rgpdf2=xreal;
                        break;
                default:  /* Return an error */
                        MLPutSymbol(stdlink, "errnan");
                        return;
        }

/* getting input rdaf2 via MathLink.  */

        dttype=MLGetType(stdlink); 
        switch (MLGetType(stdlink)) {
                case MLTKINT:   /* Mma sent integer number*/
                        MLGetInteger(stdlink, &xint);
                        rdaf2=xint;
                        break;
                case MLTKREAL:  /* Mma sent real number */
                        MLGetReal(stdlink, &xreal);
                        rdaf2=xreal;
                        break;
                default:  /* Return an error */
                        MLPutSymbol(stdlink, "errnan");
                        return;
        }

/* getting input rr2 via MathLink.  */

        dttype=MLGetType(stdlink); 
        switch (MLGetType(stdlink)) {
                case MLTKINT:   /* Mma sent integer number*/
                        MLGetInteger(stdlink, &xint);
                        rr2=xint;
                        break;
                case MLTKREAL:  /* Mma sent real number */
                        MLGetReal(stdlink, &xreal);
                        rr2=xreal;
                        break;
                default:  /* Return an error */
                        MLPutSymbol(stdlink, "errnan");
                        return;
        }

   /* putting scales in common block */

        parflt_.rgpdf2 = rgpdf2;
        parflt_.rdaf2 = rdaf2;
        parflt_.rr2 = rr2;

        parint_.nf = nf;

/* calling FORTRAN subroutine  */
        cdvemf_(&sgntr, &j, &k, &mcu);

/* returning result via MathLink  */
            MLPutFunction(stdlink, "List", 3); /* passing MB points */
              for (i = 0; i < 3; i++){
                z = mcu[i];
				MLPutFunction(stdlink, "Complex", 2);
				MLPutReal(stdlink, z.dr);
				MLPutReal(stdlink, z.di);
              }

        return;
};

int main(int argc, char *argv[]) {
           return MLMain(argc, argv);
}

