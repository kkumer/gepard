/****h* gepard/int2F1.c
*  FILE DESCRIPTION
*    Mathematica - DVEM interface (C)
*
*    $Id:$
******
*/

#include "mathlink.h"
/* #include "header.h"
* #include <string.h>
*/

struct dblcomplex {
        double dr; 
        double di;
};


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


int main(int argc, char *argv[]) {
           return MLMain(argc, argv);
}

