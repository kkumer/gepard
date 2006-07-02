/*  
    file: evol.c                           
    C wrapper which calls FORTRAN subroutine on behalf of MathLink

*/

#include "mathlink.h"
#include "gepard.h"

#pragma aux EVOLF "^"  /* For Watcom compiler, GCC ignores it */
void EVOLF(int *nf, Doublecomplex *n, 
                double *r, Doublecomplex evola[2][2][3]);

void evol(void)
{
        int mlint, nf, dttype;
        int i, j, k;
        double mlreal, mlr, mli;
        double r;
        Doublecomplex n;
        Doublecomplex evola[2][2][3];
        long int nargs; const char *fname;


/* getting input via MathLink.  */

        /* First the integer nf */
/*
        dttype=MLGetType(stdlink); 
        
        switch (MLGetType(stdlink)) {
                case MLTKINT:
                        MLGetInteger(stdlink, &nf);
                        break;
                default:
                        MLPutSymbol(stdlink, "errnan");
                        return;
        }
 */
        /* First the int nf, which mma can send variously */

        dttype=MLGetType(stdlink); 
        
        switch (MLGetType(stdlink)) {
                case MLTKINT:   /* Mma sent integer number*/
                        MLGetInteger(stdlink, &mlint);
                        nf=mlint;
                        break;
                case MLTKREAL:  /* Mma sent real number */
                        MLGetReal(stdlink, &mlreal);
                        nf = (int) mlreal;
                        break;
                default:  /* Return an error */
                        MLPutSymbol(stdlink, "errnan");
                        return;
        }


        /* Then the complex n, which mma can send variously */

        dttype=MLGetType(stdlink); 
        
        switch (MLGetType(stdlink)) {
                case MLTKINT:   /* Mma sent integer number*/
                        MLGetInteger(stdlink, &mlint);
                        n.dr=mlint;
                        n.di=0.0;
                        break;
                case MLTKREAL:  /* Mma sent real number */
                        MLGetReal(stdlink, &mlreal);
                        n.dr=mlreal;
                        n.di=0.0;
                        break;
                case MLTKFUNC:  /* Mma sent Complex[] */
                        /* FIXME: Check if it really is Complex[] */
                        /*        Treat Rational[] correctly      */
                        MLGetFunction(stdlink, &fname, &nargs);
                        MLGetReal(stdlink, &mlr);
                        MLGetReal(stdlink, &mli);
                        n.dr=mlr;
                        n.di=mli;
                        break;
                default:  /* Return an error */
                        MLPutSymbol(stdlink, "errnan");
                        return;
        }

        /* Then the double precision R, which mma can send variously */

        dttype=MLGetType(stdlink); 
        
        switch (MLGetType(stdlink)) {
                case MLTKINT:   /* Mma sent integer number*/
                        MLGetInteger(stdlink, &mlint);
                        r=mlint;
                        break;
                case MLTKREAL:  /* Mma sent real number */
                        MLGetReal(stdlink, &r);
                        break;
                default:  /* Return an error */
                        MLPutSymbol(stdlink, "errnan");
                        return;
        }

/* calling FORTRAN  function  */
        EVOLF(&nf, &n, &r, evola);

/* returning result to Mathematica via MathLink 
 */
        MLPutFunction(stdlink, "List", 3);
        for (k=0; k<3; k++) {
            /* Now comes matrix \mathcal{A}^(k) */
            MLPutFunction(stdlink, "List", 2);
            for (i=0; i<2; i++) {
                MLPutFunction(stdlink, "List", 2);
                for (j=0; j<2; j++) 
                    mlputcomplex(stdlink, evola[j][i][k]);
            }
        }
        return;
}
