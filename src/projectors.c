/*  
    file: projectors.c                           
    C wrapper which calls FORTRAN subroutine on behalf of MathLink

*/

#include "mathlink.h"
#include "gepard.h"

#pragma aux PROJECTORSF "^"  /* For Watcom compiler, GCC ignores it */
void PROJECTORSF(int *nf, Doublecomplex *n, Doublecomplex pr[2][2][2]);

void projectors(void)
{
        int mlint, nf, dttype;
        int i,j,l;
        double mlreal, mlr, mli;
        Doublecomplex n;
        Doublecomplex pr[2][2][2];
        long int nargs; const char *fname;


/* getting input via MathLink.  */

        /* First the integer nf */

        dttype=MLGetType(stdlink); 
        
        switch (MLGetType(stdlink)) {
                case MLTKINT:   /* Mma sent integer number*/
                        MLGetInteger(stdlink, &nf);
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

                       

/* calling FORTRAN  function  */
        PROJECTORSF(&nf, &n, pr);

/* returning result to Mathematica via MathLink  */
        MLPutFunction(stdlink, "List", 2);
        for (l=0; l<2; l++) {
            /* Now comes projectionmatrix P_l */
            MLPutFunction(stdlink, "List", 2);
            for (i=0; i<2; i++) {
                MLPutFunction(stdlink, "List", 2);
                for (j=0; j<2; j++) 
                    mlputcomplex(stdlink, pr[j][i][l]);
            }
        }
        return;
};

