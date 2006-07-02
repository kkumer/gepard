/*  
    file: rnnlo.c                           
    C wrapper which calls FORTRAN subroutine on behalf of MathLink

*/

#include "mathlink.h"
#include "gepard.h"

#pragma aux RNNLOF "^"  /* For Watcom compiler, GCC ignores it */
void RNNLOF(int *nf, Doublecomplex *n, 
                Doublecomplex pr[2][2][2], 
                Doublecomplex r1proj[2][2][2][2], 
                Doublecomplex r2proj[2][2][2][2]);


void rnnlo(void)
{
        int i, j, k, l;
        int mlint, nf, dttype;
        double mlreal, mlr, mli;
        Doublecomplex n;
        Doublecomplex pr[2][2][2]; 
        Doublecomplex r1proj[2][2][2][2]; 
        Doublecomplex r2proj[2][2][2][2];
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
        RNNLOF(&nf, &n, pr, r1proj, r2proj);

/* returning result to Mathematica via MathLink 
 * It's formated as a matrix 
 *            ++ +-
 *            -+ --
 * of matrices +- = P+.R1.P- etc.
 * Second such matrix is with R1->R2 and the two
 * are returned as a final list of two.
 * FIXME: Use mlputcomplex2by2matrix; return also PR, PROJ2
 */
        MLPutFunction(stdlink, "List", 2);
        for (k=0; k<2; k++) {
            MLPutFunction(stdlink, "List", 2);
            for (l=0; l<2; l++) {
                /* Now comes matrix P.R1.P */
                MLPutFunction(stdlink, "List", 2);
                for (i=0; i<2; i++) {
                    MLPutFunction(stdlink, "List", 2);
                    for (j=0; j<2; j++) 
                        mlputcomplex(stdlink, r1proj[j][i][l][k]);
                }
            }
        }
        return;
}
