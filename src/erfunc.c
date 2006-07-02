/*  
    file: erfunc.c                           
    C wrapper which calls FORTRAN subroutine on behalf of MathLink

*/

#include "mathlink.h"
#include "gepard.h"

#pragma aux ERFUNCF "^"  /* For Watcom compiler, GCC ignores it */
void ERFUNCF(int *nf, Doublecomplex *n, double *r,
                Doublecomplex lamb[2],
                Doublecomplex erfunc1[2][2],
                Doublecomplex erfunc2[2][2]);

void erfunc(void)
{
        int i, j;
        int mlint, nf, dttype;
        double mlreal, mlr, mli;
        double r;
        Doublecomplex n;
        Doublecomplex lamb[2], erfunc1[2][2], erfunc2[2][2];
        long int nargs; const char *fname;


/* getting input via MathLink.  */

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

        /* Then the double precision r, which mma can send variously */

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
        ERFUNCF(&nf, &n, &r, lamb, erfunc1, erfunc2);

/* returning result to Mathematica via MathLink 
 */
        MLPutFunction(stdlink, "List", 2);
             mlputcomplex2by2matrix(stdlink, erfunc1);
             mlputcomplex2by2matrix(stdlink, erfunc2);

        return;
}
