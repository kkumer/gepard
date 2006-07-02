/*  
    file: cdvcs.c                           
    C wrapper which calls FORTRAN subroutine on behalf of MathLink

*/

#include "mathlink.h"
#include "gepard.h"

#pragma aux MSBARF "^"  /* For Watcom compiler, GCC ignores it */
void MSBARF(int *nf, Doublecomplex *j, 
                double *rf2, double *rr2, 
                Doublecomplex bigc0[2],
                Doublecomplex bigc1[2]);

void msbar(void)
{
        int mlint, nf, dttype;
        double mlreal, mlr, mli;
        double rf2, rr2;
        Doublecomplex j;
        Doublecomplex bigc0[2], bigc1[2];
        long int nargs; const char *fname;


/* getting input via MathLink.  */

        /* First the integer nf */

/*        dttype=MLGetType(stdlink); 
        
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

        /* Then the complex j, which mma can send variously */

        dttype=MLGetType(stdlink); 
        
        switch (MLGetType(stdlink)) {
                case MLTKINT:   /* Mma sent integer number*/
                        MLGetInteger(stdlink, &mlint);
                        j.dr=mlint;
                        j.di=0.0;
                        break;
                case MLTKREAL:  /* Mma sent real number */
                        MLGetReal(stdlink, &mlreal);
                        j.dr=mlreal;
                        j.di=0.0;
                        break;
                case MLTKFUNC:  /* Mma sent Complex[] */
                        /* FIXME: Check if it really is Complex[] */
                        /*        Treat Rational[] correctly      */
                        MLGetFunction(stdlink, &fname, &nargs);
                        MLGetReal(stdlink, &mlr);
                        MLGetReal(stdlink, &mli);
                        j.dr=mlr;
                        j.di=mli;
                        break;
                default:  /* Return an error */
                        MLPutSymbol(stdlink, "errnan");
                        return;
        }

        /* Then the double precision rf2, which mma can send variously */

        dttype=MLGetType(stdlink); 
        
        switch (MLGetType(stdlink)) {
                case MLTKINT:   /* Mma sent integer number*/
                        MLGetInteger(stdlink, &mlint);
                        rf2=mlint;
                        break;
                case MLTKREAL:  /* Mma sent real number */
                        MLGetReal(stdlink, &rf2);
                        break;
                default:  /* Return an error */
                        MLPutSymbol(stdlink, "errnan");
                        return;
        }

        /* Then the double precision rr2, which mma can send variously */

        dttype=MLGetType(stdlink); 
        
        switch (MLGetType(stdlink)) {
                case MLTKINT:   /* Mma sent integer number*/
                        MLGetInteger(stdlink, &mlint);
                        rr2=mlint;
                        break;
                case MLTKREAL:  /* Mma sent real number */
                        MLGetReal(stdlink, &rr2);
                        break;
                default:  /* Return an error */
                        MLPutSymbol(stdlink, "errnan");
                        return;
        }

/* calling FORTRAN  function  */
        MSBARF(&nf, &j, &rf2, &rr2, bigc0, bigc1);

/* returning result to Mathematica via MathLink 
 */
        MLPutFunction(stdlink, "List", 2);
                MLPutFunction(stdlink, "List", 2);
                        mlputcomplex(stdlink, bigc0[0]);
                        mlputcomplex(stdlink, bigc0[1]);
                MLPutFunction(stdlink, "List", 2);
                        mlputcomplex(stdlink, bigc1[0]);
                        mlputcomplex(stdlink, bigc1[1]);

        return;
}
