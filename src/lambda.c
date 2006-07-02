/*  
    file: lambda.c                           
    C wrapper which calls FORTRAN subroutine on behalf of MathLink

    eigenvalues of LO anomalous dimensions matrix
*/

#include "mathlink.h"
#include "gepard.h"

#pragma aux LAMBDAF "^"  /* For Watcom compiler, GCC ignores it */
void LAMBDAF(int *nf, Doublecomplex *n, Doublecomplex lam[2], 
                Doublecomplex gam0[2][2]);

void lambda(void)
{
        int mlint, nf, dttype;
        double mlreal, mlr, mli;
        Doublecomplex n;
        Doublecomplex lam[2]; Doublecomplex gam0[2][2];
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
        LAMBDAF(&nf, &n, lam, gam0);

/* returning result to Mathematica via MathLink  */
        MLPutFunction(stdlink, "List", 2);
        mlputcomplex(stdlink, lam[0]);
        mlputcomplex(stdlink, lam[1]);
        return;
};

