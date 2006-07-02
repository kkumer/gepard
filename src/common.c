/****h* gepard/common.c
 *  NAME
 *    common.c
 *  DESCRIPTION
 *********/
/*  
    file: common.c                           
    C wrapper which calls FORTRAN subroutine on behalf of MathLink

    initialization of common block Harmonic
*/

#include "mathlink.h"
#include "gepard.h"

#pragma aux COMMONF "^"  /* For Watcom compiler, GCC ignores it */
void COMMONF(Doublecomplex *n, int *nf);

/****f* common.c/common
 *  NAME
 *    common
 *  SYNOPSIS
 *    void common(void)
 *****/
void common(void)
{
        int nf, mlint, dttype;
        double mlreal, mlr, mli;
        Doublecomplex n;
        long int nargs; const char *fname;


/* getting input via MathLink.  */

        /* Complex n, which mma can send variously */

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

        /* Then int nf, which mma can send variously */

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

                       

/* calling FORTRAN  function  */
        COMMONF(&n, &nf);
        mlint = 1;
        MLPutInteger(stdlink, mlint);

        return;
};

