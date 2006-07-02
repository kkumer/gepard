/*  
    file: astrong.c                           
    C wrapper which calls FORTRAN subroutine on behalf of MathLink

    running-coupling constant of QCD  / (2 pi)
    initial value as(2.5)=0.05 is hard-wired here 
*/

#include "mathlink.h"
#include "gepard.h"

#pragma aux AS2PF "^"  /* For Watcom compiler, GCC ignores it */
void AS2PF(double *, double *, double *, double *, int *, int *);

void as2p(void)
{
        int mlint, nf, ord, dttype;
        double mlreal, mlr, mli;
        double mu2;
        double as, mu20, as0;
        long int nargs; const char *fname;

        /* Receiving parameters via MathLink */
        /* double precision mu2, which mma can send variously */

        dttype=MLGetType(stdlink); 
        
        switch (MLGetType(stdlink)) {
                case MLTKINT:   /* Mma sent integer number*/
                        MLGetInteger(stdlink, &mlint);
                        mu2=mlint;
                        break;
                case MLTKREAL:  /* Mma sent real number */
                        MLGetReal(stdlink, &mu2);
                        break;
                default:  /* Return an error */
                        MLPutSymbol(stdlink, "errnan");
                        return;
        }

        /* int nf, which mma can send variously */

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

        /* int ord, which mma can send variously */

        dttype=MLGetType(stdlink); 
        
        switch (MLGetType(stdlink)) {
                case MLTKINT:   /* Mma sent integer number*/
                        MLGetInteger(stdlink, &mlint);
                        ord=mlint;
                        break;
                case MLTKREAL:  /* Mma sent real number */
                        MLGetReal(stdlink, &mlreal);
                        ord = (int) mlreal;
                        break;
                default:  /* Return an error */
                        MLPutSymbol(stdlink, "errnan");
                        return;
        }


        mu20=2.5;
        as0=0.05;
        AS2PF(&as, &mu2, &as0, &mu20, &nf, &ord); 
        MLPutReal(stdlink, as);
}

