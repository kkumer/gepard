
/* functions that facilitate communication with MathLink */

/* FIXME: stdlink should be declared? */

#include "mathlink.h"
#include "gepard.h"

void mlputcomplex(MLINK stdlink, Doublecomplex z)
{
            MLPutFunction(stdlink, "Complex", 2);
            MLPutReal(stdlink, z.dr);
            MLPutReal(stdlink, z.di);
}

void mlputcomplex2by2matrix(MLINK stdlink, Doublecomplex m[2][2])
{
        int i, j;

        MLPutFunction(stdlink, "List", 2);
        for (i=0; i<2; i++) {
            MLPutFunction(stdlink, "List", 2);
            for (j=0; j<2; j++) 
                mlputcomplex(stdlink, m[j][i]);
        }
}

