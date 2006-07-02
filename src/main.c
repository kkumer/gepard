/* Main C program for MathLink-ing */


#include "mathlink.h"
#include "gepard.h"

#pragma aux BETAF "^"  /* For Watcom compiler, GCC ignores it */
void BETAF(void);

int main(int argc, char *argv[]) 
{
           BETAF();  /* beta funct. initialization */
           return MLMain(argc, argv);
}
