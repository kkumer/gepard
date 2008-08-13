/*
 *  Flush output from fotran, so that one can easier
 *  observe output in terminal, while saving it to file
 *   */

#include <stdio.h>
/* #include <gsl/gsl_sf_gamma.h> */

void flushout_() 
{
          fflush(NULL);
} 


/*
*void clg_ (double* rez, double* imz, double* relg, double* imlg)
*{
*  gsl_sf_result lnr;
*  gsl_sf_result arg;
*  gsl_sf_lngamma_complex_e (*rez, *imz, &lnr, &arg);
*  *relg = lnr.val; 
*  *imlg = arg.val;
*  return;
*}
*/
