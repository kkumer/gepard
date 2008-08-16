/*
 *  Flush output from fotran, so that one can easier
 *  observe output in terminal, while saving it to file
 *   */

#include <stdio.h>
/* #include <gsl/gsl_sf_gamma.h> */
/* #include <math.h>
 * #include <complex.h> */

void flushout_() 
{
          fflush(NULL);
} 

/*
 *void clg_ (double* rez, double* imz, double* relg, double* imlg){
 *
 *      int K;
 *      int NMAX = 2;
 *      int NN = 2;
 *      static double EE = 2.71828182845904524;
 *      static double INVEE = 0.3678794411714423215955;
 *      static double G = 1.49614999;
 *      static double COEF[] = { 0.5613831815, 0.6055625385 };
 *      double complex TMP, X, YY, S, CLNGAMMA;
 *      double complex aux1;
 *
 *      X = *rez + *imz*I;
 *      YY = X - 1;
 *      TMP = X - 0.5;
 *      aux1 =  (TMP + G) * INVEE;
 *      TMP = cpow(aux1, TMP);
 *      S = COEF[0];
 *      for (K = 1; K != NN; ++K){
 *        YY = YY + 1;
 *        S = S + COEF[K] / YY;
 *      }
 *
 *      CLNGAMMA = clog(TMP * S);
 *
 *      *relg = creal(CLNGAMMA);
 *      *imlg = cimag(CLNGAMMA);
 *
 *}
 */

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
