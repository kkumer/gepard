#include "mathlink.h"

/* double complex type that can be simply passed to FORTRAN and returned */
struct dblcomplex {
        double dr; 
        double di;
};
typedef struct dblcomplex Doublecomplex;

/* Prototypes for external functions */

void mlputcomplex(MLINK stdlink, struct dblcomplex z);
void mlputcomplex2by2matrix(MLINK stdlink, struct dblcomplex m[2][2]);



/*
 *	  Linux/Windows translations for calling Fortran from C
 */

#if defined(__linux__) && defined(__i386__) && defined(__GNUC__)
  #define    AS2PF       as2pf_
  #define    CDVCSF      cdvcsf_
  #define    COMMONF     commonf_
  #define    ERFUNCF     erfuncf_
  #define    EVOLF       evolf_
  #define    LAMBDAF     lambdaf_
  #define    BETAF       betaf_
  #define    PROJECTORSF projectorsf_
  #define    RNNLOF      rnnlof_
  #define    MSBARF      msbarf_
#else  /* Windows Watcom compiler, the only one I tried */
  #error "Unknown system/compiler. Write gepard.h stanza."
#endif

