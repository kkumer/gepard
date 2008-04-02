
/*   "header" file with definitions */

/*
 *   Parameters
 *   ----------
 */

/* Array dimensions */
#define NPARMAX   70
#define QINDMAX   50
#define NBNDMAX    1
#define NPTSMAX  768
#define MTINDMAX 100

struct dblcomplex {
        double dr;
        double di;
};


/*
 *   Common blocks
 *   -------------
 */

/*   1. Constant parameters loaded from 'GEPARD.INI' or fixed by program */

      /* CHARACTER SCHEME*5, ANSATZ*6, PROCESS*6, FFTYPE*10 */
      /* COMMON / PARINT /  SPEED, ACC, P, NF, CZERO */
      /* DOUBLE PRECISION MU02, ASP(0:2) */
      /* COMMON / ASTRONG/  MU02, ASP */
      /* DOUBLE PRECISION Q02, RF2, RR2 */
      /* COMMON / PARFLT /  Q02, RF2, RR2 */
      /* INTEGER SPEED, ACC, P, NF, CZERO */
      /* COMMON / PARCHR /  SCHEME, ANSATZ, PROCESS, FFTYPE */
      /* COMMON / FILENAMES /  DATFILE, OUTFILE, CMDFILE */

struct{
        int speed, acc, p, nf, czero;
} parint_;
struct{
        double mu02, asp[3];
} astrong_;
struct{
        double q02, rf2, rr2;
} parflt_;
char parchr_[27];
char filenames_[60];


/*   2. Parameters from 'MINUIT.CMD' (candidates for fitting parameters) */

      /* INTEGER NPARMAX
       * PARAMETER (NPARMAX = 70)
       * DOUBLE PRECISION PAR(NPARMAX)
       * COMMON / PAR    /  PAR
       */

struct{
        double par[NPARMAX];
} par_;

/*   3. Kinematics */

      /* INTEGER QINDMAX, MTINDMAX
       *PARAMETER (QINDMAX = 50, MTINDMAX = 100)
       *INTEGER NQS, NQSDIS, MTIND, NMTS, NMTSEXP
       *DOUBLE PRECISION XI, DEL2, Q2
       *DOUBLE PRECISION QS(QINDMAX), QSDIS(QINDMAX)
       *DOUBLE PRECISION MTS(MTINDMAX), MTWG(MTINDMAX)
       *COMMON / KINEMATICS /  XI, DEL2, Q2
       *COMMON / NQS        /  NQS, NQSDIS
       *COMMON / QS         /  QS, QSDIS
       *COMMON / MT         /  MTIND, NMTS, NMTSEXP
       *COMMON / MTS        /  MTS, MTWG
       */

struct{
        double xi, del2, q2;
} kinematics_;

struct{
        int nqs, nqsdis;
} nqs_;

struct{
        double qs[QINDMAX], qsdis[QINDMAX];
} qs_;

struct{
        int mtind, nmts, nmtsexp;
} mt_;

struct{
        double mts[MTINDMAX], mtwg[MTINDMAX];
} mts_;

/*   4. Other */

/*
 **     - Mellin-Barnes integration contour points
 *      INTEGER NPTS, NPTSMAX
 *      PARAMETER (NPTSMAX = 768)
 *      DOUBLE PRECISION Y(NPTSMAX), WG(NPTSMAX)
 *      DOUBLE COMPLEX N(NPTSMAX)
 *      COMMON / CONTOUR  /  NPTS
 *      COMMON / POINTS   /  Y, WG
 *      COMMON / NPOINTS  /  N
 */

struct{
        int npts;
} contour_;

struct{
        struct dblcomplex n[NPTSMAX];
} npoints_;

/*
 **     - Values on the whole contour 
 *      DOUBLE COMPLEX HGRID(0:MTINDMAX, NPTSMAX, 2)
 *      COMMON / HGRID    /  HGRID
 */

struct{
        struct dblcomplex hgrid[2][NPTSMAX][MTINDMAX+1];
} hgrid_;

struct{
        struct dblcomplex mbgpd[2][NPTSMAX];
} mbgpd_;


/*     - Final observables */
      /* DOUBLE COMPLEX CFF(0:2)
       * COMMON / CFF      /  CFF
       * DOUBLE PRECISION F2(0:2)
       * COMMON / F2       /  F2
       */

struct{
        struct dblcomplex cff[3], cffe[3];
} cff_;

struct{
        double f2[3];
} f2_;

      /*
       *COMMON / CHISQBLK /  CHISQ
       *COMMON / NDATAPTSBLK /  NDATAPTS
       *COMMON / CHINBLK /  CHIN
       */

struct{
        double chisq[11];
} chisqblk_;

struct{
        int ndatapts[11];
} ndataptsblk_;

struct{
        int chin;
} chinblk_;

/*
 *   Prototypes
 *   -------------
 */

void cfff_();
void readpar_();
void init_();
void initgpd_();
void getmbgpdmma_();
void evolc_(long int *i, long int *j);
void mpar_(long int *id, long int *siz, char* pnam, double *vstrt, double *stp, double *lo, double *hi, int ch_len);
void mcom_(long int *siz, char* cmd, long int *ierflg, int cmd_len);
