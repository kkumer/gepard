
#include "mathlink.h"
#include "ifaces.h"

void XSunp(long int id, long int Q, long int lam, double Ee, double Ep, double xB, double Q2, double t, double phi)
{
    double res;

    res = xs(id, Q, lam, Ee, Ep, xB, Q2, t, phi);

    MLPutReal(stdlink, res);

    return;
};


int main(int argc, char *argv[]) {
           return MLMain(argc, argv);
}
