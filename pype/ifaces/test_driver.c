
#include <stdio.h>
#include "ifaces.h"

int main(void) {
    double phi, res;
    int k;

    for (k=1; k<10; k++) {
        phi = 2. - k/10.;
        res = xs(1, -1, 1, 27.6, 0., 0.1, 3.2, -0.3, phi);
        printf("HERMES: phi = %8.3f, xs = %8.3e\n", phi, res);
        res = xs(1, -1, 1, 5.0, 250., 0.1, 3.2, -0.3, phi);
        printf("EIC: phi = %8.3f, xs = %8.3e\n", phi, res);
    }
    return 0;
}

