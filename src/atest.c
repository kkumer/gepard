
#include "header.h"


main()
{  
long int evoli=1, evolj=1;

printf("--- Point B (CSBAR) ---\n");
printf("Expect NLO =  -25716.505846 + 1134612.30 I\n");

strcpy(parchr_, "12345678901DVCS  SINGLET");

readpar_();

par_.par[49] = 0.0;
par_.par[50] = 0.0;

kinematics_.q2 = 2.5;
kinematics_.xi = 0.00001;
nqs_.nqs = 1;
qs_.qs[0] = kinematics_.q2;

init_();
initgpd_();
mt_.mtind = 0;
evolc_(&evoli, &evolj);

cfff_();

printf("Result: %f + %f i\n", cff_.cff[1].dr, cff_.cff[1].di);
};

