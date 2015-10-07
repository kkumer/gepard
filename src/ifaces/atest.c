
#include "header.h"


main()
{  
long int evoli=1, evolj=1;

printf("--- Point B (CSBAR) ---\n");
printf("Expect NLO =  -25716.505846 + 1134612.30 I\n");


readpar_();
strcpy(parchr_, "CSBARHARD  DVCS  SINGLET");

parflt_.q02 = 2.5;
kinematics_.q2 = 2.5;
kinematics_.xi = 0.00001;
nqs_.nqs = 1;
qs_.qs[0] = kinematics_.q2;

init_();
kinematics_.del2 = 0.0;
evolc_(&evoli, &evolj);

cfff_();

printf("Result: %f + %f i\n", cff_.cff[1].dr, cff_.cff[1].di);
};

