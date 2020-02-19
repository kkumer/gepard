## Synopsis

**cffs4mc** - serving CFF values of KM models


## Preparation

**Prerequisites**: C and Fortran compilers. Python and Python libraries and headers (these are in package python-dev on Debian linux and in package python on Arch linux).

First make a `pygepard.so` extension library.
```sh
% make pygepard
```

Then go into pype directory and compile C-wrappers
```sh
% cd pype
% gcc -c $(python3-config --cflags) kmcffs.c && gcc kmcffs.o $(python3-config --ldflags)  -lpython3.8  -o kmcffs
% gcc -c $(python3-config --cflags) xbloop.c && gcc xbloop.o $(python3-config --ldflags)  -lpython3.8  -o xbloop
```
First one is for specifying kinematics on command line, and second one is example which loops over range of xB values.

## Using

```sh
% ./kmcffs 
Usage: kmcffs xB t Q2
% ./kmcffs 0.1 -0.3 4.0
0.100000 -0.300000 4.000000 9.049867 -1.686808 0.000000 1.900110 2.225597 0.938001 0.000000 328.094224
```

```sh
% ./xbloop 
   xB        t         Q2      ImH      ReH      ImE      ReE      ImHt     ReHt    ImEt      ReEt
0.010000 -0.300000 4.000000 116.702178 17.066079 0.000000 1.900110 3.617266 1.087262 0.000000 3436.355292 
0.032821 -0.300000 4.000000 30.989010 2.017216 0.000000 1.900110 2.878577 0.964549 0.000000 1035.007765 
0.055641 -0.300000 4.000000 17.216320 -0.347707 0.000000 1.900110 2.570554 0.940686 0.000000 603.429178 
0.078462 -0.300000 4.000000 11.781983 -1.252491 0.000000 1.900110 2.369804 0.936111 0.000000 422.899572 
0.101282 -0.300000 4.000000 8.926528 -1.705772 0.000000 1.900110 2.217930 0.938213 0.000000 323.722549 
0.124103 -0.300000 4.000000 7.189143 -1.962653 0.000000 1.900110 2.093960 0.942998 0.000000 261.019720 
0.146923 -0.300000 4.000000 6.032434 -2.116574 0.000000 1.900110 1.988022 0.948824 0.000000 217.795256 
0.169744 -0.300000 4.000000 5.213580 -2.209886 0.000000 1.900110 1.894697 0.954911 0.000000 186.193081 
0.192564 -0.300000 4.000000 4.607360 -2.264501 0.000000 1.900110 1.810697 0.960860 0.000000 162.081168 
0.215385 -0.300000 4.000000 4.142806 -2.292852 0.000000 1.900110 1.733880 0.966456 0.000000 143.078684 
[...]
```

Second output should be consistent with what you get by going to
[http://calculon.phy.hr/gpd/server/CFF-grid.html](GPD server) and just switching
from `log` to `lin`.




## License

Copyright kkumer@phy.hr
All rights reserved.
