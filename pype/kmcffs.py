#!/usr/bin/python

import sys,shelve, Model, Data

db = shelve.open('theories.db',  flag='r')
th = db['KM15']

cffs = ['ImH', 'ReH', 'ImE', 'ReE', 'ImHt', 'ReHt', 'ImEt', 'ReEt']


pt = Data.DummyPoint()

def getcffs(xB, t, Q2):
    pt.Q2 = Q2
    pt.t = t
    pt.xB = xB
    pt.xi = xB/(2.-xB)

    row = [pt.xB, pt.t, pt.Q2]
    for cff in cffs:
        th.m.g.newcall = 1
        row.append(getattr(th.m, cff)(pt))

    return row


if __name__ == '__main__':
    usage = """

  kmcffs.py  ModelID  xB  t  Q2

returns 8 DVCS CFFs Output is:

xB  t  Q2  ImH  ReH  ImE  ReE  ImHt   ReHt  ImEt  ReEt


ModelID is one of  
   1 KM09a - arXiv:0904.0458 fit without Hall A data,
   2 KM09b - arXiv:0904.0458 fit with Hall A harmonics ratio,
   3 KM10  - arXiv:1105.0899 fit with Hall A harmonics
   4 KM10a - arXiv:1105.0899 fit without Hall A data
   5 KM10b - arXiv:1105.0899 fit with Hall A harmonics ratio
   6 KMM12 - arXiv:1301.1230 fit with Hall A harmonics and polarized target
   7 KM15  - arXiv:1512.09014 fit now includes 2015 CLAS and Hall A data

Example:
    ./kmcffs.py  KM15  0.1  -0.3  4.0

(Output:)

0.100000 -0.300000 4.000000 9.049867 -1.686808 0.000000 1.900110 2.225597 0.938001 0.000000 328.094224
"""

    try:
        modelname = sys.argv[1]
        numargs = [float(k) for k in sys.argv[2:]]
    except:
        sys.stdout.write(usage)
        sys.exit(1)
    if len(numargs) != 3:
        sys.stdout.write(usage)
        sys.exit(1)
    else:
        th = db[modelname]
        lst = getcffs(*numargs)
        sys.stdout.write((11* "{:.6f} "+"\n").format(*lst))
    sys.stdout.flush()
