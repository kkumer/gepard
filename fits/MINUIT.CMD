SET TITLE
  ---->  Fitting to H1/ZEUS DVCS and DIS data  <-----
PARAMETERS
1  'NG        '   0.4d0     .1d0
2  'NSEA      '   0.05d0     .1d0
3  'MG        '   1.0d0     .1d0 
4  'MSEA      '   1.0d0     .1d0
5  'ALPHA0G   '   1.2d0     .1d0
6  'ALPHA0SEA '   0.4d0     .1d0
7  'Q20       '   1.0d0     .1d0

fixed 3 4 5 6 7
migrad
stop
