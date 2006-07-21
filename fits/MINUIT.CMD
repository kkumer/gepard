0   <--   1=INTERACTIVE   0=BATCH
SET TITLE
       Fitting to H1/ZEUS DVCS and DIS data
PARAMETERS
1  'NG        '   0.6d0     .1d0
2  'NSEA      '   0.1d0     .1d0
3  'MG        '   1.0d0     .1d0 
4  'MSEA      '   1.0d0     .1d0
5  'ALPHA0G   '   -0.08d0     .1d0
6  'ALPHA0SEA '   1.24d0     .1d0
7  'Q20       '   1.0d0     .1d0

fixed 3 4 7 
migrad
cali 3
