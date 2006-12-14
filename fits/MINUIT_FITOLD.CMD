0   <--   1=INTERACTIVE   0=BATCH  (Must be 0 on Windows!)
SET TITLE
       Fitting using ansatz: FITOLD
PARAMETERS
1  'Q02       '   7.0       .1  
2  'AS0       '   0.0536    .1  
3  'MU02      '   2.5       .1  
11 'NS        '   0.2       .1  
12 'AL0S      '   1.1       .1   
13 'ALPS      '   0.25      .1  
14 'MS        '   1.0       .1  
21 'NG        '   0.5       .1  
22 'AL0G      '   1.0       .1   
23 'ALPG      '   0.25      .1  
24 'MG        '   1.0       .1  


fix  1 2 3 13 23

migrad
cali 3
