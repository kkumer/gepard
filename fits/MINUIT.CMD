0   <--   1=INTERACTIVE   0=BATCH  (Must be 0 on Windows!)
SET TITLE
       Fitting to ansatz: FIT
PARAMETERS
11 'NS        '   0.18      .1  
12 'AL0S      '   1.1       .1   
13 'ALPS      '   0.15      .1  
14 'M02S      '   1.4       .1
15 'DELM2S    '   0.00      .1  
16 'PS        '   3.0       .1  
19 'SKEWS     '   1.0       .1  
21 'NG        '   0.5       .1  
22 'AL0G      '   1.20      .1   
23 'ALPG      '   0.15      .1  
24 'M02G      '   1.1       .1
25 'DELM2G    '   0.00      .1  
26 'PG        '   2.0       .1  
29 'SKEWG     '   1.0       .1  

fix  13 14 15 16 
fix  23 24 25 26

migrad

cali 3
