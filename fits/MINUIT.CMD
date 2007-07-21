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
21 'NG        '   0.5       .1  
22 'AL0G      '   1.20      .1   
23 'ALPG      '   0.15      .1  
24 'M02G      '   1.1       .1
25 'DELM2G    '   0.00      .1  
26 'PG        '   2.0       .1  
31 'NU        '   0.0       .1  
32 'AL0U      '   0.5       .1   
33 'ALPU      '   1.0       .1  
34 'M02U      '   3.53      .1  
35 'DELM2U    '   0.88      .1  
36 'PU        '   1.0       .1  
41 'ND        '   0.0       .1  
42 'AL0D      '   0.5       .1   
43 'ALPD      '   1.0       .1  
44 'M02D      '   3.53      .1  
45 'DELM2D    '   0.88      .1  
46 'PD        '   1.0       .1  
50 'SECQ      '   0.0       .1  
51 'SECG      '   0.0       .1  


fix  13 15 16
fix  23 25 26
fix  31 32 33 34 35 36
fix  41 42 43 44 45 46
fix  50 51

migrad
cali 3
