0   <--   1=INTERACTIVE   0=BATCH 
SET TITLE
       Fitting to ansatz: FIT (3 PWs)
PARAMETERS
11 'NS        '   0.15204   .1   
12 'AL0S      '   1.1575    .1   
13 'ALPS      '   0.15      .1  
14 'M02S      '   0.4416    .1 
15 'DELM2S    '   0.00      .1
16 'PS        '   2.0       .1
17 'SECS      '  -0.449165  .1
18 'KAPS      '   0.0       .1
19 'SKEWS     '   0.0       .1
21 'NG        '   0.5       .1 
22 'AL0G      '   1.2473    .1
23 'ALPG      '   0.15      .1 
24 'M02G      '   0.7       .1
25 'DELM2G    '   0.00      .1
26 'PG        '   2.0       .1
27 'SECG      '  -2.370178  .1 
28 'KAPG      '   0.0       .1 
29 'SKEWG     '   0.0       .1 
32 'THIS      '   0.091884  .1 
42 'THIG      '   0.844465  .1
                                
fix  11 12 13    15 16 17 18 19  32
fix  21 22 23 24 25 26    28 29 
                                
cali 3                          
                                