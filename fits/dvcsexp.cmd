0   <--   1=INTERACTIVE   0=BATCH 
SET TITLE
       Fitting to ansatz: FIT
PARAMETERS
11 'NS        '   0.1520386 .1  
12 'AL0S      '   1.157507  .1   
13 'ALPS      '   0.00      .1  
14 'M02S      '   0.2       .1
15 'DELM2S    '   0.00      .1  
16 'PS        '   3.0       .1  
17 'SECS      '   0.0       .1  
18 'KAPS      '   0.0       .1  
19 'SKEWS     '   0.0       .1  
21 'NG        '   0.5       .1  
22 'AL0G      '   1.2473157 .1   
23 'ALPG      '   0.00      .1  
24 'M02G      '   0.2159827 .1
25 'DELM2G    '   0.00      .1  
26 'PG        '   2.0       .1  
27 'SECG      '   0.0       .1  
28 'KAPG      '   0.0       .1  
29 'SKEWG     '   0.0       .1  
32 'THIS      '   0.0       .1  
42 'THIG      '   0.0       .1  
111 'NSE      '   0.0       .1


fix  11 12 13    15 16 18 19 32
fix  21 22 23 24 25 26 28 29 42
fix 111

migrad

cali 3
