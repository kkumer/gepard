
parstemplate = ('NS', 'alS', 'alpS', 'MS', 'rS', 'bS', 'Nv', 'alv', 'alpv', 'Mv', 'rv', 'bv', 
        'C', 'MC', 'tNv', 'tMv', 'trv', 'tbv')

# DM's GLO fit to BCAcos1 and ALUIsin1 (HERMES) and ALUa (CLAS)

DMGLOHt={'NS': 1.5,
        'alS': 1.13,  
        'alpS': 0.15,  
        'MS': 0.707107,   
        'rS': 1.0,       
        'bS': 3.07387,  
        'Nv': 1.35,   
        'alv': 0.43,    
        'alpv': 0.85, 
        'Mv': 0.742178, 
        'rv': 0.729705, 
        'bv': 1.00418,  
        'C': 2.1928,
       'MC': 1.71754, 
       'tNv': 0.6, 
       'tMv': 1.71063, 
       'trv': 0.700422,   
       'tbv': 1.07523   
       }

DMGLOHttpl = tuple([DMGLOHt[key] for key in parstemplate])

# Same, but without \tilde{H}

DMGLO = {'C': 1.1207235939592248, 
 'MC': 1.2216488253146187,        
 'MS': 0.7071067811865476,        
 'Mv': 0.6829863538482371,        
 'NS': 1.5,                           
 'Nv': 1.3500000000000001,            
 'alS': 1.1299999999999999,           
 'alpS': 0.14999999999999999,         
 'alpv': 0.84999999999999998,         
 'alv': 0.42999999999999999,          
 'bS': 2.254514866291967,             
 'bv': 0.5,                           
 'rS': 1.0,                           
 'rv': 0.6842807978805212,            
 'tMv': 2.69,           
 'tNv': 0.0,                          
 'tbv': 3.2560699999999998,           
 'trv': 5.9792300000000003}           

DMGLOtpl = tuple([DMGLO[key] for key in parstemplate])

# As DMGLO, but with scan-minimization of variable params

DMGLOB =  {'C': 1.4937499999999984,
 'MC': 1.092,
 'MS': 0.707107,
 'Mv': 1.0,
 'NS': 1.5,
 'Nv': 1.3500000000000001,
 'alS': 1.1299999999999999,
 'alpS': 0.14999999999999999,
 'alpv': 0.84999999999999998,
 'alv': 0.42999999999999999,
 'bS': 3.1900000000000004,
 'bv': 0.5,
 'rS': 1.0,
 'rv': 0.71250000000000002,
 'tMv': 2.7,
 'tNv': 0.0,
 'tbv': 3.2560699999999998,
 'trv': 5.9792300000000003}

DMGLOBtpl = tuple([DMGLOB[key] for key in parstemplate])


# As DMGLOB, but with global-minimization with MINUIT afterwards

DMGLOKK = {'C': 1.4056450660676154,
 'MC': 1.05,
 'MS': 0.707107,
 'Mv': 0.95,
 'NS': 1.5,
 'Nv': 1.3500000000000001,
 'alS': 1.1299999999999999,
 'alpS': 0.14999999999999999,
 'alpv': 0.84999999999999998,
 'alv': 0.42999999999999999,
 'bS': 3.2999814723772429,
 'bv': 0.5,
 'rS': 1.0,
 'rv': 0.73469278205761201,
 'tMv': 2.7,
 'tNv': 0.0,
 'tbv': 3.2560699999999998,
 'trv': 5.9792300000000003}


DMGLOKKtpl = tuple([DMGLOKK[key] for key in parstemplate])

# KK's GLO fit to BCAcos1 and ALUIsin1 (HERMES) and ALUa (CLAS), without \tilde{H}
# P(chi-square, d.o.f) = P(28.31, 30) = 0.5539

KKGLO = {'C': 1.7463470186048999,
 'MC': 0.707,
 'MS': 0.707107,
 'Mv': 1.41,
 'NS': 1.5,
 'Nv': 1.3500000000000001,
 'alS': 1.1299999999999999,
 'alpS': 0.14999999999999999,
 'alpv': 0.84999999999999998,
 'alv': 0.42999999999999999,
 'bS': 0.50000002008444855,
 'bv': 0.50000013966239731,
 'rS': 1.0,
 'rv': 0.39057585298799602,
 'tMv': 2.7,
 'tNv': 0.0,
 'tbv': 3.2560699999999998,
 'trv': 5.9792300000000003}

KKGLOtpl = tuple([KKGLO[key] for key in parstemplate])

# DM's GLO1 fit to BCAcos1 and ALUIsin1 (HERMES) and ALUa (CLAS) and b1ovb0 (HALL-A)

DMGLO1={'NS': 1.5,
        'alS': 1.13,  
        'alpS': 0.15,  
        'MS': 0.707107,   
        'rS': 1.0,       
        'bS': 2.00203,  
        'Nv': 1.35,   
        'alv': 0.43,    
        'alpv': 0.85, 
        'Mv': 1.01097, 
        'rv': 0.496383, 
        'bv': 2.15682,  
        'C': 6.90484,
       'MC': 1.33924, 
       'tNv': 0.6, 
       'tMv': 2.69667, 
       'trv': 5.97923,   
       'tbv': 3.25607   
       }

DMGLO1tpl = tuple([DMGLO1[key] for key in parstemplate])

DMGLO1B={'NS': 1.5,
        'alS': 1.13,  
        'alpS': 0.15,  
        'MS': 0.707107,   
        'rS': 1.0,       
        'bS': 2.15,  
        'Nv': 1.35,   
        'alv': 0.43,    
        'alpv': 0.85, 
        'Mv': 0.898, 
        'rv': 0.496383, 
        'bv': 2.15682,  
        'C': 6.81,
       'MC': 1.33924, 
       'tNv': 0.6, 
       'tMv': 2.736, 
       'trv': 5.97923,   
       'tbv': 3.39   
       }

# pype fit to same data, with same parameters released
# ncalls =  3016
# P(chi-square, d.o.f) = P(14.90, 31) = 0.9935

DMGLO1KK = {'C': 6.3735252628367105,
 'MC': 1.42,
 'MS': 0.707107,
 'Mv': 1.53,
 'NS': 1.5,
 'Nv': 1.3500000000000001,
 'alS': 1.1299999999999999,
 'alpS': 0.14999999999999999,
 'alpv': 0.84999999999999998,
 'alv': 0.42999999999999999,
 'bS': 1.5000005888078514,
 'bv': 1.200000189885863,
 'rS': 1.0,
 'rv': 0.30903306714870188,
 'tMv': 3.0,
 'tNv': 0.59999999999999998,
 'tbv': 3.6300100917776876,
 'trv': 6.1340651857666719}



DMGLO1KKtpl = tuple([DMGLO1KK[key] for key in parstemplate])

# pype fit, adding data[22] i.e. phi-dependent HALL-A  at t=-0.23 GeV^2
# P(chi-square, d.o.f) = P(89.92, 73) = 0.0871


pfit1 = {
  'C': 7.076250109725728,
 'MC': 0.964,
 'MS': 0.707107,
 'Mv': 0.45,
 'NS': 1.5,
 'Nv': 1.3500000000000001,
 'alS': 1.1299999999999999,
 'alpS': 0.14999999999999999,
 'alpv': 0.84999999999999998,
 'alv': 0.42999999999999999,
 'bS': 0.44927870053766178,
 'bv': 1.0492124373622569,
 'rS': 1.0,
 'rv': 0.42506259798769952,
 'tMv': 2.24,
 'tNv': 0.59999999999999998,
 'tbv': 4.0822813556178339,
 'trv': 8.611140873616101}



pfit1tpl = tuple([pfit1[key] for key in parstemplate])

pfit2 = {
  'C': 6.8422440523143369,
 'MC': 1.0103180741897759,
 'MS': 0.70710700000000004,
 'Mv': 0.20409699243778689,
 'NS': 1.5,
 'Nv': 1.3500000000000001,
 'alS': 1.1299999999999999,
 'alpS': 0.14999999999999999,
 'alpv': 0.84999999999999998,
 'alv': 0.42999999999999999,
 'bS': 1.0077537498166014,
 'bv': 1.6074343029487501,
 'rS': 1.0,
 'rv': 0.69999988135649893,
 'tMv': 8.999998478702004,
 'tNv': 0.59999999999999998,
 'tbv': 2.8835965728359838,
 'trv': 6.9999999127837178}


# after fixing 1+eps
#{0.9479,0.8,0.4481,3.0879,-0.2924,0.5002,0.,0.8,3.}
#{rv,Mv,bv,bs,d,Mc,rtv,Mtv,btv}
# chisq/d.o.f  = 32/31  (npt=36)

DMepsGLO = {
  'C': 0.24,
 'MC': 0.5002,
 'MS': 0.7071067811865476,        
 'Mv': 0.8,
 'NS': 1.5,                           
 'Nv': 1.3500000000000001,            
 'alS': 1.1299999999999999,           
 'alpS': 0.14999999999999999,         
 'alpv': 0.84999999999999998,         
 'alv': 0.42999999999999999,          
 'bS': 3.0879,
 'bv': 0.4481,
 'rS': 1.0,                           
 'rv': 0.9479,
 'tMv': 0.8,
 'tNv': 0.0,                          
 'tbv': 3.0,
 'trv': 0.0}           

#{1.1064,0.8,2.398,4.5845,-6.042,1.5238,4.8966,0.8,1.5}
#{rv,Mv,bv,bs,d,Mc,rtv,Mtv,btv}
# chisq/??? = 0.97
DMepsGLO1 = {
  'C': 6.042,
 'MC': 1.5238,
 'MS': 0.7071067811865476,        
 'Mv': 0.8,
 'NS': 1.5,                           
 'Nv': 1.3500000000000001,            
 'alS': 1.1299999999999999,           
 'alpS': 0.14999999999999999,         
 'alpv': 0.84999999999999998,         
 'alv': 0.42999999999999999,          
 'bS': 4.5845,
 'bv': 2.398,
 'rS': 1.0,                           
 'rv': 1.1064,
 'tMv': 0.8,
 'tNv': 0.6,                          
 'tbv': 1.5, 
 'trv': 4.8966}           


# Hybrid fit to DVCSpoints and GLOpoints
# ('M02S','SECS','SECG', 'rv', 'bv', 'C', 'MC', 'trv', 'tbv')
hy = {
 'AL0G': 1.2473167010704711,
 'AL0S': 1.1575060246398083,
 'ALPG': 0.14999999999999999,
 'ALPS': 0.14999999999999999,
 'BETG': 6.0,
 'BETS': 8.0,
 'DELM2G': 0.0,
 'DELM2S': 0.0,
 'KAPG': 0.0,
 'KAPS': 0.0,
 'M02G': 0.69999999999999996,
 'M02S': 0.46405126256515294,
 'NG': 0.5,
 'NS': 0.15203911208796006,
 'PG': 2.0,
 'PS': 2.0,
 'SECG': -0.69409006621234526,
 'SECS': -0.16659753378728048,
 'SKEWG': 0.0,
 'SKEWS': 0.0,
 'THIG': 0.0,
 'THIS': 0.0,
 'C': 9.9999999997984226,
 'MC': 0.84923046744683672,
 'MS': 0.70699999999999996,
 'Mres': 1.0,
 'Mv': 1.0,
 'Nres': 0.0,
 'Nsea': 0.0,
 'Nv': 1.3500000000000001,
 'alS': 1.1299999999999999,
 'alpS': 0.14999999999999999,
 'alpv': 0.84999999999999998,
 'alv': 0.42999999999999999,
 'bS': 2.0,
 'bres': 2.0,
 'bv': 0.40000005719394938,
 'rS': 1.0,
 'rv': 0.55167159092742724,
 'tMv': 2.7000000000000002,
 'tNv': 0.59999999999999998,
 'tbv': 3.1830461073709944,
 'trv': 7.999963926692665}

# Hybrid fit to DVCSpoints and GLO1points
# ('M02S','SECS','SECG', 'rv', 'bv', 'C', 'MC', 'trv', 'tbv')
hy1 = {
 'AL0G': 1.2473167010704711,
 'AL0S': 1.1575060246398083,
 'ALPG': 0.14999999999999999,
 'ALPS': 0.14999999999999999,
 'BETG': 6.0,
 'BETS': 8.0,
 'DELM2G': 0.0,
 'DELM2S': 0.0,
 'KAPG': 0.0,
 'KAPS': 0.0,
 'M02G': 0.69999999999999996,
 'M02S': 0.45108533906274095,
 'NG': 0.5,
 'NS': 0.15203911208796006,
 'PG': 2.0,
 'PS': 2.0,
 'SECG': -0.63062801768125887,
 'SECS': -0.17423580228386584,
 'SKEWG': 0.0,
 'SKEWS': 0.0,
 'THIG': 0.0,
 'THIS': 0.0,
 'C': 6.4994425806649616,
 'MC': 1.0577160575834483,
 'MS': 0.70699999999999996,
 'Mres': 1.0,
 'Mv': 1.0,
 'Nres': 0.0,
 'Nsea': 0.0,
 'Nv': 1.3500000000000001,
 'alS': 1.1299999999999999,
 'alpS': 0.14999999999999999,
 'alpv': 0.84999999999999998,
 'alv': 0.42999999999999999,
 'bS': 2.0,
 'bres': 2.0,
 'bv': 0.40000000011514253,
 'rS': 1.0,
 'rv': 0.63432083446847365,
 'tMv': 2.7000000000000002,
 'tNv': 0.59999999999999998,
 'tbv': 3.9488870772139104,
 'trv': 6.5439416177483896}

# Hybrid fit + 3rd PW to DVCSpoints and GLO1points  chisq = 96/129
# ('M02S','SECS','SECG', 'THIS', 'THIG', 'rv', 'bv', 'C', 'MC', 'trv', 'tbv')
hy1THI = {
 'AL0G': 1.2473167010704711,
 'AL0S': 1.1575060246398083,
 'ALPG': 0.14999999999999999,
 'ALPS': 0.14999999999999999,
 'DELM2G': 0.0,
 'DELM2S': 0.0,
 'KAPG': 0.0,
 'KAPS': 0.0,
 'M02G': 0.69999999999999996,
 'M02S': 0.49754317018981614,
 'NG': 0.5,
 'NS': 0.15203911208796006,
 'PG': 2.0,
 'PS': 2.0,
 'SECG': -2.5151319493485427,
 'SECS': -0.46005118719187721,
 'SKEWG': 0.0,
 'SKEWS': 0.0,
 'THIG': 0.89157575591751848,
 'THIS': 0.093517989519796618,
 'C': 4.868157888008362,
 'MC': 1.9147345466937296,
 'MS': 0.70699999999999996,
 'Mv': 1.0,
 'Nsea': 0.0,
 'Nv': 1.3500000000000001,
 'alS': 1.1299999999999999,
 'alpS': 0.14999999999999999,
 'alpv': 0.84999999999999998,
 'alv': 0.42999999999999999,
 'bS': 2.0,
 'bv': 0.40309038577130119,
 'rS': 1.0,
 'rv': 0.84353216408918019,
 'tMv': 2.7000000000000002,
 'tNv': 0.59999999999999998,
 'tbv': 4.9999998875540301,
 'trv': 6.9128424208463581}

# HERMES BCA uglyish. Otherwise OK.
hy1_17_THI = {
 'AL0G': 1.2473167010704711,
 'AL0S': 1.1575060246398083,
 'ALPG': 0.14999999999999999,
 'ALPS': 0.14999999999999999,
 'DELM2G': 0.0,
 'DELM2S': 0.0,
 'KAPG': 0.0,
 'KAPS': 0.0,
 'M02G': 0.69999999999999996,
 'M02S': 0.49500902065346486,
 'NG': 0.5,
 'NS': 0.15203911208796006,
 'PG': 2.0,
 'PS': 2.0,
 'SECG': -2.5869643934477677,
 'SECS': -0.47960950381037593,
 'SKEWG': 0.0,
 'SKEWS': 0.0,
 'THIG': 0.93133198156955843,
 'THIS': 0.099433539629814882,
 'C': 5.2760506533706852,
 'MC': 1.6923428827202578,
 'MS': 0.70699999999999996,
 'Mv': 1.0,
 'Nsea': 0.0,
 'Nv': 1.3500000000000001,
 'alS': 1.1299999999999999,
 'alpS': 0.14999999999999999,
 'alpv': 0.84999999999999998,
 'alv': 0.42999999999999999,
 'bS': 2.0,
 'bv': 0.78374274705991609,
 'rS': 1.0,
 'rv': 0.94219056266652901,
 'tMv': 2.7000000000000002,
 'tNv': 0.59999999999999998,
 'tbv': 4.2753844806174115,
 'trv': 7.9999635209873254}

# Hybrid fit + 3rd PW to DVCSpoints + data[48] + ALTGLO1points  chisq = 273/220
# ALTGLO1points = data[5] + data[25] + data[32] + HAD17 + HA17
# ('M02S','SECS','SECG', 'THIS', 'THIG', 'rv', 'bv', 'C', 'MC', 'trv', 'tbv')
# evolution cut at Q2 = 2 Gev^2
allTHI = {
 'AL0G': 1.2473167010704711,
 'AL0S': 1.1575060246398083,
 'ALPG': 0.14999999999999999,
 'ALPS': 0.14999999999999999,
 'DELM2G': 0.0,
 'DELM2S': 0.0,
 'KAPG': 0.0,
 'KAPS': 0.0,
 'M02G': 0.69999999999999996,
 'M02S': 0.44638594996060565,
 'NG': 0.5,
 'NS': 0.15203911208796006,
 'PG': 2.0,
 'PS': 2.0,
 'SECG': -2.3086829433781451,
 'SECS': -0.44159594968364457,
 'SKEWG': 0.0,
 'SKEWS': 0.0,
 'THIG': 0.81237709898019395,
 'THIS': 0.089236254203313298,
  'C': 8.5915831306978632,
 'MC': 0.92825390641864269,
 'MS': 0.70699999999999996,
 'Mv': 0.80000000000000004,
 'Nsea': 0.0,
 'Nv': 1.3500000000000001,
 'alS': 1.1299999999999999,
 'alpS': 0.14999999999999999,
 'alpv': 0.84999999999999998,
 'alv': 0.42999999999999999,
 'bS': 2.0,
 'bv': 1.6123630754854807,
 'rS': 1.0,
 'rv': 0.84881129590414339,
 'tMv': 1.0,
 'tNv': 0.59999999999999998,
 'tbv': 1.7612546082213929,
 'trv': 6.3096149789061649}


# DVCSpoints + data[48] + ALTGLOpoints
# cutq2 = 0.5 !!!  P(133.06, 155) = 0.8983
ALTGLO = {
'fix_THIG': False, 'KAPS': 0.0, 'fix_ALPG': True, 'fix_DELM2S': True,
'fix_SECG': False, 'KAPG': 0.0, 'fix_SECS': False, 'fix_DELM2G': True, 'PS':
2.0, 'SECG': -2.9559676899183551, 'fix_AL0S': True, 'fix_THIS': False,
'fix_KAPG': True, 'SECS': 0.46024259945701673, 'fix_M02G': True, 'NG': 0.5,
'fix_ALPS': True, 'M02S': 0.52112685191431951, 'PG': 2.0, 'ALPG':
0.14999999999999999, 'fix_M02S': False, 'NS': 0.15203911208796006, 'fix_AL0G':
True, 'fix_PG': True, 'AL0S': 1.1575060246398083, 'limit_M02S':
(0.29999999999999999, 1.5), 'fix_PS': True, 'SKEWG': 0.0, 'AL0G':
1.2473167010704711, 'SKEWS': 0.0, 'ALPS': 0.14999999999999999, 'M02G':
0.69999999999999996, 'fix_NG': True, 'THIS': -0.18555650737608792, 'fix_SKEWS':
True, 'DELM2S': 0.0, 'fix_SKEWG': True, 'fix_NS': True, 'limit_M02G':
(0.29999999999999999, 1.5), 'fix_KAPS': True, 'THIG': 0.94543535937507428,
'DELM2G': 0.0,
'tMv': 0.80000000000000004, 'rS': 1.0, 'limit_tbv': (0.40000000000000002, 5.0),
'fix_Mv': False, 'alv': 0.42999999999999999, 'fix_bv': False, 'Nv':
1.3500000000000001, 'fix_trv': True, 'limit_C': (-10.0, 10.0), 'rv':
0.88380402754305054, 'fix_tNv': True, 'fix_MC': False, 'Nsea': 0.0, 'fix_MS':
True, 'limit_trv': (0.0, 8.0), 'limit_rv': (0.0, 8.0), 'alS':
1.1299999999999999, 'fix_bS': True, 'alpS': 0.14999999999999999, 'C':
1.7215515470676213, 'tNv': 0.0, 'limit_Mv': (0.40000000000000002, 1.5),
'limit_tMv': (0.40000000000000002, 2.0), 'limit_bS': (0.40000000000000002,
    5.0), 'fix_alpv': True, 'tbv': 3.0, 'fix_alv': True, 'bv':
0.40000080960497564, 'Mv': 1.4999998470926843, 'fix_rS': True, 'fix_Nv': True,
'fix_tbv': True, 'alpv': 0.84999999999999998, 'MC': 1.9999991700781385,
'fix_Nsea': True, 'fix_alpS': True, 'fix_rv': False, 'limit_bv':
(0.40000000000000002, 5.0), 'limit_MC': (0.40000000000000002, 2.0), 'fix_alS':
True, 'MS': 0.70699999999999996, 'bS': 2.0, 'fix_C': False, 'fix_tMv': True,
'trv': 6.0}
# cutq2 = 2.0 !!!  P(143.10, 155) = 0.7441
ALTGLOcut = {
'fix_THIG': False, 'KAPS': 0.0, 'fix_ALPG': True, 'fix_DELM2S': True,
'fix_SECG': False, 'KAPG': 0.0, 'fix_SECS': False, 'fix_DELM2G': True, 'PS':
2.0, 'SECG': -2.230076249645792, 'fix_AL0S': True, 'fix_THIS': False,
'fix_KAPG': True, 'SECS': -0.41810500770532288, 'fix_M02G': True, 'NG': 0.5,
'fix_ALPS': True, 'M02S': 0.48812880836448502, 'PG': 2.0, 'ALPG':
0.14999999999999999, 'fix_M02S': False, 'NS': 0.15203911208796006, 'fix_AL0G':
True, 'fix_PG': True, 'AL0S': 1.1575060246398083, 'limit_M02S':
(0.29999999999999999, 1.5), 'fix_PS': True, 'SKEWG': 0.0, 'AL0G':
1.2473167010704711, 'SKEWS': 0.0, 'ALPS': 0.14999999999999999, 'M02G':
0.69999999999999996, 'fix_NG': True, 'THIS': 0.080530830450221372, 'fix_SKEWS':
True, 'DELM2S': 0.0, 'fix_SKEWG': True, 'fix_NS': True, 'limit_M02G':
(0.29999999999999999, 1.5), 'fix_KAPS': True, 'THIG': 0.75029383674047267,
'DELM2G': 0.0,
'tMv': 0.80000000000000004, 'rS': 1.0, 'limit_tbv':
(0.40000000000000002, 5.0), 'fix_Mv': False, 'alv':
0.42999999999999999, 'fix_bv': False, 'Nv': 1.3500000000000001,
'fix_trv': True, 'limit_C': (-10.0, 10.0), 'rv': 0.95475017635119874,
'fix_tNv': True, 'fix_MC': False, 'Nsea': 0.0, 'fix_MS': True,
'limit_trv': (0.0, 8.0), 'limit_rv': (0.0, 8.0), 'alS':
1.1299999999999999, 'fix_bS': True, 'alpS': 0.14999999999999999, 'C':
1.5142035872884954, 'tNv': 0.0, 'limit_Mv': (0.40000000000000002, 1.5),
'limit_tMv': (0.40000000000000002, 2.0), 'limit_bS':
(0.40000000000000002, 5.0), 'fix_alpv': True, 'tbv': 3.0, 'fix_alv':
True, 'bv': 0.40000000046746681, 'Mv': 1.0800853752231196, 'fix_rS':
True, 'fix_Nv': True, 'fix_tbv': True, 'alpv': 0.84999999999999998,
'MC': 1.9999999688352545, 'fix_Nsea': True, 'fix_alpS': True, 'fix_rv':
False, 'limit_bv': (0.40000000000000002, 5.0), 'limit_MC':
(0.40000000000000002, 2.0), 'fix_alS': True, 'MS': 0.70699999999999996,
'bS': 2.0, 'fix_C': False, 'fix_tMv': True, 'trv': 6.0}


# Fit to DVCSpoints+data[48]+ALTGLO2points   P(197.43, 178) = 0.1516
# ALTGLO2points = data[5] + data[25] + data[32][18:] + HAD17[::2] + HA17[::2]
# HAD17 = utils.select(data[33], criteria=['Q2 == 1.5', 't == -0.17'])
# HA17 = utils.select(data[34], criteria=['t == -0.17'])
# HA28 = utils.select(data[34], criteria=['t == -0.28'])
#  ('M02S','SECS','SECG', 'THIS', 'THIG', 'rv', 'bv', 'C', 'MC', 'trv', 'tbv')
ALTGLO2 = { 
'fix_THIG': False, 'KAPS': 0.0, 'fix_ALPG': True, 'fix_DELM2S': True,
        'fix_SECG': False, 'KAPG': 0.0, 'fix_SECS': False, 'fix_DELM2G': True,
        'PS': 2.0, 'SECG': -2.6551353141118836, 'fix_AL0S': True, 'fix_THIS':
        False, 'fix_KAPG': True, 'SECS': 0.086814446170543652, 'fix_M02G':
        True, 'NG': 0.5, 'fix_ALPS': True, 'M02S': 0.50518802930909179, 'PG':
        2.0, 'ALPG': 0.14999999999999999, 'fix_M02S': False, 'NS':
        0.15203911208796006, 'fix_AL0G': True, 'fix_PG': True, 'AL0S':
        1.1575060246398083, 'limit_M02S': (0.29999999999999999, 1.5), 'fix_PS':
        True, 'SKEWG': 0.0, 'AL0G': 1.2473167010704711, 'SKEWS': 0.0, 'ALPS':
        0.14999999999999999, 'M02G': 0.69999999999999996, 'fix_NG': True,
        'THIS': -0.072323341948540149, 'fix_SKEWS': True, 'DELM2S': 0.0,
        'fix_SKEWG': True, 'fix_NS': True, 'limit_M02G': (0.29999999999999999,
            1.5), 'fix_KAPS': True, 'THIG': 0.87056201398571653, 'DELM2G': 0.0,
'tMv': 0.80000000000000004, 'rS': 1.0, 'limit_tbv':
        (0.40000000000000002, 5.0), 'fix_Mv': True, 'alv': 0.42999999999999999,
        'fix_bv': False, 'Nv': 1.3500000000000001, 'fix_trv': False, 'limit_C':
        (-10.0, 10.0), 'rv': 0.57894112816586984, 'fix_tNv': True, 'fix_MC':
        False, 'Nsea': 0.0, 'fix_MS': True, 'limit_trv': (0.0, 8.0),
        'limit_rv': (0.0, 8.0), 'alS': 1.1299999999999999, 'fix_bS': True,
        'alpS': 0.14999999999999999, 'C': 9.999987252629829, 'tNv':
        0.59999999999999998, 'limit_Mv': (0.40000000000000002, 1.5),
        'limit_tMv': (0.40000000000000002, 2.0), 'limit_bS':
        (0.40000000000000002, 5.0), 'fix_alpv': True, 'tbv':
        2.7484545877275171, 'fix_alv': True, 'bv': 0.40000055736013956, 'Mv':
        0.80000000000000004, 'fix_rS': True, 'fix_Nv': True, 'fix_tbv': False,
        'alpv': 0.84999999999999998, 'MC': 0.70580723238549492, 'fix_Nsea':
        True, 'fix_alpS': True, 'fix_rv': False, 'limit_bv':
        (0.40000000000000002, 5.0), 'limit_MC': (0.40000000000000002, 2.0),
        'fix_alS': True, 'MS': 0.70699999999999996, 'bS': 2.0, 'fix_C': False,
        'fix_tMv': True, 'trv': 7.871587331745177}


# as ALTGLO2 but, starting from hy1THI  chisq/dof = 195/176
ALTGLO2B = {
'fix_THIG': False, 'KAPS': 0.0, 'fix_ALPG': True, 'fix_DELM2S': True,
  'fix_SECG': False, 'KAPG': 0.0, 'fix_SECS': False, 'fix_DELM2G': True,
  'PS': 2.0, 'SECG': -2.6338224450174748, 'fix_AL0S': True, 'fix_THIS':
  False, 'fix_KAPG': True, 'SECS': -0.056898607234204653, 'fix_M02G':
  True, 'NG': 0.5, 'fix_ALPS': True, 'M02S': 0.50156492842554723, 'PG':
  2.0, 'ALPG': 0.14999999999999999, 'fix_M02S': False, 'NS':
  0.15203911208796006, 'fix_AL0G': True, 'fix_PG': True, 'AL0S':
  1.1575060246398083, 'limit_M02S': (0.29999999999999999, 1.5), 'fix_PS':
  True, 'SKEWG': 0.0, 'AL0G': 1.2473167010704711, 'SKEWS': 0.0, 'ALPS':
  0.14999999999999999, 'M02G': 0.69999999999999996, 'fix_NG': True,
  'THIS': -0.028774918377538185, 'fix_SKEWS': True, 'DELM2S': 0.0,
  'fix_SKEWG': True, 'fix_NS': True, 'limit_M02G': (0.29999999999999999,
      1.5), 'fix_KAPS': True, 'THIG': 0.88660385756192783, 'DELM2G': 0.0,
'tMv': 0.72300855154777954, 'rS': 1.0, 'limit_tbv':
  (0.40000000000000002, 5.0), 'fix_Mv': False, 'alv':
  0.42999999999999999, 'fix_bv': False, 'Nv': 1.3500000000000001,
  'fix_trv': False, 'limit_C': (-10.0, 10.0), 'rv': 0.48857257980811708,
  'fix_tNv': True, 'fix_MC': False, 'Nsea': 0.0, 'fix_MS': True,
  'limit_trv': (0.0, 8.0), 'limit_rv': (0.0, 8.0), 'alS':
  1.1299999999999999, 'fix_bS': True, 'alpS': 0.14999999999999999, 'C':
  9.9999973166591687, 'tNv': 0.59999999999999998, 'limit_Mv':
  (0.40000000000000002, 1.5), 'limit_tMv': (0.40000000000000002, 2.0),
  'limit_bS': (0.40000000000000002, 5.0), 'fix_alpv': True, 'tbv':
  2.5968261726412112, 'fix_alv': True, 'bv': 0.43965295355048556, 'Mv':
  1.4999996661106532, 'fix_rS': True, 'fix_Nv': True, 'fix_tbv': False,
  'alpv': 0.84999999999999998, 'MC': 0.73870986675252859, 'fix_Nsea':
  True, 'fix_alpS': True, 'fix_rv': False, 'limit_bv':
  (0.40000000000000002, 5.0), 'limit_MC': (0.40000000000000002, 2.0),
  'fix_alS': True, 'MS': 0.70699999999999996, 'bS': 2.0, 'fix_C': False,
  'fix_tMv': False, 'trv': 7.9999982364753528
}

# as ALTGLO2 but starting from hy1THI AND keeping SECS fixed  chisq/dof = 196/177
ALTGLO2C = {
'fix_THIG': False, 'KAPS': 0.0, 'fix_ALPG': True, 'fix_DELM2S': True,
'fix_SECG': False, 'KAPG': 0.0, 'fix_SECS': True, 'fix_DELM2G': True, 'PS':
2.0, 'SECG': -2.5440702464099867, 'fix_AL0S': True, 'fix_THIS': False,
'fix_KAPG': True, 'SECS': -0.44159594968364457, 'fix_M02G': True, 'NG': 0.5,
'fix_ALPS': True, 'M02S': 0.49720249608913697, 'PG': 2.0, 'ALPG':
0.14999999999999999, 'fix_M02S': False, 'NS': 0.15203911208796006, 'fix_AL0G':
True, 'fix_PG': True, 'AL0S': 1.1575060246398083, 'limit_M02S':
(0.29999999999999999, 1.5), 'fix_PS': True, 'SKEWG': 0.0, 'AL0G':
1.2473167010704711, 'SKEWS': 0.0, 'ALPS': 0.14999999999999999, 'M02G':
0.69999999999999996, 'fix_NG': True, 'THIS': 0.087583283215139196, 'fix_SKEWS':
True, 'DELM2S': 0.0, 'fix_SKEWG': True, 'fix_NS': True, 'limit_M02G':
(0.29999999999999999, 1.5), 'fix_KAPS': True, 'THIG': 0.90870318050138843,
'DELM2G': 0.0,
'tMv': 0.715584941355635, 'rS': 1.0, 'limit_tbv': (0.40000000000000002, 5.0),
'fix_Mv': False, 'alv': 0.42999999999999999, 'fix_bv': False, 'Nv':
1.3500000000000001, 'fix_trv': False, 'limit_C': (-10.0, 10.0), 'rv':
0.45374522143316343, 'fix_tNv': True, 'fix_MC': False, 'Nsea': 0.0, 'fix_MS':
True, 'limit_trv': (0.0, 8.0), 'limit_rv': (0.0, 8.0), 'alS':
1.1299999999999999, 'fix_bS': True, 'alpS': 0.14999999999999999, 'C':
9.9999999527870074, 'tNv': 0.59999999999999998, 'limit_Mv':
(0.40000000000000002, 1.5), 'limit_tMv': (0.40000000000000002, 2.0),
'limit_bS': (0.40000000000000002, 5.0), 'fix_alpv': True, 'tbv':
2.5636206794234195, 'fix_alv': True, 'bv': 0.40002654542014571, 'Mv':
1.499999483194582, 'fix_rS': True, 'fix_Nv': True, 'fix_tbv': False, 'alpv':
0.84999999999999998, 'MC': 0.73974790427385917, 'fix_Nsea': True, 'fix_alpS':
True, 'fix_rv': False, 'limit_bv': (0.40000000000000002, 5.0), 'limit_MC':
(0.40000000000000002, 2.0), 'fix_alS': True, 'MS': 0.70699999999999996, 'bS':
2.0, 'fix_C': False, 'fix_tMv': False, 'trv': 7.9999998367016341}
  
# chisq/dof = 1.36
#{"rv" -> 1.24, "Mv" -> 0.863, "bv" -> 0.531, "bs" -> 5.77, "d" -> -2.,
#"Mc" -> 0.28, "rtv" -> 0.559, "Mtv" -> 0.42, "btv" -> 1.17, "rpi" -> 4.75}
DMPP = {
  'C': 1.7794,
 'MC': 2.1333,
 'MS': 0.7071067811865476,        
 'Mv': 0.8,
 'NS': 1.5,                           
 'Nv': 1.3500000000000001,            
 'alS': 1.1299999999999999,           
 'alpS': 0.14999999999999999,         
 'alpv': 0.84999999999999998,         
 'alv': 0.42999999999999999,          
 'bS': 3.6652,
 'bv': 0.7411,
 'rS': 1.0,                           
 'rv': 1.0084,
 'tMv': 1.8928,
 'tNv': 0.6,                          
 'tbv': 1.0168,
 'trv': 0.8715,
 'rpi': 4.1352,
 'Mpi': 1.4634}

