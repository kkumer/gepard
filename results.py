
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
# Or is it 108/141?
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
KM10a = ALTGLO
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

# chisq/dof = 200/172 = 0.07
# total uncut unp (my CLAS)
# 15 params DRPP + gepard sea
KKunpfull = {'tMv': 3.9999968268035802, 'rS': 1.0, 'limit_tbv':
    (0.40000000000000002, 5.0), 'rpi': 4.1889940690558589, 'fix_Mv': True,
    'fix_Mpi': False, 'alv': 0.42999999999999999, 'fix_bv': False, 'Mpi':
    0.40000343322644777, 'Nv': 1.3500000000000001, 'fix_trv': False, 'limit_C':
    (-10.0, 10.0), 'rv': 0.54769484516465106, 'fix_tNv': True, 'fix_MC': False,
    'Nsea': 0.0, 'fix_MS': True, 'limit_trv': (0.0, 8.0), 'limit_rv': (0.0,
        8.0), 'alS': 1.1299999999999999, 'fix_bS': True, 'alpS':
    0.14999999999999999, 'limit_Mpi': (0.40000000000000002, 4.0), 'C':
    7.907821229500307, 'tNv': 0.59999999999999998, 'limit_Mv':
    (0.40000000000000002, 4.0), 'limit_tMv': (0.40000000000000002, 4.0),
    'limit_bS': (0.40000000000000002, 5.0), 'fix_alpv': True, 'tbv':
    0.75365091725160882, 'fix_alv': True, 'limit_rpi': (-8, 8.0), 'bv':
    0.77540603437335276, 'Mv': 3.9999999240306829, 'fix_rS': True, 'fix_Nv':
    True, 'fix_tbv': False, 'fix_rpi': False, 'alpv': 0.84999999999999998,
    'MC': 0.98980908992269756, 'fix_Nsea': True, 'fix_alpS': True, 'fix_rv':
    False, 'limit_bv': (0.40000000000000002, 5.0), 'limit_MC':
    (0.40000000000000002, 4.0), 'fix_alS': True, 'MS': 0.70699999999999996,
    'bS': 2.0, 'fix_C': False, 'fix_tMv': True, 'trv': 4.1649860987845884,
'fix_THIG': False, 'KAPS': 0.0, 'fix_ALPG': True,
'fix_DELM2S': True, 'fix_SECG': False, 'KAPG': 0.0, 'fix_SECS':
False, 'fix_DELM2G': True, 'PS': 2.0, 'SECG': -2.7222253774417884,
'fix_AL0S': True, 'fix_THIS': False, 'fix_KAPG': True, 'SECS':
0.28377840745999455, 'fix_M02G': True, 'NG': 0.5, 'fix_ALPS': True,
'M02S': 0.51824628465962441, 'PG': 2.0, 'ALPG':
0.14999999999999999, 'fix_M02S': False, 'NS': 0.15203911208796006,
'fix_AL0G': True, 'fix_PG': True, 'AL0S': 1.1575060246398083,
'limit_M02S': (0.29999999999999999, 1.5), 'fix_PS': True, 'SKEWG':
0.0, 'AL0G': 1.2473167010704711, 'SKEWS': 0.0, 'ALPS':
0.14999999999999999, 'M02G': 0.69999999999999996, 'fix_NG': True,
'THIS': -0.13244402186832566, 'fix_SKEWS': True, 'DELM2S': 0.0,
'fix_SKEWG': True, 'fix_NS': True, 'limit_M02G':
(0.29999999999999999, 1.5), 'fix_KAPS': True, 'THIG':
0.86655265321494723, 'DELM2G': 0.0}


# chisq/dof = 194/163 = 0.05
# unp + TSA1 cut Q2min=1.6 GeV^2
# 15 params DRPP + gepard sea
KKunpTSAcut = {
'tMv': 1.6005309774046621, 'rS': 1.0, 'limit_tbv':
(0.40000000000000002, 5.0), 'rpi': 3.4505551335493436, 'fix_Mv': False,
'fix_Mpi': False, 'alv': 0.42999999999999999, 'fix_bv': False, 'Mpi':
3.9999999903396257, 'Nv': 1.3500000000000001, 'fix_trv': False,
'limit_C': (-10.0, 10.0), 'rv': 0.96621876134621099, 'fix_tNv': True,
'fix_MC': False, 'Nsea': 0.0, 'fix_MS': True, 'limit_trv': (0.0, 8.0),
'limit_rv': (0.0, 8.0), 'alS': 1.1299999999999999, 'fix_bS': True,
'alpS': 0.14999999999999999, 'limit_Mpi': (0.40000000000000002, 4.0),
'C': 1.4557762594377621, 'tNv': 0.59999999999999998, 'limit_Mv':
(0.40000000000000002, 4.0), 'limit_tMv': (0.40000000000000002, 4.0),
'limit_bS': (0.40000000000000002, 5.0), 'fix_alpv': True, 'tbv':
0.4000000000269115, 'fix_alv': True, 'limit_rpi': (-8, 8.0), 'bv':
0.40000000000095121, 'Mv': 3.9999999723207176, 'fix_rS': True,
'fix_Nv': True, 'fix_tbv': False, 'fix_rpi': False, 'alpv':
0.84999999999999998, 'MC': 3.9999999991266391, 'fix_Nsea': True,
'fix_alpS': True, 'fix_rv': False, 'limit_bv': (0.40000000000000002,
    5.0), 'limit_MC': (0.40000000000000002, 4.0), 'fix_alS': True,
'MS': 0.70699999999999996, 'bS': 2.0, 'fix_C': False, 'fix_tMv': False,
'trv': 0.89909150018653072,
'fix_THIG': False, 'KAPS': 0.0, 'fix_ALPG': True, 'fix_DELM2S': True,
'fix_SECG': False, 'KAPG': 0.0, 'fix_SECS': False, 'fix_DELM2G': True,
'PS': 2.0, 'SECG': -3.7004671516071683, 'fix_AL0S': True, 'fix_THIS':
False, 'fix_KAPG': True, 'SECS': 1.2202011017848651, 'fix_M02G': True,
'NG': 0.5, 'fix_ALPS': True, 'M02S': 0.5325145484619549, 'PG': 2.0,
'ALPG': 0.14999999999999999, 'fix_M02S': False, 'NS':
0.15203911208796006, 'fix_AL0G': True, 'fix_PG': True, 'AL0S':
1.1575060246398083, 'limit_M02S': (0.29999999999999999, 1.5), 'fix_PS':
True, 'SKEWG': 0.0, 'AL0G': 1.2473167010704711, 'SKEWS': 0.0, 'ALPS':
0.14999999999999999, 'M02G': 0.69999999999999996, 'fix_NG': True,
'THIS': -0.41512377901430425, 'fix_SKEWS': True, 'DELM2S': 0.0,
'fix_SKEWG': True, 'fix_NS': True, 'limit_M02G': (0.29999999999999999,
    1.5), 'fix_KAPS': True, 'THIG': 1.1755426941813287, 'DELM2G': 0.0}

# P(131.94, 160) = 0.9488
# H1ZEUSpoints+UNP5points
# 15 params DRPP + gepard sea
#
#   Mv = 4.00 +- 3.33 (edge)
#   rv = 0.62 +- 0.06
#   bv = 0.40 +- 0.67
#    C = 8.78 +- 0.98
#   MC = 0.97 +- 0.11
#  tMv = 0.88 +- 0.24
#  trv = 7.76 +- 1.39
#  tbv = 2.05 +- 0.40
#  rpi = 3.54 +- 1.77
#  Mpi = 0.73 +- 0.37
# M02S = 0.51 +- 0.02
# SECS = 0.28 +- 0.02
# THIS = -0.13 +- 0.01
# SECG = -2.79 +- 0.12
# THIG = 0.90 +- 0.05
KKunp5 = {
'tMv': 0.88386035557556641, 'rS': 1.0, 'limit_tbv': (0.40000000000000002, 5.0),
'rpi': 3.5355742996824659, 'fix_Mv': False, 'fix_Mpi': False, 'alv':
0.42999999999999999, 'fix_bv': False, 'Mpi': 0.72673907227274226, 'Nv':
1.3500000000000001, 'fix_trv': False, 'limit_C': (-10.0, 10.0), 'rv':
0.61981779615528687, 'fix_tNv': True, 'fix_MC': False, 'Nsea': 0.0, 'fix_MS':
True, 'limit_trv': (0.0, 8.0), 'limit_rv': (0.0, 8.0), 'alS':
1.1299999999999999, 'fix_bS': True, 'alpS': 0.14999999999999999, 'limit_Mpi':
(0.40000000000000002, 4.0), 'C': 8.7771705279055432, 'tNv':
0.59999999999999998, 'limit_Mv': (0.40000000000000002, 4.0), 'limit_tMv':
(0.40000000000000002, 4.0), 'limit_bS': (0.40000000000000002, 5.0), 'fix_alpv':
True, 'tbv': 2.0496515976221952, 'fix_alv': True, 'limit_rpi': (-8, 8.0), 'bv':
0.40375381018570683, 'Mv': 3.999996566773552, 'fix_rS': True, 'fix_Nv': True,
'fix_tbv': False, 'fix_rpi': False, 'alpv': 0.84999999999999998, 'MC':
0.9746453693853796, 'fix_Nsea': True, 'fix_alpS': True, 'fix_rv': False,
'limit_bv': (0.40000000000000002, 5.0), 'limit_MC': (0.40000000000000002, 4.0),
'fix_alS': True, 'MS': 0.70699999999999996, 'bS': 2.0, 'fix_C': False,
'fix_tMv': False, 'trv': 7.7592155354071064,
'fix_THIG': False, 'KAPS': 0.0, 'fix_ALPG': True, 'fix_DELM2S': True,
'fix_SECG': False, 'KAPG': 0.0, 'fix_SECS': False, 'fix_DELM2G': True,
'PS': 2.0, 'SECG': -2.787206139384161, 'fix_AL0S': True, 'fix_THIS':
False, 'fix_KAPG': True, 'SECS': 0.27800486286895021, 'fix_M02G': True,
'NG': 0.5, 'fix_ALPS': True, 'M02S': 0.51318802297356658, 'PG': 2.0,
'ALPG': 0.14999999999999999, 'fix_M02S': False, 'NS':
0.15203911208796006, 'fix_AL0G': True, 'fix_PG': True, 'AL0S':
1.1575060246398083, 'limit_M02S': (0.29999999999999999, 1.5), 'fix_PS':
True, 'SKEWG': 0.0, 'AL0G': 1.2473167010704711, 'SKEWS': 0.0, 'ALPS':
0.14999999999999999, 'M02G': 0.69999999999999996, 'fix_NG': True,
'THIS': -0.12999801345190121, 'fix_SKEWS': True, 'DELM2S': 0.0,
'fix_SKEWG': True, 'fix_NS': True, 'limit_M02G': (0.29999999999999999,
1.5), 'fix_KAPS': True, 'THIG': 0.89633675296304127, 'DELM2G': 0.0}
KM10 = KKunp5

# P(196.41, 172) = 0.0978
# H1ZEUSpoints+UNP5points+TSA1points
# 
#   Mv = 4.00 +- 3.54 (edge)
#   rv = 1.07 +- 0.04
#   bv = 0.40 +- 0.02 (edge)
#    C = 1.05 +- 0.30
#   MC = 4.00 +- 3.38 (edge)
#  tMv = 1.32 +- 2.26
#  trv = 0.82 +- 0.19
#  tbv = 0.40 +- 0.16 (edge)
#  rpi = 3.38 +- 0.16
#  Mpi = 4.00 +- 2.33 (edge)
# M02S = 0.54 +- 0.02
# SECS = 1.49 +- 0.02
# THIS = -0.50 +- 0.01
# SECG = -3.34 +- 0.12
# THIG = 0.94 +- 0.05
KKunpTSA1 = {
'tMv': 1.3156138534382904, 'rS': 1.0, 'limit_tbv':
    (0.40000000000000002, 5.0), 'rpi': 3.3825826231632483, 'fix_Mv': False,
    'fix_Mpi': False, 'alv': 0.42999999999999999, 'fix_bv': False, 'Mpi':
    3.999996566773552, 'Nv': 1.3500000000000001, 'fix_trv': False, 'limit_C':
    (-10.0, 10.0), 'rv': 1.073018581431247, 'fix_tNv': True, 'fix_MC': False,
    'Nsea': 0.0, 'fix_MS': True, 'limit_trv': (0.0, 8.0), 'limit_rv': (0.0,
        8.0), 'alS': 1.1299999999999999, 'fix_bS': True, 'alpS':
    0.14999999999999999, 'limit_Mpi': (0.40000000000000002, 4.0), 'C':
    1.0467670088984438, 'tNv': 0.59999999999999998, 'limit_Mv':
    (0.40000000000000002, 4.0), 'limit_tMv': (0.40000000000000002, 4.0),
    'limit_bS': (0.40000000000000002, 5.0), 'fix_alpv': True, 'tbv':
    0.40000438690046103, 'fix_alv': True, 'limit_rpi': (-8, 8.0), 'bv':
    0.40000438690046103, 'Mv': 3.9999999070113295, 'fix_rS': True, 'fix_Nv':
    True, 'fix_tbv': False, 'fix_rpi': False, 'alpv': 0.84999999999999998,
    'MC': 3.9999983844050448, 'fix_Nsea': True, 'fix_alpS': True, 'fix_rv':
    False, 'limit_bv': (0.40000000000000002, 5.0), 'limit_MC':
    (0.40000000000000002, 4.0), 'fix_alS': True, 'MS': 0.70699999999999996,
    'bS': 2.0, 'fix_C': False, 'fix_tMv': False, 'trv': 0.82031179907880603,
'fix_THIG': False, 'KAPS': 0.0, 'fix_ALPG': True, 'fix_DELM2S': True,
    'fix_SECG': False, 'KAPG': 0.0, 'fix_SECS': False, 'fix_DELM2G': True,
    'PS': 2.0, 'SECG': -3.3409075255124683, 'fix_AL0S': True, 'fix_THIS':
    False, 'fix_KAPG': True, 'SECS': 1.4901840478689927, 'fix_M02G': True,
    'NG': 0.5, 'fix_ALPS': True, 'M02S': 0.53728944647200971, 'PG': 2.0,
    'ALPG': 0.14999999999999999, 'fix_M02S': False, 'NS': 0.15203911208796006,
    'fix_AL0G': True, 'fix_PG': True, 'AL0S': 1.1575060246398083, 'limit_M02S':
    (0.29999999999999999, 1.5), 'fix_PS': True, 'SKEWG': 0.0, 'AL0G':
    1.2473167010704711, 'SKEWS': 0.0, 'ALPS': 0.14999999999999999, 'M02G':
    0.69999999999999996, 'fix_NG': True, 'THIS': -0.49729314108803641,
    'fix_SKEWS': True, 'DELM2S': 0.0, 'fix_SKEWG': True, 'fix_NS': True,
    'limit_M02G': (0.29999999999999999, 1.5), 'fix_KAPS': True, 'THIG':
    0.94501600044624157, 'DELM2G': 0.0}


# P(168.82, 156) = 0.2283
# utils.select(H1ZEUSpoints+UNP5points+TSA1points, criteria=['Q2>=1.6'])
#
#   Mv = 4.00 +- 3.58 (edge)
#   rv = 1.03 +- 0.04
#   bv = 0.40 +- 0.03 (edge)
#    C = 1.23 +- 0.33
#   MC = 4.00 +- 3.36 (edge)
#  tMv = 1.03 +- 0.82
#  trv = 0.92 +- 0.23
#  tbv = 0.40 +- 0.33 (edge)
#  rpi = 3.38 +- 0.17
#  Mpi = 4.00 +- 2.35 (edge)
# M02S = 0.52 +- 0.02
# SECS = 0.57 +- 0.03
# THIS = -0.22 +- 0.01
# SECG = -3.30 +- 0.18
# THIG = 1.09 +- 0.09
KKunpTSA1cut16 = {
'tMv': 1.0256490851163997, 'rS': 1.0, 'limit_tbv':
        (0.40000000000000002, 5.0), 'rpi': 3.3811271921721318, 'fix_Mv': False,
        'fix_Mpi': False, 'alv': 0.42999999999999999, 'fix_bv': False, 'Mpi':
        3.999996566773552, 'Nv': 1.3500000000000001, 'fix_trv': False,
        'limit_C': (-10.0, 10.0), 'rv': 1.0274081985087875, 'fix_tNv': True,
        'fix_MC': False, 'Nsea': 0.0, 'fix_MS': True, 'limit_trv': (0.0, 8.0),
        'limit_rv': (0.0, 8.0), 'alS': 1.1299999999999999, 'fix_bS': True,
        'alpS': 0.14999999999999999, 'limit_Mpi': (0.40000000000000002, 4.0),
        'C': 1.2299651857374361, 'tNv': 0.59999999999999998, 'limit_Mv':
        (0.40000000000000002, 4.0), 'limit_tMv': (0.40000000000000002, 4.0),
        'limit_bS': (0.40000000000000002, 5.0), 'fix_alpv': True, 'tbv':
        0.40000438690046103, 'fix_alv': True, 'limit_rpi': (-8, 8.0), 'bv':
        0.40000438690046103, 'Mv': 3.9999941114628839, 'fix_rS': True,
        'fix_Nv': True, 'fix_tbv': False, 'fix_rpi': False, 'alpv':
        0.84999999999999998, 'MC': 3.9999999029164206, 'fix_Nsea': True,
        'fix_alpS': True, 'fix_rv': False, 'limit_bv': (0.40000000000000002,
            5.0), 'limit_MC': (0.40000000000000002, 4.0), 'fix_alS': True,
        'MS': 0.70699999999999996, 'bS': 2.0, 'fix_C': False, 'fix_tMv': False,
        'trv': 0.92231447360569696,
'fix_THIG': False, 'KAPS': 0.0,
        'fix_ALPG': True, 'fix_DELM2S': True, 'fix_SECG': False,
        'KAPG': 0.0, 'fix_SECS': False, 'fix_DELM2G': True, 'PS': 2.0,
        'SECG': -3.2970087421085648, 'fix_AL0S': True, 'fix_THIS':
        False, 'fix_KAPG': True, 'SECS': 0.56506108003883404,
        'fix_M02G': True, 'NG': 0.5, 'fix_ALPS': True, 'M02S':
        0.52222295960573129, 'PG': 2.0, 'ALPG': 0.14999999999999999,
        'fix_M02S': False, 'NS': 0.15203911208796006, 'fix_AL0G': True,
        'fix_PG': True, 'AL0S': 1.1575060246398083, 'limit_M02S':
        (0.29999999999999999, 1.5), 'fix_PS': True, 'SKEWG': 0.0,
        'AL0G': 1.2473167010704711, 'SKEWS': 0.0, 'ALPS':
        0.14999999999999999, 'M02G': 0.69999999999999996, 'fix_NG':
        True, 'THIS': -0.21683562189055644, 'fix_SKEWS': True,
        'DELM2S': 0.0, 'fix_SKEWG': True, 'fix_NS': True, 'limit_M02G':
        (0.29999999999999999, 1.5), 'fix_KAPS': True, 'THIG':
        1.0922767800778619, 'DELM2G': 0.0}

# 86/127
# DVCSpoints+GLOpoints
# t.m.release_parameters('M02S', 'SECS', 'SECG', 'THIS', 'THIG', 'rv', 'bv', 'Mv', 'C', 'MC', 'trv', 'tbv', 'tMv')
# mDRsea = Model.ComptonModelDRsea()
#m = Model.HybridDipole(mGepard, mDRsea)
#t = Approach.hotfixedBMK(m)
KM10aCand1 = {
'tMv': 0.72866940517442902, 'rS': 1.0, 'limit_tbv': (0.40000000000000002, 5.0),
'fix_Mv': False, 'alv': 0.42999999999999999, 'fix_bv': False, 'Nv':
1.3500000000000001, 'fix_trv': True, 'limit_C': (-10.0, 10.0), 'rv':
0.49042198904829348, 'fix_tNv': True, 'fix_MC': False, 'Nsea': 0.0, 'fix_MS':
True, 'limit_trv': (0.0, 8.0), 'limit_rv': (0.0, 8.0), 'alS':
1.1299999999999999, 'fix_bS': True, 'alpS': 0.14999999999999999, 'C':
9.1754359638461835, 'tNv': 0.59999999999999998, 'limit_Mv':
(0.40000000000000002, 1.5), 'limit_tMv': (0.40000000000000002, 2.0),
'limit_bS': (0.40000000000000002, 5.0), 'fix_alpv': True, 'tbv':
3.1820059611233771, 'fix_alv': True, 'bv': 0.40000252441118156, 'Mv':
1.4999996092274657, 'fix_rS': True, 'fix_Nv': True, 'fix_tbv': True, 'alpv':
0.84999999999999998, 'MC': 0.74710776531169398, 'fix_Nsea': True, 'fix_alpS':
True, 'fix_rv': False, 'limit_bv': (0.40000000000000002, 5.0), 'limit_MC':
(0.40000000000000002, 2.0), 'fix_alS': True, 'MS': 0.70699999999999996, 'bS':
2.0, 'fix_C': False, 'fix_tMv': True, 'trv': 7.9999646375514626,
'fix_THIG': False, 'KAPS': 0.0, 'fix_ALPG': True, 'fix_DELM2S': True,
'fix_SECG': False, 'KAPG': 0.0, 'fix_SECS': False, 'fix_DELM2G': True, 'PS':
2.0, 'SECG': -2.0577533684495162, 'fix_AL0S': True, 'fix_THIS': False,
'fix_KAPG': True, 'SECS': -1.0722853854418892, 'fix_M02G': True, 'NG': 0.5,
'fix_ALPS': True, 'M02S': 0.48024712029919264, 'PG': 2.0, 'ALPG':
0.14999999999999999, 'fix_M02S': False, 'NS': 0.15203911208796006, 'fix_AL0G':
True, 'fix_PG': True, 'AL0S': 1.1575060246398083, 'limit_M02S':
(0.29999999999999999, 1.5), 'fix_PS': True, 'SKEWG': 0.0, 'AL0G':
1.2473167010704711, 'SKEWS': 0.0, 'ALPS': 0.14999999999999999, 'M02G':
0.69999999999999996, 'fix_NG': True, 'THIS': 0.27563608304019077, 'fix_SKEWS':
True, 'DELM2S': 0.0, 'fix_SKEWG': True, 'fix_NS': True, 'limit_M02G':
(0.29999999999999999, 1.5), 'fix_KAPS': True, 'THIG': 0.8023193873056319,
'DELM2G': 0.0}

# Hybrid fit + 3rd PW 
KM10b = {
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
 'C': 5.4259,
 'MC': 1.3305,
 'MS': 0.70699999999999996,
 'Mv': 0.8,
 'Nsea': 0.0,
 'Nv': 1.3500000000000001,
 'alS': 1.1299999999999999,
 'alpS': 0.14999999999999999,
 'alpv': 0.84999999999999998,
 'alv': 0.42999999999999999,
 'bS': 2.0,
 'bv': 0.7706,
 'rS': 1.0,
 'rv': 0.8081,
 'tMv': 0.8,
 'tNv': 0.59999999999999998,
 'tbv': 1.0,
 'trv': 3.2931,
'rpi': 4.0201, 
'Mpi': 1.5369}

KM10cov = {('MC', 'Mv'): 0.00046177257426809835, ('trv', 'SECG'):
-0.00013803569735782616, ('C', 'tMv'): 0.0076173475635255179, ('SECS', 'THIS'):
-9.2666461860067169e-05, ('Mv', 'Mv'): 0.030136677951808178, ('SECG', 'MC'):
-0.00043428438612909588, ('M02S', 'bv'): 2.5533094568802685e-07, ('rpi', 'Mv'):
0.0058237665801370939, ('Mv', 'tbv'): 0.01047298222366057, ('MC', 'rv'):
0.0012919571164652522, ('trv', 'C'): 0.024580425390947118, ('bv', 'M02S'):
2.5533094568802685e-07, ('C', 'Mpi'): 0.033533953127871163, ('rpi', 'THIS'):
-4.4512209185522965e-05, ('THIG', 'tMv'): -0.00036813684294659244, ('tMv',
'Mv'): 0.0007058811798656881, ('SECG', 'rpi'): -0.0059901485556498749, ('SECG',
'bv'): 5.5867312893547216e-07, ('tbv', 'rv'): 0.016138173428212971, ('Mpi',
'rv'): -0.0066859466473994369, ('C', 'M02S'): -0.00016835768404170823, ('MC',
'SECG'): -0.00043428438612909588, ('tbv', 'tMv'): 0.018425480149743462, ('rpi',
'Mpi'): -1.1410704284544626, ('THIS', 'Mv'): -5.4510004097628973e-06, ('bv',
'MC'): 2.7414566836294294e-05, ('bv', 'THIG'): 6.1210208869426482e-07, ('SECS',
'SECS'): 0.00031741183284160766, ('rv', 'Mpi'): -0.0066859466473994369,
('THIS', 'MC'): -3.8633960535166563e-06, ('rpi', 'MC'): -0.084847904114066064,
('trv', 'trv'): 0.14004092999796483, ('THIS', 'M02S'): -1.2750854559960807e-05,
('rv', 'MC'): 0.0012919571164652522, ('trv', 'bv'): -6.45996002760603e-05,
('bv', 'C'): 0.00049573984495388098, ('SECS', 'Mpi'): -4.9995881049174269e-05,
('rpi', 'tMv'): -0.11809057869871689, ('tbv', 'THIG'): -0.00018332566290391429,
('THIG', 'C'): -0.0037465604426767055, ('MC', 'bv'): 2.7414566836294294e-05,
('M02S', 'THIG'): 0.0001240871438872065, ('C', 'rv'): -0.01913536327848691,
('tbv', 'Mpi'): -0.0094233508455311803, ('THIS', 'C'): 7.3981105852644218e-05,
('SECG', 'THIS'): 1.8170497729357869e-05, ('MC', 'tMv'):
0.00089859952414649837, ('C', 'THIS'): 7.3981105852644218e-05, ('SECG', 'trv'):
-0.00013803569735782616, ('trv', 'tbv'): 0.040071096500902022, ('THIG',
'SECG'): -0.0056477913742646072, ('Mpi', 'THIG'): -0.0012774339177200611,
('Mpi', 'MC'): 0.016774499369639435, ('M02S', 'SECG'): -0.00065274600517088483,
('Mpi', 'THIS'): 1.9452061103007188e-05, ('bv', 'bv'): 0.00032781113091329291,
('bv', 'THIS'): 3.3483224714976356e-08, ('M02S', 'rv'): 1.6966102046237698e-05,
('Mpi', 'Mpi'): 0.25944010783314697, ('SECG', 'C'): 0.0062260289244321773,
('THIG', 'M02S'): 0.0001240871438872065, ('THIG', 'trv'):
1.9259745629808435e-05, ('THIS', 'rpi'): -4.4512209185522965e-05, ('tMv',
'THIS'): 4.578506313793805e-07, ('rv', 'THIS'): -9.9195203113994761e-06,
('tMv', 'trv'): -0.017112636595746616, ('Mpi', 'tMv'): 0.023924151157681827,
('bv', 'tMv'): -4.0455184030838398e-05, ('tMv', 'rpi'): -0.11809057869871689,
('bv', 'SECS'): -3.8554773538420096e-07, ('tbv', 'Mv'): 0.01047298222366057,
('SECG', 'tbv'): -0.00038966065116295923, ('Mv', 'THIS'):
-5.4510004097628973e-06, ('Mv', 'bv'): -1.1999255463837836e-05, ('tMv',
'SECS'): 1.0454724955636238e-05, ('tbv', 'rpi'): -0.077531796761015792, ('MC',
'THIS'): -3.8633960535166563e-06, ('trv', 'SECS'): 3.4042952793082553e-06,
('SECS', 'C'): -1.4277914986737694e-05, ('bv', 'Mpi'): 0.00016645045552772029,
('M02S', 'C'): -0.00016835768404170823, ('THIG', 'THIG'): 0.002596968207828061,
('MC', 'trv'): 0.0010161713431604599, ('Mpi', 'Mv'): -0.0053185636966613676,
('THIG', 'rv'): 0.00010520521114085837, ('trv', 'rpi'): -0.048045475169797416,
('rpi', 'rv'): 0.0031353718480231782, ('SECG', 'rv'): -0.00036243927535111489,
('M02S', 'SECS'): -5.0430389084205586e-05, ('tbv', 'tbv'):
0.099753313287388568, ('THIG', 'Mv'): 0.00049377196297537256, ('SECG', 'Mv'):
-0.00087216225032152686, ('MC', 'SECS'): -1.2520130272335073e-05, ('THIG',
'bv'): 6.1210208869426482e-07, ('M02S', 'Mpi'): 0.00022850003340474631, ('Mv',
'C'): -0.00054377989359903124, ('C', 'MC'): -0.036258588372149528, ('tbv',
'M02S'): 0.00012051361385844407, ('rv', 'rpi'): 0.0031353718480231782, ('MC',
'tbv'): 0.0061058836885243732, ('THIS', 'tbv'): -1.9388200839129332e-05, ('rv',
'tMv'): 0.0018211902078971432, ('Mv', 'THIG'): 0.00049377196297537256, ('Mv',
'MC'): 0.00046177257426809835, ('rv', 'bv'): 2.635554115731397e-07, ('SECS',
'rpi'): 0.00012872560338281992, ('THIG', 'tbv'): -0.00018332566290391429,
('tMv', 'MC'): 0.00089859952414649837, ('tbv', 'trv'): 0.040071096500902022,
('Mv', 'Mpi'): -0.0053185636966613676, ('tMv', 'rv'): 0.0018211902078971432,
('tbv', 'THIS'): -1.9388200839129332e-05, ('trv', 'rv'): 0.004166816278100764,
('C', 'tbv'): -0.047093749916571746, ('THIS', 'THIS'): 2.8871715759127388e-05,
('SECG', 'Mpi'): 0.002219969458552852, ('THIG', 'SECS'):
-5.4630238622194277e-06, ('Mv', 'SECS'): 1.7414693760166068e-05, ('rpi',
'SECS'): 0.00012872560338281992, ('Mv', 'tMv'): 0.0007058811798656881, ('C',
'bv'): 0.00049573984495388098, ('bv', 'tbv'): -0.00033210010411835868, ('SECG',
'THIG'): -0.0056477913742646072, ('SECS', 'bv'): -3.8554773538420096e-07,
('Mpi', 'trv'): 0.015250653525779433, ('rpi', 'tbv'): -0.077531796761015792,
('Mv', 'rpi'): 0.0058237665801370939, ('tbv', 'SECG'): -0.00038966065116295923,
('tMv', 'bv'): -4.0455184030838398e-05, ('bv', 'Mv'): -1.1999255463837836e-05,
('Mpi', 'bv'): 0.00016645045552772029, ('rv', 'C'): -0.01913536327848691,
('Mpi', 'rpi'): -1.1410704284544626, ('rpi', 'M02S'): -0.00083139401293003622,
('tMv', 'M02S'): -1.6647261533993509e-05, ('trv', 'M02S'):
9.9243841643798775e-05, ('C', 'trv'): 0.024580425390947118, ('tMv', 'THIG'):
-0.00036813684294659244, ('SECS', 'rv'): 3.3520217872798061e-05, ('Mpi',
'SECG'): 0.002219969458552852, ('MC', 'C'): -0.036258588372149528, ('THIS',
'SECS'): -9.2666461860067169e-05, ('tbv', 'bv'): -0.00033210010411835868, ('C',
'C'): 0.36952914866155484, ('SECS', 'trv'): 3.4042952793082553e-06, ('rpi',
'trv'): -0.048045475169797416, ('THIG', 'rpi'): 0.0036935618641351651, ('trv',
'MC'): 0.0010161713431604599, ('rpi', 'rpi'): 5.5622375900778049, ('Mpi', 'C'):
0.033533953127871163, ('MC', 'MC'): 0.0066379815677804297, ('rv', 'trv'):
0.004166816278100764, ('MC', 'THIG'): 0.00020618032904151748, ('M02S', 'M02S'):
0.00060817238394075682, ('trv', 'Mpi'): 0.015250653525779433, ('THIG', 'THIS'):
9.9533010440212809e-06, ('tbv', 'MC'): 0.0061058836885243732, ('tMv', 'C'):
0.0076173475635255179, ('tbv', 'SECS'): 7.6333160918176288e-05, ('SECG',
'M02S'): -0.00065274600517088483, ('trv', 'tMv'): -0.017112636595746616, ('C',
'SECG'): 0.0062260289244321773, ('rpi', 'bv'): -0.00035898280427913557, ('tbv',
'C'): -0.047093749916571746, ('rpi', 'THIG'): 0.0036935618641351651, ('SECS',
'tbv'): 7.6333160918176288e-05, ('M02S', 'rpi'): -0.00083139401293003622,
('tMv', 'SECG'): 0.00057355472202981147, ('THIS', 'trv'):
-8.059915186543189e-06, ('THIG', 'Mpi'): -0.0012774339177200611, ('M02S',
'tMv'): -1.6647261533993509e-05, ('SECS', 'Mv'): 1.7414693760166068e-05, ('MC',
'rpi'): -0.084847904114066064, ('Mpi', 'SECS'): -4.9995881049174269e-05, ('rv',
'rv'): 0.0045218537938113276, ('THIS', 'SECG'): 1.8170497729357869e-05,
('SECS', 'SECG'): -6.9724574897603845e-05, ('THIG', 'MC'):
0.00020618032904151748, ('THIS', 'bv'): 3.3483224714976356e-08, ('Mv', 'trv'):
0.0015379082074404202, ('Mv', 'M02S'): -5.2515916348167587e-05, ('rv', 'SECG'):
-0.00036243927535111489, ('rv', 'Mv'): -0.0019548263784784792, ('tMv', 'tMv'):
0.01917967709116021, ('THIS', 'THIG'): 9.9533010440212809e-06, ('MC', 'M02S'):
0.00012825992107221143, ('rpi', 'C'): 0.08204047012647972, ('M02S', 'THIS'):
-1.2750854559960807e-05, ('bv', 'trv'): -6.45996002760603e-05, ('C', 'rpi'):
0.08204047012647972, ('SECG', 'SECG'): 0.012916128216219243, ('tMv', 'Mpi'):
0.023924151157681827, ('Mpi', 'tbv'): -0.0094233508455311803, ('trv', 'THIG'):
1.9259745629808435e-05, ('C', 'Mv'): -0.00054377989359903124, ('Mv', 'SECG'):
-0.00087216225032152686, ('C', 'SECS'): -1.4277914986737694e-05, ('trv', 'Mv'):
0.0015379082074404202, ('M02S', 'MC'): 0.00012825992107221143, ('rv', 'M02S'):
1.6966102046237698e-05, ('Mv', 'rv'): -0.0019548263784784792, ('SECS', 'tMv'):
1.0454724955636238e-05, ('SECG', 'SECS'): -6.9724574897603845e-05, ('SECS',
'M02S'): -5.0430389084205586e-05, ('M02S', 'trv'): 9.9243841643798775e-05,
('M02S', 'Mv'): -5.2515916348167587e-05, ('MC', 'Mpi'): 0.016774499369639435,
('THIS', 'tMv'): 4.578506313793805e-07, ('C', 'THIG'): -0.0037465604426767055,
('rv', 'SECS'): 3.3520217872798061e-05, ('THIS', 'rv'):
-9.9195203113994761e-06, ('SECS', 'MC'): -1.2520130272335073e-05, ('SECS',
'THIG'): -5.4630238622194277e-06, ('rv', 'tbv'): 0.016138173428212971, ('trv',
'THIS'): -8.059915186543189e-06, ('rv', 'THIG'): 0.00010520521114085837, ('bv',
'rv'): 2.635554115731397e-07, ('M02S', 'tbv'): 0.00012051361385844407, ('SECG',
'tMv'): 0.00057355472202981147, ('bv', 'SECG'): 5.5867312893547216e-07, ('Mpi',
'M02S'): 0.00022850003340474631, ('rpi', 'SECG'): -0.0059901485556498749,
('bv', 'rpi'): -0.00035898280427913557, ('tMv', 'tbv'): 0.018425480149743462,
('THIS', 'Mpi'): 1.9452061103007188e-05}

herm = KM10.copy()
herm.update({
   'Mv' : 4,
   'rv' : 0.740366,
   'bv' : 0.403212,
   'MC' : 0.400134,
  'tMv' : 3.99618,
  'trv' : 0.309523})

lowHt = KM10.copy()
lowHt.update({
   'Mv' : 1.2384,
   'rv' : 0.818009,
   'bv' : 0.400002,
    'C' : 3.06796,
   'MC' : 1.68725,
  'tMv' : 3.99497,
  'trv' : 2,
  'tbv' : 0.400119,
  'rpi' : 3.78526,
  'Mpi' : 3.95165 })

zeroHt = KM10.copy()
zeroHt.update({
   'Mv' : 4.0,
   'rv' : 1.12062,
   'bv' : 0.400009,
    'C' : 0.41788,
   'MC' : 3.96858,
  'trv' : 0,
  'rpi' : 2.93138,
  'Mpi' : 4
})


#######    DM12  ##########

#  P(chi-square, d.o.f) = P(1887.05, 1687) = 0.0004

DM12 = {"NS" : 0.152, "AL0S" : 1.1575, "ALPS" : 0.1028, "M02S" : 0.1814, "DELM2S" : 0., 
"PS" : 2., "SECS" : 0.5142, "THIS" : -0.2106, "SKEWS" : 0., "AL0G" : 1.2473, 
"ALPG" : -0.023, "M02G" : 0.2272, "DELM2G" : 0., "PG" : 2., "SECG" : -4.587, "THIG" : 1.7607, 
"SKEWG" : 0., "KAPS" : 1.5265, "EAL0S" : 1.1754, "EALPS" : 0.0239, "EM02S" : 0.1808,
"EDELM2S" : 0., "EPS" : 2., "ESECS" : 0.5235, "ETHIS" : -0.2194, "ESKEWS" : 0., 
"EAL0G" : 1.8486, "EALPG" : 0.05, "EM02G" : 0.1487, "EDELM2G" : 0., "EPG" : 2., 
"ESECG" : -4.6366, "ETHIG" : 1.8638, "ESKEWG" : 0.,
'fix_ALPS' : False, 'fix_M02S' : False, 'fix_SECS' : False, 'fix_THIS' : False, 'fix_KAPS' : False,
'fix_ALPG' : False, 'fix_M02G' : False, 'fix_SECG' : False, 'fix_THIG' : False, 'fix_EAL0S' : False,
'fix_EALPS' : False, 'fix_ESECS' : False, 'fix_ETHIS' : False, 'fix_EAL0G' : False, 'fix_ESECG' : False}


#     Estimate    Standard Error    t Statistic    P-Value
#  ALPS    0.1028    0.0110955    9.26852    5.55507*10^-20
#  M02S    0.1814    0.00382209    47.4512    3.230049385373141*10^-313
#  SECS    0.5142    0.0532896    9.64923    1.75017*10^-21
#  THIS    -0.2106    0.0150207    -14.0196    2.50685*10^-42
#  KAPS    1.526    0.18243    8.36746    1.21567*10^-16
#  ALPG    -0.02300    0.0718885    -0.319978    0.749024
#  M02G    0.2272    0.0234185    9.70242    1.06926*10^-21
#  SECG    -4.587    0.212031    -21.6339    7.82618*10^-92
#  THIG    1.761    0.104934    16.7794    1.4333*10^-58
#  EAL0S    1.175    0.0152958    76.8444    3.489664919791352*10^-554
#  EALPS    0.02395    0.017757    1.34866    0.177628
#  EM02S    0.1808    0.00464265    38.9539    1.74197*10^-237
#  ESECS    0.5235    0.0919149    5.69577    1.44549*10^-8
#  ETHIS    -0.2194    0.0247518    -8.86231    1.9407*10^-18
#  EAL0G    1.849    0.0399845    46.2329    2.01468*10^-302
#  EM02G    0.1487    0.350896    0.423722    0.671822
#  ESECG    -4.637    0.00502219    -923.221    3.081045734257441*10^-2289


#######    EIC12A  ##########
#
# P(chi-square, d.o.f) = P(1746.46, 1685) = 0.1451
#
#    ALPS =    0.083 +- 0.007
#    M02S =    0.177 +- 0.002
#    SECS =    0.515 +- 0.003
#    THIS =   -0.211 +- 0.001
#    KAPS =    1.327 +- 0.058
#    ALPG =   -0.020 +- 0.059
#    M02G =    0.199 +- 0.011
#    SECG =   -4.562 +- 0.037
#    THIG =    1.750 +- 0.018
#   EAL0S =    1.176 +- 0.006
#   EALPS =    0.056 +- 0.008
#   EM02S =    0.192 +- 0.003
#   ESECS =    0.509 +- 0.009
#   ETHIS =   -0.210 +- 0.003
#   EAL0G =    1.360 +- 0.004
#   EM02G =    0.102 +- 0.023
#   ETHIG =    1.753 +- 0.079
#
EIC12A = {'EAL0G': 1.3597809599163202, 'KAPS': 1.327224163631622, 'fix_ALPG':
False, 'EDELM2S': 0.0, 'EPS': 2.0, 'EAL0S': 1.1763542985774482, 'KAPG': 0.0,
'fix_ETHIG': False, 'EPG': 2.0, 'EDELM2G': 0.0, 'SECG': -4.5621938107115172,
'EKAPG': 0.0, 'ESKEWG': 0.0, 'fix_KAPG': True, 'M02S': 0.17704540271398184,
'NG': 0.5, 'fix_KAPS': False, 'SECS': 0.51482258400203484, 'EKAPS': 0.0, 'NS':
0.152, 'fix_PG': True, 'ALPS': 0.083014840550278871, 'fix_EM02G': False,
'fix_DELM2S': True, 'SKEWG': 0.0, 'fix_EM02S': False, 'SKEWS': 0.0, 'fix_THIS':
False, 'fix_ESKEWS': True, 'fix_ETHIS': False, 'fix_NG': True, 'THIS':
-0.21074242644214614, 'fix_EDELM2S': True, 'fix_NS': True, 'fix_THIG': False,
'fix_DELB': True, 'THIG': 1.7500648350951769, 'fix_EDELM2G': True, 'fix_EPG':
True, 'fix_ALPS': False, 'ETHIG': 1.7528113682623814, 'ESECS':
0.50878056876778166, 'fix_EAL0G': False, 'fix_SECG': False, 'ETHIS':
-0.21010325491275958, 'ESECG': -4.8055000000000003, 'fix_EPS': True,
'fix_EAL0S': False, 'fix_SECS': False, 'PS': 2.0, 'EALPG':
0.050000000000000003, 'EM02S': 0.19166116334516523, 'fix_AL0S': True,
'fix_EALPG': True, 'fix_M02G': False, 'EALPS': 0.056182226959520791, 'ESKEWS':
0.0, 'EM02G': 0.10171629960468281, 'PG': 2.0, 'fix_EALPS': False, 'limit_M02G':
(0.10000000000000001, 1.5), 'fix_M02S': False, 'fix_DELM2G': True, 'fix_AL0G':
True, 'AL0S': 1.1575, 'DELB': 0.0, 'fix_ESECS': False, 'fix_PS': True,
'limit_EM02S': (0.10000000000000001, 1.5), 'fix_ESECG': True, 'AL0G':
1.2473000000000001, 'limit_EM02G': (0.10000000000000001, 1.5), 'limit_M02S':
(0.10000000000000001, 1.5), 'M02G': 0.19924174280855703, 'fix_SKEWG': True,
'DELM2S': 0.0, 'fix_ESKEWG': True, 'ALPG': -0.019731836326246754, 'fix_SKEWS':
True, 'DELM2G': 0.0}

EIC12Acov = {('M02S', 'ALPS'): 1.2049728351022623e-05, ('M02G', 'M02S'):
-2.9565281733383333e-06, ('M02S', 'EM02S'): -6.8332107408400988e-07, ('M02G',
'SECS'): 1.9033698235126133e-06, ('ETHIS', 'ESECS'): -1.7411549640212335e-05,
('KAPS', 'M02G'): 1.6464900790739317e-06, ('ALPG', 'ALPG'):
0.0035116111387995799, ('ETHIS', 'THIS'): -2.4640020190640667e-09, ('EM02S',
'SECS'): 1.8078129998064364e-07, ('ESECS', 'ETHIG'): -0.00010211361516690027,
('M02S', 'EAL0G'): -1.8514050217649475e-09, ('ESECS', 'ALPG'):
1.8490173980185694e-05, ('EM02G', 'EM02G'): 0.00051672237229165514, ('M02S',
'SECG'): -4.745158805249578e-07, ('KAPS', 'ETHIS'): -2.5180410037418932e-05,
('KAPS', 'M02S'): -9.4843054364071431e-06, ('SECS', 'THIS'):
-2.4306310957311147e-06, ('THIS', 'ETHIS'): -2.4640020190640667e-09, ('ETHIS',
'M02S'): -3.3290869189852085e-07, ('M02G', 'EAL0S'): -7.313835738865201e-06,
('ETHIS', 'KAPS'): -2.5180410037418932e-05, ('EM02G', 'ALPS'):
-1.4057044485689222e-06, ('EALPS', 'ESECS'): -1.5361193565530699e-06, ('EALPS',
'EAL0S'): 2.2466520081652913e-05, ('EALPS', 'ETHIS'): -1.9152788339573876e-07,
('EM02S', 'ALPS'): -7.2792728334303704e-07, ('SECG', 'ESECS'):
5.9318383636448337e-07, ('ETHIG', 'KAPS'): 0.0019457609831708706, ('EAL0G',
'ETHIG'): 1.2670313408270592e-05, ('ETHIG', 'THIS'): -8.7255221522515036e-07,
('THIS', 'ALPG'): -4.7323487391975163e-07, ('EAL0S', 'THIG'):
-4.4538707734219515e-07, ('SECS', 'SECS'): 8.2369868661379153e-06, ('ETHIG',
'ALPS'): 7.7943614833933693e-06, ('ALPG', 'EAL0S'): -3.6638828911868013e-06,
('KAPS', 'EM02S'): -5.849900532357199e-05, ('EALPS', 'KAPS'):
-7.8818922805272128e-05, ('THIS', 'M02S'): -2.5384644758036288e-08, ('ESECS',
'EAL0G'): -3.3365901151833993e-09, ('EAL0S', 'THIS'): -1.4802679856148335e-09,
('ESECS', 'ETHIS'): -1.7411549640212335e-05, ('M02G', 'M02G'):
0.00011147622812902296, ('THIG', 'ESECS'): 2.3640209772637901e-06, ('SECS',
'EM02G'): 1.8517588878172949e-07, ('SECG', 'M02G'): -5.0235857310916075e-05,
('EAL0G', 'EAL0S'): 7.9024413845848361e-09, ('SECS', 'ETHIS'):
-4.0658996100351956e-08, ('EAL0S', 'EALPS'): 2.2466520081652913e-05, ('M02S',
'THIG'): -9.2955974671248785e-07, ('EM02G', 'SECS'): 1.8517588878172949e-07,
('EM02G', 'ETHIG'): 0.0002256103891042714, ('EALPS', 'THIG'):
3.3823184782853714e-06, ('THIS', 'ETHIG'): -8.7255221522515036e-07, ('EAL0G',
'THIS'): -8.4252454083836033e-10, ('M02G', 'EAL0G'): -2.7617780434527381e-08,
('THIG', 'ETHIG'): 1.8128617497470098e-05, ('SECG', 'THIS'):
5.4002275631286939e-07, ('ETHIS', 'SECG'): 1.5742845925344904e-07, ('EALPS',
'EAL0G'): 4.0193568936800894e-08, ('ETHIS', 'SECS'): -4.0658996100351956e-08,
('EM02G', 'THIG'): -1.4884490666510246e-06, ('EALPS', 'THIS'):
7.6152998254088232e-08, ('SECG', 'EALPS'): -1.9158732399360614e-07, ('EM02S',
'THIG'): 4.744168825993999e-07, ('ALPS', 'M02G'): 1.7096505001995378e-06,
('ALPS', 'EAL0S'): -1.0122760501659253e-06, ('SECG', 'EM02G'):
-2.1865021836215072e-06, ('EAL0G', 'EM02G'): -9.3068659501089688e-08, ('EALPS',
'M02S'): -1.0549415005855865e-06, ('EM02S', 'THIS'): 5.3502365745368942e-08,
('SECS', 'ALPS'): 2.4723711991996906e-07, ('SECG', 'ALPS'):
-2.6154101761916735e-07, ('ALPS', 'ALPS'): 4.4277971390091834e-05, ('EALPS',
'EM02S'): 1.8251169649282824e-05, ('KAPS', 'SECG'): -1.0137689075722377e-05,
('ALPG', 'THIS'): -4.7323487391975163e-07, ('THIG', 'EM02S'):
4.744168825993999e-07, ('EM02G', 'ESECS'): 9.3151865307933914e-06, ('SECS',
'KAPS'): 2.0436430564805467e-07, ('EAL0S', 'ETHIG'): 0.00013325785093884094,
('M02S', 'M02G'): -2.9565281733383333e-06, ('THIG', 'THIG'):
0.00033146228764424401, ('ALPS', 'ETHIG'): 7.7943614833933693e-06, ('EAL0S',
'EM02G'): -6.2756025716729734e-07, ('M02S', 'THIS'): -2.5384644758036288e-08,
('EM02S', 'M02G'): -1.1499695710174902e-06, ('ALPG', 'ETHIG'):
-0.00014555641578088338, ('M02G', 'ETHIS'): 1.1376294641247202e-06, ('EAL0S',
'KAPS'): -0.00014000591671517752, ('EAL0G', 'EM02S'): -2.4572994676998672e-08,
('EM02G', 'SECG'): -2.1865021836215072e-06, ('ALPG', 'EAL0G'):
-4.1470793415819951e-07, ('M02S', 'EALPS'): -1.0549415005855865e-06, ('KAPS',
'ESECS'): -0.00013736875338772595, ('M02G', 'SECG'): -5.0235857310916075e-05,
('EM02G', 'ETHIS'): 3.9375283292099267e-06, ('EM02S', 'EM02S'):
7.4748739010752092e-06, ('KAPS', 'EM02G'): -0.00020263025042151941, ('ALPG',
'ALPS'): -0.00011574299843498496, ('ETHIS', 'ALPG'): 6.892547178279151e-06,
('ALPS', 'ESECS'): -3.6494408320709715e-07, ('SECS', 'ESECS'):
-1.1798822636329908e-07, ('THIS', 'EM02G'): 4.4951481696254473e-08, ('M02S',
'SECS'): -7.245521493290007e-08, ('ALPG', 'EALPS'): 9.4879291882298886e-05,
('EAL0S', 'SECG'): 8.7156634838624941e-08, ('EM02S', 'ETHIG'):
-9.380520424250625e-06, ('THIG', 'SECG'): -0.00065872317628020793, ('SECG',
'SECS'): -1.9043456547186915e-06, ('EALPS', 'M02G'): -3.3211362171056837e-06,
('ETHIG', 'THIG'): 1.8128617497470098e-05, ('EAL0G', 'M02G'):
-2.7617780434527381e-08, ('ALPS', 'ALPG'): -0.00011574299843498496, ('SECS',
'M02G'): 1.9033698235126133e-06, ('EM02G', 'ALPG'): -5.5011545742614754e-05,
('ETHIG', 'SECS'): -2.3836193423509363e-06, ('EAL0S', 'EAL0G'):
7.9024413845848361e-09, ('ALPS', 'M02S'): 1.2049728351022623e-05, ('EALPS',
'ALPS'): 1.4905561138850678e-06, ('ALPG', 'EM02S'): 2.1981840788377481e-05,
('ALPG', 'M02S'): -1.7205340112925764e-05, ('KAPS', 'THIG'):
-2.6966199493273926e-05, ('EAL0G', 'ETHIS'): -6.8743692935029391e-09, ('ALPG',
'ETHIS'): 6.892547178279151e-06, ('M02G', 'ETHIG'): -5.799064274595644e-05,
('ETHIS', 'EALPS'): -1.9152788339573876e-07, ('M02G', 'THIS'):
1.1257432436600145e-07, ('EM02S', 'SECG'): 9.1520062515918034e-07, ('M02S',
'ETHIS'): -3.3290869189852085e-07, ('THIG', 'EAL0S'): -4.4538707734219515e-07,
('EM02S', 'KAPS'): -5.849900532357199e-05, ('ALPS', 'KAPS'):
7.2149292466479026e-07, ('EALPS', 'EM02G'): 6.456761624656562e-07, ('KAPS',
'EAL0G'): 1.1530158789106907e-06, ('SECS', 'EALPS'): 2.0134196359805206e-07,
('ALPG', 'KAPS'): -0.00028910986837330199, ('THIS', 'THIS'):
7.4899095215078674e-07, ('EM02S', 'EALPS'): 1.8251169649282824e-05, ('EM02G',
'M02G'): -1.278625983885186e-05, ('ETHIS', 'M02G'): 1.1376294641247202e-06,
('THIS', 'M02G'): 1.1257432436600145e-07, ('EALPS', 'SECS'):
2.0134196359805206e-07, ('ALPG', 'ESECS'): 1.8490173980185694e-05, ('THIG',
'SECS'): -1.0985319644281822e-06, ('EM02S', 'M02S'): -6.8332107408400988e-07,
('EAL0S', 'ETHIS'): -2.8363080669008082e-06, ('EM02G', 'EM02S'):
6.0137580340902593e-07, ('ESECS', 'THIS'): -1.0557792663767025e-08, ('THIG',
'M02S'): -9.2955974671248785e-07, ('SECG', 'THIG'): -0.00065872317628020793,
('SECS', 'EAL0G'): -2.1729999762342794e-09, ('M02S', 'KAPS'):
-9.4843054364071431e-06, ('ETHIG', 'SECG'): -2.9603459919588895e-06, ('SECS',
'EM02S'): 1.8078129998064364e-07, ('ESECS', 'M02S'): -1.1230524750169296e-06,
('ALPS', 'EAL0G'): 3.5016994373529228e-08, ('ETHIG', 'EM02G'):
0.0002256103891042714, ('THIG', 'M02G'): 1.3050488470452282e-06, ('ALPS',
'THIS'): 9.1903381561329149e-08, ('KAPS', 'THIS'): -1.555897252700987e-07,
('SECG', 'ETHIS'): 1.5742845925344904e-07, ('ALPG', 'SECG'):
-0.00028553151777308276, ('THIG', 'ETHIS'): 9.4680678095343999e-07, ('ESECS',
'SECS'): -1.1798822636329908e-07, ('EM02S', 'ESECS'): -1.052999341306095e-06,
('ALPS', 'EM02S'): -7.2792728334303704e-07, ('EM02G', 'THIS'):
4.4951481696254473e-08, ('ESECS', 'ESECS'): 8.4478203104939876e-05, ('EM02G',
'KAPS'): -0.00020263025042151941, ('M02S', 'ETHIG'): 7.077377601802516e-06,
('KAPS', 'SECS'): 2.0436430564805467e-07, ('THIS', 'SECS'):
-2.4306310957311147e-06, ('ETHIG', 'EAL0S'): 0.00013325785093884094, ('M02S',
'EM02G'): -9.2048452441504519e-07, ('M02G', 'EALPS'): -3.3211362171056837e-06,
('ALPG', 'M02G'): 0.000230947076082261, ('KAPS', 'ALPG'):
-0.00028910986837330199, ('EALPS', 'SECG'): -1.9158732399360614e-07, ('EALPS',
'EALPS'): 5.7297070250105663e-05, ('SECG', 'KAPS'): -1.0137689075722377e-05,
('EAL0G', 'SECG'): 1.3324607953161108e-07, ('THIG', 'EALPS'):
3.3823184782853714e-06, ('M02S', 'M02S'): 4.3913116774174167e-06, ('EAL0S',
'SECS'): 3.542236718998718e-08, ('ALPS', 'SECS'): 2.4723711991996906e-07,
('SECG', 'EM02S'): 9.1520062515918034e-07, ('THIG', 'THIS'):
4.9203644533953894e-08, ('EAL0S', 'ALPG'): -3.6638828911868013e-06, ('ESECS',
'EM02S'): -1.052999341306095e-06, ('EAL0G', 'THIG'): 2.9915638641773397e-08,
('ETHIG', 'M02S'): 7.077377601802516e-06, ('ALPS', 'THIG'):
-3.0658431055579528e-06, ('SECG', 'M02S'): -4.745158805249578e-07, ('SECG',
'ETHIG'): -2.9603459919588895e-06, ('EAL0G', 'ESECS'): -3.3365901151833993e-09,
('ESECS', 'SECG'): 5.9318383636448337e-07, ('THIG', 'EM02G'):
-1.4884490666510246e-06, ('EAL0G', 'M02S'): -1.8514050217649475e-09, ('EALPS',
'ETHIG'): 7.7515922531079596e-05, ('ESECS', 'ALPS'): -3.6494408320709715e-07,
('SECS', 'ALPG'): -1.437759980737935e-06, ('KAPS', 'ETHIG'):
0.0019457609831708706, ('ALPS', 'EALPS'): 1.4905561138850678e-06, ('KAPS',
'EALPS'): -7.8818922805272128e-05, ('THIG', 'ALPS'): -3.0658431055579528e-06,
('EM02G', 'EAL0S'): -6.2756025716729734e-07, ('EALPS', 'ALPG'):
9.4879291882298886e-05, ('THIS', 'SECG'): 5.4002275631286939e-07, ('ETHIS',
'EM02G'): 3.9375283292099267e-06, ('THIS', 'EALPS'): 7.6152998254088232e-08,
('EAL0G', 'EAL0G'): 1.830643554523766e-05, ('EM02G', 'M02S'):
-9.2048452441504519e-07, ('EAL0G', 'KAPS'): 1.1530158789106907e-06, ('ETHIS',
'ETHIS'): 7.0531804902457707e-06, ('KAPS', 'ALPS'): 7.2149292466479026e-07,
('SECG', 'ALPG'): -0.00028553151777308276, ('M02G', 'THIG'):
1.3050488470452282e-06, ('THIS', 'ALPS'): 9.1903381561329149e-08, ('KAPS',
'KAPS'): 0.0033181160246602615, ('EM02S', 'EAL0G'): -2.4572994676998672e-08,
('EM02S', 'ETHIS'): -1.0099045625278657e-07, ('SECS', 'SECG'):
-1.9043456547186915e-06, ('KAPS', 'EAL0S'): -0.00014000591671517752, ('ETHIG',
'ALPG'): -0.00014555641578088338, ('THIS', 'THIG'): 4.9203644533953894e-08,
('EAL0S', 'M02S'): 2.2455669994562932e-06, ('THIS', 'ESECS'):
-1.0557792663767025e-08, ('SECG', 'SECG'): 0.0013733832412664052, ('THIS',
'EAL0S'): -1.4802679856148335e-09, ('SECG', 'EAL0S'): 8.7156634838624941e-08,
('ETHIG', 'M02G'): -5.799064274595644e-05, ('ETHIG', 'ETHIS'):
-4.1821613737621859e-05, ('EM02S', 'ALPG'): 2.1981840788377481e-05, ('ESECS',
'KAPS'): -0.00013736875338772595, ('M02G', 'ALPG'): 0.000230947076082261,
('ESECS', 'THIG'): 2.3640209772637901e-06, ('M02G', 'EM02S'):
-1.1499695710174902e-06, ('SECS', 'ETHIG'): -2.3836193423509363e-06, ('EM02G',
'EALPS'): 6.456761624656562e-07, ('ALPG', 'EM02G'): -5.5011545742614754e-05,
('EAL0S', 'M02G'): -7.313835738865201e-06, ('ETHIG', 'EALPS'):
7.7515922531079596e-05, ('ALPS', 'ETHIS'): -7.859514369008109e-08, ('ETHIS',
'ALPS'): -7.859514369008109e-08, ('THIS', 'KAPS'): -1.555897252700987e-07,
('SECS', 'M02S'): -7.245521493290007e-08, ('EAL0S', 'ESECS'):
-1.0159434888145642e-05, ('ETHIG', 'EAL0G'): 1.2670313408270592e-05, ('THIG',
'KAPS'): -2.6966199493273926e-05, ('M02S', 'ALPG'): -1.7205340112925764e-05,
('ETHIG', 'ETHIG'): 0.0062641875611631776, ('EAL0G', 'ALPS'):
3.5016994373529228e-08, ('M02G', 'EM02G'): -1.278625983885186e-05, ('ESECS',
'EM02G'): 9.3151865307933914e-06, ('M02G', 'KAPS'): 1.6464900790739317e-06,
('EAL0S', 'EM02S'): 9.1710498122126999e-06, ('ETHIS', 'ETHIG'):
-4.1821613737621859e-05, ('EAL0G', 'ALPG'): -4.1470793415819951e-07, ('EM02S',
'EAL0S'): 9.1710498122126999e-06, ('ALPS', 'EM02G'): -1.4057044485689222e-06,
('M02G', 'ALPS'): 1.7096505001995378e-06, ('SECS', 'EAL0S'):
3.542236718998718e-08, ('EAL0S', 'ALPS'): -1.0122760501659253e-06, ('EM02S',
'EM02G'): 6.0137580340902593e-07, ('SECS', 'THIG'): -1.0985319644281822e-06,
('ETHIS', 'EAL0S'): -2.8363080669008082e-06, ('THIS', 'EM02S'):
5.3502365745368942e-08, ('EAL0G', 'EALPS'): 4.0193568936800894e-08, ('ESECS',
'EALPS'): -1.5361193565530699e-06, ('THIS', 'EAL0G'): -8.4252454083836033e-10,
('SECG', 'EAL0G'): 1.3324607953161108e-07, ('THIG', 'EAL0G'):
2.9915638641773397e-08, ('ETHIG', 'ESECS'): -0.00010211361516690027, ('M02S',
'EAL0S'): 2.2455669994562932e-06, ('EAL0G', 'SECS'): -2.1729999762342794e-09,
('ETHIS', 'EAL0G'): -6.8743692935029391e-09, ('M02G', 'ESECS'):
3.0276276208270805e-06, ('ETHIS', 'THIG'): 9.4680678095343999e-07, ('ESECS',
'M02G'): 3.0276276208270805e-06, ('ALPS', 'SECG'): -2.6154101761916735e-07,
('M02S', 'ESECS'): -1.1230524750169296e-06, ('EAL0S', 'EAL0S'):
4.1090996254539085e-05, ('ETHIS', 'EM02S'): -1.0099045625278657e-07, ('THIG',
'ALPG'): 0.00011944083264526104, ('EM02G', 'EAL0G'): -9.3068659501089688e-08,
('ALPG', 'SECS'): -1.437759980737935e-06, ('ESECS', 'EAL0S'):
-1.0159434888145642e-05, ('ETHIG', 'EM02S'): -9.380520424250625e-06, ('ALPG',
'THIG'): 0.00011944083264526104}


#######    EIC12B  ##########
#
# P(chi-square, d.o.f) = P(1752.77, 1685) = 0.1223

#        |       EIC12B               EIC12A               DM12
# -------+-----------------------------------------------------------------
#   ALPS |    0.081 +- 0.007       0.083 +- 0.007    ALPS     0.102 +- 0.011
#   M02S |    0.177 +- 0.002       0.177 +- 0.002    M02S     0.181 +- 0.003
#   SECS |   -0.076 +- 0.003       0.515 +- 0.003    SECS     0.514 +- 0.053
#   THIS |   -0.033 +- 0.001      -0.211 +- 0.001    THIS    -0.210 +- 0.015
#   KAPS |    0.553 +- 0.030       1.327 +- 0.058    KAPS     1.526 +- 0.182
#   ALPG |    0.127 +- 0.062      -0.020 +- 0.059    ALPG    -0.023 +- 0.071
#   M02G |    0.204 +- 0.014       0.199 +- 0.011    M02G     0.227 +- 0.023
#   SECG |   -4.186 +- 0.041      -4.562 +- 0.037    SECG    -4.587 +- 0.212
#   THIG |    1.683 +- 0.020       1.750 +- 0.018    THIG     1.761 +- 0.104
#  EAL0S |    1.166 +- 0.007       1.176 +- 0.006    EAL0S    1.175 +- 0.015
#  EALPS |    0.058 +- 0.008       0.056 +- 0.008    EALPS    0.023 +- 0.017
#  EM02S |    0.193 +- 0.003       0.192 +- 0.003    EM02S    0.180 +- 0.004
#  ESECS |    0.017 +- 0.023       0.509 +- 0.009    ESECS    0.523 +- 0.091
#  ETHIS |    0.006 +- 0.007      -0.210 +- 0.003    ETHIS   -0.219 +- 0.024
#  EAL0G |    1.281 +- 0.034       1.360 +- 0.004    EAL0G    1.849 +- 0.039
#  EM02G |    0.185 +- 0.008       0.102 +- 0.023    EM02G    0.148 +- 0.350
#  ESECG |   -4.808  (fixed)      -4.805  (fixed)    ESECG   -4.637 +- 0.005
#  ETHIG |    0.097 +- 0.377       1.753 +- 0.079    ETHIG    1.8638 (fixed)




EIC12B = {'EAL0G': 1.2805971215351604, 'KAPS': 0.55255578605233191, 'fix_ALPG':
False, 'EDELM2S': 0.0, 'EPS': 2.0, 'EAL0S': 1.1657208144370803, 'KAPG': 0.0,
'fix_ETHIG': False, 'EPG': 2.0, 'EDELM2G': 0.0, 'SECG': -4.1855350243434275,
'EKAPG': 0.0, 'ESKEWG': 0.0, 'fix_KAPG': True, 'M02S': 0.17683972235685824,
'NG': 0.5, 'fix_KAPS': False, 'SECS': -0.07619257067625522, 'EKAPS': 0.0, 'NS':
0.152, 'fix_PG': True, 'ALPS': 0.080506506379079482, 'fix_EM02G': False,
'fix_DELM2S': True, 'SKEWG': 0.0, 'fix_EM02S': False, 'SKEWS': 0.0, 'fix_THIS':
False, 'fix_ESKEWS': True, 'fix_ETHIS': False, 'fix_NG': True, 'THIS':
-0.032658600865168509, 'fix_EDELM2S': True, 'fix_NS': True, 'fix_THIG': False,
'fix_DELB': True, 'THIG': 1.6826549379699176, 'fix_EDELM2G': True, 'fix_EPG':
True, 'fix_ALPS': False, 'ETHIG': 0.09654407713620787, 'ESECS':
0.01742574168690041, 'fix_EAL0G': False, 'fix_SECG': False, 'ETHIS':
0.0056243056996434643, 'ESECG': -4.8082000000000003, 'fix_EPS': True,
'fix_EAL0S': False, 'fix_SECS': False, 'PS': 2.0, 'EALPG':
0.050000000000000003, 'EM02S': 0.19297940466612939, 'fix_AL0S': True,
'fix_EALPG': True, 'fix_M02G': False, 'EALPS': 0.057835219358367523, 'ESKEWS':
0.0, 'EM02G': 0.18467777779536093, 'PG': 2.0, 'fix_EALPS': False, 'limit_M02G':
(0.10000000000000001, 1.5), 'fix_M02S': False, 'fix_DELM2G': True, 'fix_AL0G':
True, 'AL0S': 1.1575, 'DELB': 0.0, 'fix_ESECS': False, 'fix_PS': True,
'limit_EM02S': (0.10000000000000001, 1.5), 'fix_ESECG': True, 'AL0G':
1.2473000000000001, 'limit_EM02G': (0.10000000000000001, 1.5), 'limit_M02S':
(0.10000000000000001, 1.5), 'M02G': 0.20405296690561009, 'fix_SKEWG': True,
'DELM2S': 0.0, 'fix_ESKEWG': True, 'ALPG': 0.12690487080060397, 'fix_SKEWS':
True, 'DELM2G': 0.0}

EIC12Bcov = {('M02S', 'ALPS'): 1.260979878641184e-05, ('M02G', 'M02S'):
-6.5242491086238194e-06, ('M02S', 'EM02S'): -5.1303296996099959e-07, ('M02G',
'SECS'): -6.0218982462590493e-08, ('ETHIS', 'ESECS'): -6.660367820400771e-05,
('KAPS', 'M02G'): -5.2297903143786912e-06, ('ALPG', 'ALPG'):
0.003868572961556415, ('ETHIS', 'THIS'): -1.9843412644694561e-08, ('EM02S',
'SECS'): 2.4616065853973514e-07, ('ESECS', 'ETHIG'): 0.00027752646009619741,
('M02S', 'EAL0G'): -1.7859093569898691e-06, ('ESECS', 'ALPG'):
-1.8825468573201988e-05, ('EM02G', 'EM02G'): 6.534063645925377e-05, ('M02S',
'SECG'): 1.7400489818729947e-07, ('KAPS', 'ETHIS'): -7.6893240273754205e-05,
('KAPS', 'M02S'): -3.6385344412308079e-06, ('SECS', 'THIS'):
-2.0386199628974772e-06, ('THIS', 'ETHIS'): -1.9843412644694561e-08, ('ETHIS',
'M02S'): -1.2209024429379266e-06, ('M02G', 'EAL0S'): -9.2499980098943874e-06,
('ETHIS', 'KAPS'): -7.6893240273754205e-05, ('EM02G', 'ALPS'):
-7.458458518836197e-07, ('EALPS', 'ESECS'): -8.6746179578284716e-06, ('EALPS',
'EAL0S'): 2.2669496445206561e-05, ('EALPS', 'ETHIS'): -2.6556059728779129e-06,
('EM02S', 'ALPS'): 6.467644794030242e-07, ('SECG', 'ESECS'):
2.6442913256888534e-06, ('ETHIG', 'KAPS'): 0.0033328686666509, ('EAL0G',
'ETHIG'): 0.010946422734973156, ('ETHIG', 'THIS'): -1.6342539646336474e-06,
('THIS', 'ALPG'): -4.2760808896701591e-07, ('EAL0S', 'THIG'):
3.7276510139084841e-07, ('SECS', 'SECS'): 6.9551067026797656e-06, ('ETHIG',
'ALPS'): 0.00010491700145781693, ('ALPG', 'EAL0S'): 4.1998646095354774e-06,
('KAPS', 'EM02S'): -9.7389873958406264e-06, ('EALPS', 'KAPS'):
-1.1173150291810436e-05, ('THIS', 'M02S'): -3.3426110533137095e-08, ('ESECS',
'EAL0G'): 8.8021887011745172e-05, ('EAL0S', 'THIS'): 8.3871753797471618e-10,
('ESECS', 'ETHIS'): -6.660367820400771e-05, ('M02G', 'M02G'):
0.0001959032717565426, ('THIG', 'ESECS'): 1.92440310355236e-06, ('SECS',
'EM02G'): -2.9855704706382502e-07, ('SECG', 'M02G'): -5.8617946586914997e-05,
('EAL0G', 'EAL0S'): -8.0872135851053057e-05, ('SECS', 'ETHIS'):
-6.4759200375389983e-08, ('EAL0S', 'EALPS'): 2.2669496445206561e-05, ('M02S',
'THIG'): -6.4644968205716274e-07, ('EM02G', 'SECS'): -2.9855704706382502e-07,
('EM02G', 'ETHIG'): 0.0014453855447171162, ('EALPS', 'THIG'):
1.7368744677116121e-06, ('THIS', 'ETHIG'): -1.6342539646336474e-06, ('EAL0G',
'THIS'): 1.2107471905007029e-07, ('M02G', 'EAL0G'): 3.859316213756672e-05,
('THIG', 'ETHIG'): 0.00015540424646918536, ('SECG', 'THIS'):
4.0575912805264915e-07, ('ETHIS', 'SECG'): 1.034437963665292e-06, ('EALPS',
'EAL0G'): 3.5285256752309806e-06, ('ETHIS', 'SECS'): -6.4759200375389983e-08,
('EM02G', 'THIG'): 3.0243020500196741e-06, ('EALPS', 'THIS'):
1.0111458140306299e-07, ('SECG', 'EALPS'): 3.8965258006337444e-06, ('EM02S',
'THIG'): -2.9707448919391737e-07, ('ALPS', 'M02G'): -2.57934039316493e-05,
('ALPS', 'EAL0S'): -9.8968307219750694e-07, ('SECG', 'EM02G'):
1.2531138588759759e-05, ('EAL0G', 'EM02G'): 4.6601093421252054e-05, ('EALPS',
'M02S'): -9.0218299403484985e-07, ('EM02S', 'THIS'): 7.9480175783448388e-08,
('SECS', 'ALPS'): 2.6790189734449376e-07, ('SECG', 'ALPS'):
-8.7736264435156026e-07, ('ALPS', 'ALPS'): 4.445920858052714e-05, ('EALPS',
'EM02S'): 1.9265179002736986e-05, ('KAPS', 'SECG'): -2.2565826151084938e-06,
('ALPG', 'THIS'): -4.2760808896701591e-07, ('THIG', 'EM02S'):
-2.9707448919391737e-07, ('EM02G', 'ESECS'): 1.8469179442784673e-05, ('SECS',
'KAPS'): -1.3270562742701384e-07, ('EAL0S', 'ETHIG'): -0.00078018630011961328,
('M02S', 'M02G'): -6.5242491086238194e-06, ('THIG', 'THIG'):
0.00041027689145104926, ('ALPS', 'ETHIG'): 0.00010491700145781693, ('EAL0S',
'EM02G'): -1.2874606312116068e-06, ('M02S', 'THIS'): -3.3426110533137095e-08,
('EM02S', 'M02G'): -2.783581538605611e-07, ('ALPG', 'ETHIG'):
-0.0024510676007550303, ('M02G', 'ETHIS'): 3.128664338984604e-06, ('EAL0S',
'KAPS'): -3.2799494962942556e-05, ('EAL0G', 'EM02S'): -7.8736649776381057e-07,
('EM02G', 'SECG'): 1.2531138588759759e-05, ('ALPG', 'EAL0G'):
-0.00023722724696855613, ('M02S', 'EALPS'): -9.0218299403484985e-07, ('KAPS',
'ESECS'): -0.00025398927479318022, ('M02G', 'SECG'): -5.8617946586914997e-05,
('EM02G', 'ETHIS'): 6.836501170289595e-06, ('EM02S', 'EM02S'):
9.2269136914048872e-06, ('KAPS', 'EM02G'): -6.6415456189136351e-06, ('ALPG',
'ALPS'): -0.00012871478432552314, ('ETHIS', 'ALPG'): -4.7056839835685231e-06,
('ALPS', 'ESECS'): -3.1740012536633365e-07, ('SECS', 'ESECS'):
-2.0121674558731996e-07, ('THIS', 'EM02G'): -1.2647423712460107e-07, ('M02S',
'SECS'): -1.1158288162583633e-07, ('ALPG', 'EALPS'): -8.4073779935221309e-06,
('EAL0S', 'SECG'): 5.3375163828298914e-07, ('EM02S', 'ETHIG'):
-9.410922980137059e-05, ('THIG', 'SECG'): -0.00080522582934712105, ('SECG',
'SECS'): -1.6234446819554406e-06, ('EALPS', 'M02G'): -2.4101094586187489e-06,
('ETHIG', 'THIG'): 0.00015540424646918536, ('EAL0G', 'M02G'):
3.859316213756672e-05, ('ALPS', 'ALPG'): -0.00012871478432552314, ('SECS',
'M02G'): -6.0218982462590493e-08, ('EM02G', 'ALPG'): -6.0493406816738566e-05,
('ETHIG', 'SECS'): -3.715569591008904e-06, ('EAL0S', 'EAL0G'):
-8.0872135851053057e-05, ('ALPS', 'M02S'): 1.260979878641184e-05, ('EALPS',
'ALPS'): 5.471090455228193e-06, ('ALPG', 'EM02S'): -2.9574085356854936e-06,
('ALPG', 'M02S'): -1.5366601770891061e-05, ('KAPS', 'THIG'):
9.1683269202221572e-07, ('EAL0G', 'ETHIS'): 2.6012268288907672e-05, ('ALPG',
'ETHIS'): -4.7056839835685231e-06, ('M02G', 'ETHIG'): 0.00029737699870826921,
('ETHIS', 'EALPS'): -2.6556059728779129e-06, ('M02G', 'THIS'):
1.134385361897796e-07, ('EM02S', 'SECG'): 3.8498005168720357e-07, ('M02S',
'ETHIS'): -1.2209024429379266e-06, ('THIG', 'EAL0S'): 3.7276510139084841e-07,
('EM02S', 'KAPS'): -9.7389873958406264e-06, ('ALPS', 'KAPS'):
-1.2961900529615779e-06, ('EALPS', 'EM02G'): 9.6945089225221332e-06, ('KAPS',
'EAL0G'): -1.4217173450886373e-05, ('SECS', 'EALPS'): 2.9089737545948461e-07,
('ALPG', 'KAPS'): -2.1301776011306715e-05, ('THIS', 'THIS'):
6.3150623482119438e-07, ('EM02S', 'EALPS'): 1.9265179002736986e-05, ('EM02G',
'M02G'): -1.751679451684225e-05, ('ETHIS', 'M02G'): 3.128664338984604e-06,
('THIS', 'M02G'): 1.134385361897796e-07, ('EALPS', 'SECS'):
2.9089737545948461e-07, ('ALPG', 'ESECS'): -1.8825468573201988e-05, ('THIG',
'SECS'): -1.3398426270970084e-06, ('EM02S', 'M02S'): -5.1303296996099959e-07,
('EAL0S', 'ETHIS'): -1.396425222511639e-05, ('EM02G', 'EM02S'):
-4.3946765356726235e-06, ('ESECS', 'THIS'): -6.1176682395655221e-08, ('THIG',
'M02S'): -6.4644968205716274e-07, ('SECG', 'THIG'): -0.00080522582934712105,
('SECS', 'EAL0G'): 3.9313690929845491e-07, ('M02S', 'KAPS'):
-3.6385344412308079e-06, ('ETHIG', 'SECG'): 0.00011474450665029275, ('SECS',
'EM02S'): 2.4616065853973514e-07, ('ESECS', 'M02S'): -3.9214056008025818e-06,
('ALPS', 'EAL0G'): 1.1969287332789032e-05, ('ETHIG', 'EM02G'):
0.0014453855447171162, ('THIG', 'M02G'): 2.9742272198131247e-05, ('ALPS',
'THIS'): 6.8002496537966799e-08, ('KAPS', 'THIS'): -3.3916140165418734e-08,
('SECG', 'ETHIS'): 1.034437963665292e-06, ('ALPG', 'SECG'):
-0.00015435010660870226, ('THIG', 'ETHIS'): 5.9432422799379101e-07, ('ESECS',
'SECS'): -2.0121674558731996e-07, ('EM02S', 'ESECS'): -1.037581320726209e-05,
('ALPS', 'EM02S'): 6.467644794030242e-07, ('EM02G', 'THIS'):
-1.2647423712460107e-07, ('ESECS', 'ESECS'): 0.00052102401192583828, ('EM02G',
'KAPS'): -6.6415456189136351e-06, ('M02S', 'ETHIG'): -4.3258434880276861e-05,
('KAPS', 'SECS'): -1.3270562742701384e-07, ('THIS', 'SECS'):
-2.0386199628974772e-06, ('ETHIG', 'EAL0S'): -0.00078018630011961328, ('M02S',
'EM02G'): -6.4675861975468266e-07, ('M02G', 'EALPS'): -2.4101094586187489e-06,
('ALPG', 'M02G'): 0.00078404998800953734, ('KAPS', 'ALPG'):
-2.1301776011306715e-05, ('EALPS', 'SECG'): 3.8965258006337444e-06, ('EALPS',
'EALPS'): 5.8268505412828015e-05, ('SECG', 'KAPS'): -2.2565826151084938e-06,
('EAL0G', 'SECG'): -4.3791724542021899e-06, ('THIG', 'EALPS'):
1.7368744677116121e-06, ('M02S', 'M02S'): 4.8758794284386603e-06, ('EAL0S',
'SECS'): 8.5026420243886942e-09, ('ALPS', 'SECS'): 2.6790189734449376e-07,
('SECG', 'EM02S'): 3.8498005168720357e-07, ('THIG', 'THIS'):
-1.4326130913442591e-08, ('EAL0S', 'ALPG'): 4.1998646095354774e-06, ('ESECS',
'EM02S'): -1.037581320726209e-05, ('EAL0G', 'THIG'): 1.2972205895813857e-06,
('ETHIG', 'M02S'): -4.3258434880276861e-05, ('ALPS', 'THIG'):
-4.862500965018276e-06, ('SECG', 'M02S'): 1.7400489818729947e-07, ('SECG',
'ETHIG'): 0.00011474450665029275, ('EAL0G', 'ESECS'): 8.8021887011745172e-05,
('ESECS', 'SECG'): 2.6442913256888534e-06, ('THIG', 'EM02G'):
3.0243020500196741e-06, ('EAL0G', 'M02S'): -1.7859093569898691e-06, ('EALPS',
'ETHIG'): 0.00022361377446369356, ('ESECS', 'ALPS'): -3.1740012536633365e-07,
('SECS', 'ALPG'): -3.8493844910861319e-06, ('KAPS', 'ETHIG'):
0.0033328686666509, ('ALPS', 'EALPS'): 5.471090455228193e-06, ('KAPS',
'EALPS'): -1.1173150291810436e-05, ('THIG', 'ALPS'): -4.862500965018276e-06,
('EM02G', 'EAL0S'): -1.2874606312116068e-06, ('EALPS', 'ALPG'):
-8.4073779935221309e-06, ('THIS', 'SECG'): 4.0575912805264915e-07, ('ETHIS',
'EM02G'): 6.836501170289595e-06, ('THIS', 'EALPS'): 1.0111458140306299e-07,
('EAL0G', 'EAL0G'): 0.0011372238711895204, ('EM02G', 'M02S'):
-6.4675861975468266e-07, ('EAL0G', 'KAPS'): -1.4217173450886373e-05, ('ETHIS',
'ETHIS'): 4.801113948643847e-05, ('KAPS', 'ALPS'): -1.2961900529615779e-06,
('SECG', 'ALPG'): -0.00015435010660870226, ('M02G', 'THIG'):
2.9742272198131247e-05, ('THIS', 'ALPS'): 6.8002496537966799e-08, ('KAPS',
'KAPS'): 0.00089973628524091221, ('EM02S', 'EAL0G'): -7.8736649776381057e-07,
('EM02S', 'ETHIS'): -3.1925351651332098e-06, ('SECS', 'SECG'):
-1.6234446819554406e-06, ('KAPS', 'EAL0S'): -3.2799494962942556e-05, ('ETHIG',
'ALPG'): -0.0024510676007550303, ('THIS', 'THIG'): -1.4326130913442591e-08,
('EAL0S', 'M02S'): 2.6208346919670487e-06, ('THIS', 'ESECS'):
-6.1176682395655221e-08, ('SECG', 'SECG'): 0.0016703029400358576, ('THIS',
'EAL0S'): 8.3871753797471618e-10, ('SECG', 'EAL0S'): 5.3375163828298914e-07,
('ETHIG', 'M02G'): 0.00029737699870826921, ('ETHIG', 'ETHIS'):
-1.5796090542597436e-05, ('EM02S', 'ALPG'): -2.9574085356854936e-06, ('ESECS',
'KAPS'): -0.00025398927479318022, ('M02G', 'ALPG'): 0.00078404998800953734,
('ESECS', 'THIG'): 1.92440310355236e-06, ('M02G', 'EM02S'):
-2.783581538605611e-07, ('SECS', 'ETHIG'): -3.715569591008904e-06, ('EM02G',
'EALPS'): 9.6945089225221332e-06, ('ALPG', 'EM02G'): -6.0493406816738566e-05,
('EAL0S', 'M02G'): -9.2499980098943874e-06, ('ETHIG', 'EALPS'):
0.00022361377446369356, ('ALPS', 'ETHIS'): -6.4957744787557963e-08, ('ETHIS',
'ALPS'): -6.4957744787557963e-08, ('THIS', 'KAPS'): -3.3916140165418734e-08,
('SECS', 'M02S'): -1.1158288162583633e-07, ('EAL0S', 'ESECS'):
-4.3578747678479637e-05, ('ETHIG', 'EAL0G'): 0.010946422734973156, ('THIG',
'KAPS'): 9.1683269202221572e-07, ('M02S', 'ALPG'): -1.5366601770891061e-05,
('ETHIG', 'ETHIG'): 0.14211991206545635, ('EAL0G', 'ALPS'):
1.1969287332789032e-05, ('M02G', 'EM02G'): -1.751679451684225e-05, ('ESECS',
'EM02G'): 1.8469179442784673e-05, ('M02G', 'KAPS'): -5.2297903143786912e-06,
('EAL0S', 'EM02S'): 1.1054841115607334e-05, ('ETHIS', 'ETHIG'):
-1.5796090542597436e-05, ('EAL0G', 'ALPG'): -0.00023722724696855613, ('EM02S',
'EAL0S'): 1.1054841115607334e-05, ('ALPS', 'EM02G'): -7.458458518836197e-07,
('M02G', 'ALPS'): -2.57934039316493e-05, ('SECS', 'EAL0S'):
8.5026420243886942e-09, ('EAL0S', 'ALPS'): -9.8968307219750694e-07, ('EM02S',
'EM02G'): -4.3946765356726235e-06, ('SECS', 'THIG'): -1.3398426270970084e-06,
('ETHIS', 'EAL0S'): -1.396425222511639e-05, ('THIS', 'EM02S'):
7.9480175783448388e-08, ('EAL0G', 'EALPS'): 3.5285256752309806e-06, ('ESECS',
'EALPS'): -8.6746179578284716e-06, ('THIS', 'EAL0G'): 1.2107471905007029e-07,
('SECG', 'EAL0G'): -4.3791724542021899e-06, ('THIG', 'EAL0G'):
1.2972205895813857e-06, ('ETHIG', 'ESECS'): 0.00027752646009619741, ('M02S',
'EAL0S'): 2.6208346919670487e-06, ('EAL0G', 'SECS'): 3.9313690929845491e-07,
('ETHIS', 'EAL0G'): 2.6012268288907672e-05, ('M02G', 'ESECS'):
7.7052225983942358e-06, ('ETHIS', 'THIG'): 5.9432422799379101e-07, ('ESECS',
'M02G'): 7.7052225983942358e-06, ('ALPS', 'SECG'): -8.7736264435156026e-07,
('M02S', 'ESECS'): -3.9214056008025818e-06, ('EAL0S', 'EAL0S'):
4.9497423901356871e-05, ('ETHIS', 'EM02S'): -3.1925351651332098e-06, ('THIG',
'ALPG'): 0.00014920217524762942, ('EM02G', 'EAL0G'): 4.6601093421252054e-05,
('ALPG', 'SECS'): -3.8493844910861319e-06, ('ESECS', 'EAL0S'):
-4.3578747678479637e-05, ('ETHIG', 'EM02S'): -9.410922980137059e-05, ('ALPG',
'THIG'): 0.00014920217524762942}
