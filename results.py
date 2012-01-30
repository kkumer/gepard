
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

DM12 = {
"NS" : 0.152, "AL0S" : 1.1575, "ALPS" : 0.1, 
"M02S" : 0.17857142857142858, "DELM2S" : 0., "PS" : 2., 
"SECS" : 0.5127, "THIS" : -0.21, 
"SKEWS" : 0., "AL0G" : 1.2473, "ALPG" : 0.1, "M02G" : 0.25, "DELM2G" : 0., "PG" : 2., 
"SECG" : -4.8055, "THIG" : 1.8638, "SKEWG" : 0., 
"KAPS" : 1.5, "EAL0S" : 1.1575, "EALPS" : 0.02, "EM02S" : 1./6.6, "EDELM2S" : 0., 
"EPS" : 2., "ESECS" : 0.5127, "ETHIS" : -0.21, 
"ESKEWS" : 0., "EAL0G" : 1.2473, "EALPG" : 0.05, "EM02G" : 0.2, "EDELM2G" : 0., "EPG" : 2., 
"ESECG" : -4.8055, "ETHIG" : 1.8638, "ESKEWG" : 0.}

DM12B = {"NS" : 0.152, "AL0S" : 1.1575, "ALPS" : 0.101, "M02S" : 0.1799, 
 "DELM2S" : 0., "PS" : 2., "SECS" : 0.5181, "THIS" : -0.2113, 
 "SKEWS" : 0., "AL0G" : 1.2473, "ALPG" : 0.0727, "M02G" : 0.2401, 
 "DELM2G" : 0., "PG" : 2., "SECG" : -4.6367, "THIG" : 1.7804, 
 "SKEWG" : 0., "KAPS" : 1.3607, "EAL0S" : 1.1726, 
 "EALPS" : 0.0178, "EM02S" : 0.1515, "EDELM2S" : 0., "EPS" : 2., 
 "ESECS" : 0.4588, "ETHIS" : -0.1941, "ESKEWS" : 0., 
 "EAL0G" : 1.3516, "EALPG" : 0.05, "EM02G" : 0.2, "EDELM2G" : 0., 
 "EPG" : 2., "ESECG" : -4.8082, "ETHIG" : 1.8638, "ESKEWG" : 0.}

EIC12 = {'EAL0G': 1.3597814352989468, 'KAPS': 1.4826933110382694, 'fix_ALPG': False,
'EDELM2S': 0.0, 'EPS': 2.0, 'EAL0S': 1.1563710165663585, 'KAPG': 0.0,
'fix_ETHIG': False, 'EPG': 2.0, 'EDELM2G': 0.0, 'SECG': -4.80030098989329,
'EKAPG': 0.0, 'ESKEWG': 0.0, 'fix_KAPG': True, 'M02S': 0.17937862461764434,
'NG': 0.5, 'fix_KAPS': False, 'SECS': 0.51226256188813835, 'EKAPS': 0.0, 'NS':
0.152, 'fix_PG': True, 'ALPS': 0.089373302356791487, 'fix_EM02G': False,
'fix_DELM2S': True, 'SKEWG': 0.0, 'fix_EM02S': False, 'SKEWS': 0.0, 'fix_THIS':
False, 'fix_ESKEWS': True, 'fix_ETHIS': False, 'fix_NG': True, 'THIS':
-0.21016444056375719, 'fix_EDELM2S': True, 'fix_NS': True, 'fix_THIG': False,
'fix_DELB': True, 'THIG': 1.8656572821913338, 'fix_EDELM2G': True, 'fix_EPG':
True, 'fix_ALPS': False, 'ETHIG': 1.796777834019575, 'ESECS':
0.51079535661179543, 'fix_EAL0G': False, 'fix_SECG': False, 'ETHIS':
-0.21059338446198711, 'ESECG': -4.8055000000000003, 'fix_EPS': True,
'fix_EAL0S': False, 'fix_SECS': False, 'PS': 2.0, 'EALPG':
0.050000000000000003, 'EM02S': 0.15151515151515155, 'fix_AL0S': True,
'fix_EALPG': True, 'fix_M02G': False, 'EALPS': 0.025321317404584983, 'ESKEWS':
0.0, 'EM02G': 0.20000000000000004, 'PG': 2.0, 'fix_EALPS': False, 'limit_M02G':
(0.10000000000000001, 1.5), 'fix_M02S': False, 'fix_DELM2G': True, 'fix_AL0G':
True, 'AL0S': 1.1575, 'DELB': 0.0, 'fix_ESECS': False, 'fix_PS': True,
'limit_EM02S': (0.10000000000000001, 1.5), 'fix_ESECG': True, 'AL0G':
1.2473000000000001, 'limit_EM02G': (0.10000000000000001, 1.5), 'limit_M02S':
(0.10000000000000001, 1.5), 'M02G': 0.19693384526258259, 'fix_SKEWG': True,
'DELM2S': 0.0, 'fix_ESKEWG': True, 'ALPG': -0.039805409303853505, 'fix_SKEWS':
True, 'DELM2G': 0.0}

EIC12cov = {('M02S', 'ALPS'): -1.521013508245142e-07, ('M02G', 'M02S'):
-8.1279473828786811e-07, ('M02S', 'EM02S'): 0.0, ('M02G', 'SECS'):
1.3561294171657433e-06, ('ETHIS', 'ESECS'): 1.7491483939595893e-08, ('KAPS',
'M02G'): 1.1738921189923695e-05, ('ALPG', 'ALPG'): 0.00035978721176466784,
('ETHIS', 'THIS'): -6.0797412024120346e-09, ('EM02S', 'SECS'): 0.0, ('ESECS',
'ETHIG'): 5.155730523607306e-06, ('M02S', 'EAL0G'): -3.8621911523573612e-09,
('ESECS', 'ALPG'): 1.1416623328303995e-06, ('EM02G', 'EM02G'):
0.26000000000000006, ('M02S', 'SECG'): 9.9430841783265254e-08, ('KAPS',
'ETHIS'): -1.0470598944252407e-06, ('KAPS', 'M02S'): -4.2413128071918736e-07,
('SECS', 'THIS'): 6.08051419300745e-09, ('THIS', 'ETHIS'):
-6.0797412024120346e-09, ('ETHIS', 'M02S'): -7.9804887592991765e-09, ('M02G',
'EAL0S'): -1.4607629894898336e-06, ('ETHIS', 'KAPS'): -1.0470598944252407e-06,
('EM02G', 'ALPS'): 0.0, ('EALPS', 'ESECS'): 5.0979307037748514e-07, ('EALPS',
'EAL0S'): 1.5084632635774174e-06, ('EALPS', 'ETHIS'): 1.9021730648414777e-07,
('EM02S', 'ALPS'): 0.0, ('SECG', 'ESECS'): 8.1249337700340916e-07, ('ETHIG',
'KAPS'): -4.6025488131970195e-05, ('EAL0G', 'ETHIG'): 8.651261809031964e-07,
('ETHIG', 'THIS'): -1.6710129677247377e-07, ('THIS', 'ALPG'):
-2.0944531693854201e-07, ('EAL0S', 'THIG'): 8.5274493541677246e-07, ('SECS',
'SECS'): 1.368636382520591e-07, ('ETHIG', 'ALPS'): 4.2504182584925699e-06,
('ALPG', 'EAL0S'): 2.7705144957565145e-06, ('KAPS', 'EM02S'): 0.0, ('EALPS',
'KAPS'): -5.0707307468705406e-06, ('THIS', 'M02S'): -7.2779696000472363e-09,
('ESECS', 'EAL0G'): 1.6120423127841767e-08, ('EAL0S', 'THIS'):
-4.5441474101194452e-08, ('ESECS', 'ETHIS'): 1.7491483939595893e-08, ('M02G',
'M02G'): 8.3207258240174245e-05, ('THIG', 'ESECS'): 2.3150285229816843e-07,
('SECS', 'EM02G'): 0.0, ('SECG', 'M02G'): -2.0453950035592996e-05, ('EAL0G',
'EAL0S'): 7.0076954873384941e-08, ('SECS', 'ETHIS'): -1.6331320322652216e-08,
('EAL0S', 'EALPS'): 1.5084632635774174e-06, ('M02S', 'THIG'):
3.1279039374527087e-08, ('EM02G', 'SECS'): 0.0, ('EM02G', 'ETHIG'): 0.0,
('EALPS', 'THIG'): 1.3088073710247217e-06, ('THIS', 'ETHIG'):
-1.6710129677247377e-07, ('EAL0G', 'THIS'): 6.8966289500986728e-10, ('M02G',
'EAL0G'): 3.9196135510126979e-07, ('THIG', 'ETHIG'): 4.4714252515195068e-06,
('SECG', 'THIS'): -2.2622796924209372e-07, ('ETHIS', 'SECG'):
3.2279231295898578e-07, ('EALPS', 'EAL0G'): 2.5535659974614765e-08, ('ETHIS',
'SECS'): -1.6331320322652216e-08, ('EM02G', 'THIG'): 0.0, ('EALPS', 'THIS'):
-8.5440008282716879e-08, ('SECG', 'EALPS'): 3.8961440993238121e-06, ('EM02S',
'THIG'): 0.0, ('ALPS', 'M02G'): 1.4359224420333175e-05, ('ALPS', 'EAL0S'):
-4.800957442453362e-07, ('SECG', 'EM02G'): 0.0, ('EAL0G', 'EM02G'): 0.0,
('EALPS', 'M02S'): 1.2221346658311379e-07, ('EM02S', 'THIS'): 0.0, ('SECS',
'ALPS'): 3.4950636427360906e-07, ('SECG', 'ALPS'): -3.014844665917385e-06,
('ALPS', 'ALPS'): 4.546476481332087e-06, ('EALPS', 'EM02S'): 0.0, ('KAPS',
'SECG'): -1.3388869764770864e-05, ('ALPG', 'THIS'): -2.0944531693854201e-07,
('THIG', 'EM02S'): 0.0, ('EM02G', 'ESECS'): 0.0, ('SECS', 'KAPS'):
4.8176197414287816e-07, ('EAL0S', 'ETHIG'): 1.4808752765086284e-05, ('M02S',
'M02G'): -8.1279473828786811e-07, ('THIG', 'THIG'): 3.9809323782685245e-06,
('ALPS', 'ETHIG'): 4.2504182584925699e-06, ('EAL0S', 'EM02G'): 0.0, ('M02S',
'THIS'): -7.2779696000472363e-09, ('EM02S', 'M02G'): 0.0, ('ALPG', 'ETHIG'):
4.061777554572628e-05, ('M02G', 'ETHIS'): -1.5645339041813079e-07, ('EAL0S',
'KAPS'): -7.9008316901597777e-06, ('EAL0G', 'EM02S'): 0.0, ('EM02G', 'SECG'):
0.0, ('ALPG', 'EAL0G'): -1.7857710566950838e-07, ('M02S', 'EALPS'):
1.2221346658311379e-07, ('KAPS', 'ESECS'): -2.8598706302377214e-06, ('M02G',
'SECG'): -2.0453950035592996e-05, ('EM02G', 'ETHIS'): 0.0, ('EM02S', 'EM02S'):
0.13893480257116625, ('KAPS', 'EM02G'): 0.0, ('ALPG', 'ALPS'):
-7.1972114718848861e-06, ('ETHIS', 'ALPG'): 4.3707938936467013e-07, ('ALPS',
'ESECS'): -2.3965071820922081e-08, ('SECS', 'ESECS'): -4.3603850879636102e-08,
('THIS', 'EM02G'): 0.0, ('M02S', 'SECS'): -2.1587221004910876e-08, ('ALPG',
'EALPS'): 7.0052840588288252e-06, ('EAL0S', 'SECG'): 2.9681814636767001e-06,
('EM02S', 'ETHIG'): 0.0, ('THIG', 'SECG'): 2.625687795606487e-06, ('SECG',
'SECS'): -6.3946744943539416e-07, ('EALPS', 'M02G'): -8.5952224552298213e-06,
('ETHIG', 'THIG'): 4.4714252515195068e-06, ('EAL0G', 'M02G'):
3.9196135510126979e-07, ('ALPS', 'ALPG'): -7.1972114718848861e-06, ('SECS',
'M02G'): 1.3561294171657433e-06, ('EM02G', 'ALPG'): 0.0, ('ETHIG', 'SECS'):
-3.6013818449051294e-07, ('EAL0S', 'EAL0G'): 7.0076954873384941e-08, ('ALPS',
'M02S'): -1.521013508245142e-07, ('EALPS', 'ALPS'): -1.8176402553586655e-06,
('ALPG', 'EM02S'): 0.0, ('ALPG', 'M02S'): 1.4275390526229874e-07, ('KAPS',
'THIG'): -3.8296814440104735e-06, ('EAL0G', 'ETHIS'): 6.8483968659596611e-09,
('ALPG', 'ETHIS'): 4.3707938936467013e-07, ('M02G', 'ETHIG'):
2.0828003995404611e-05, ('ETHIS', 'EALPS'): 1.9021730648414777e-07, ('M02G',
'THIS'): 4.8170615531297017e-07, ('EM02S', 'SECG'): 0.0, ('M02S', 'ETHIS'):
-7.9804887592991765e-09, ('THIG', 'EAL0S'): 8.5274493541677246e-07, ('EM02S',
'KAPS'): 0.0, ('ALPS', 'KAPS'): 4.2289644797571366e-06, ('EALPS', 'EM02G'):
0.0, ('KAPS', 'EAL0G'): -3.510778024291086e-07, ('SECS', 'EALPS'):
-2.340529992704832e-07, ('ALPG', 'KAPS'): -1.3934242048266493e-05, ('THIS',
'THIS'): 1.3997371503287312e-08, ('EM02S', 'EALPS'): 0.0, ('EM02G', 'M02G'):
0.0, ('ETHIS', 'M02G'): -1.5645339041813079e-07, ('THIS', 'M02G'):
4.8170615531297017e-07, ('EALPS', 'SECS'): -2.340529992704832e-07, ('ALPG',
'ESECS'): 1.1416623328303995e-06, ('THIG', 'SECS'): -2.555178090391354e-07,
('EM02S', 'M02S'): 0.0, ('EAL0S', 'ETHIS'): 1.3718727952139344e-07, ('EM02G',
'EM02S'): 0.0, ('ESECS', 'THIS'): -1.6066088325636097e-08, ('THIG', 'M02S'):
3.1279039374527087e-08, ('SECG', 'THIG'): 2.625687795606487e-06, ('SECS',
'EAL0G'): 2.5736373354658603e-09, ('M02S', 'KAPS'): -4.2413128071918736e-07,
('ETHIG', 'SECG'): 1.8917717808263846e-05, ('SECS', 'EM02S'): 0.0, ('ESECS',
'M02S'): -2.7822034508540003e-08, ('ALPS', 'EAL0G'): 5.5343493256667721e-08,
('ETHIG', 'EM02G'): 0.0, ('THIG', 'M02G'): -8.6482378999837254e-06, ('ALPS',
'THIS'): 1.2044919766881962e-07, ('KAPS', 'THIS'): 1.9272024989328865e-07,
('SECG', 'ETHIS'): 3.2279231295898578e-07, ('ALPG', 'SECG'):
2.4115528065997216e-05, ('THIG', 'ETHIS'): 9.1667053921050091e-08, ('ESECS',
'SECS'): -4.3603850879636102e-08, ('EM02S', 'ESECS'): 0.0, ('ALPS', 'EM02S'):
0.0, ('EM02G', 'THIS'): 0.0, ('ESECS', 'ESECS'): 9.5474326522483904e-07,
('EM02G', 'KAPS'): 0.0, ('M02S', 'ETHIG'): -3.9261224522352571e-07, ('KAPS',
'SECS'): 4.8176197414287816e-07, ('THIS', 'SECS'): 6.08051419300745e-09,
('ETHIG', 'EAL0S'): 1.4808752765086284e-05, ('M02S', 'EM02G'): 0.0, ('M02G',
'EALPS'): -8.5952224552298213e-06, ('ALPG', 'M02G'): -3.0552028541221329e-05,
('KAPS', 'ALPG'): -1.3934242048266493e-05, ('EALPS', 'SECG'):
3.8961440993238121e-06, ('EALPS', 'EALPS'): 3.2550464734298223e-06, ('SECG',
'KAPS'): -1.3388869764770864e-05, ('EAL0G', 'SECG'): 3.0834373147791672e-08,
('THIG', 'EALPS'): 1.3088073710247217e-06, ('M02S', 'M02S'):
6.2362552502198413e-08, ('EAL0S', 'SECS'): -1.167954519937184e-07, ('ALPS',
'SECS'): 3.4950636427360906e-07, ('SECG', 'EM02S'): 0.0, ('THIG', 'THIS'):
-8.88154307904622e-08, ('EAL0S', 'ALPG'): 2.7705144957565145e-06, ('ESECS',
'EM02S'): 0.0, ('EAL0G', 'THIG'): -9.6714551617384486e-09, ('ETHIG', 'M02S'):
-3.9261224522352571e-07, ('ALPS', 'THIG'): -1.1100767579423132e-06, ('SECG',
'M02S'): 9.9430841783265254e-08, ('SECG', 'ETHIG'): 1.8917717808263846e-05,
('EAL0G', 'ESECS'): 1.6120423127841767e-08, ('ESECS', 'SECG'):
8.1249337700340916e-07, ('THIG', 'EM02G'): 0.0, ('EAL0G', 'M02S'):
-3.8621911523573612e-09, ('EALPS', 'ETHIG'): 6.7207370068910775e-06, ('ESECS',
'ALPS'): -2.3965071820922081e-08, ('SECS', 'ALPG'): -5.1164860952341495e-07,
('KAPS', 'ETHIG'): -4.6025488131970195e-05, ('ALPS', 'EALPS'):
-1.8176402553586655e-06, ('KAPS', 'EALPS'): -5.0707307468705406e-06, ('THIG',
'ALPS'): -1.1100767579423132e-06, ('EM02G', 'EAL0S'): 0.0, ('EALPS', 'ALPG'):
7.0052840588288252e-06, ('THIS', 'SECG'): -2.2622796924209372e-07, ('ETHIS',
'EM02G'): 0.0, ('THIS', 'EALPS'): -8.5440008282716879e-08, ('EAL0G', 'EAL0G'):
9.8766671656054909e-07, ('EM02G', 'M02S'): 0.0, ('EAL0G', 'KAPS'):
-3.510778024291086e-07, ('ETHIS', 'ETHIS'): 9.7101435963308939e-08, ('KAPS',
'ALPS'): 4.2289644797571366e-06, ('SECG', 'ALPG'): 2.4115528065997216e-05,
('M02G', 'THIG'): -8.6482378999837254e-06, ('THIS', 'ALPS'):
1.2044919766881962e-07, ('KAPS', 'KAPS'): 9.2778909027417343e-05, ('EM02S',
'EAL0G'): 0.0, ('EM02S', 'ETHIS'): 0.0, ('SECS', 'SECG'):
-6.3946744943539416e-07, ('KAPS', 'EAL0S'): -7.9008316901597777e-06, ('ETHIG',
'ALPG'): 4.061777554572628e-05, ('THIS', 'THIG'): -8.88154307904622e-08,
('EAL0S', 'M02S'): -9.7331124050211067e-09, ('THIS', 'ESECS'):
-1.6066088325636097e-08, ('SECG', 'SECG'): 2.0847143234661476e-05, ('THIS',
'EAL0S'): -4.5441474101194452e-08, ('SECG', 'EAL0S'): 2.9681814636767001e-06,
('ETHIG', 'M02G'): 2.0828003995404611e-05, ('ETHIG', 'ETHIS'):
1.9266583823151722e-06, ('EM02S', 'ALPG'): 0.0, ('ESECS', 'KAPS'):
-2.8598706302377214e-06, ('M02G', 'ALPG'): -3.0552028541221329e-05, ('ESECS',
'THIG'): 2.3150285229816843e-07, ('M02G', 'EM02S'): 0.0, ('SECS', 'ETHIG'):
-3.6013818449051294e-07, ('EM02G', 'EALPS'): 0.0, ('ALPG', 'EM02G'): 0.0,
('EAL0S', 'M02G'): -1.4607629894898336e-06, ('ETHIG', 'EALPS'):
6.7207370068910775e-06, ('ALPS', 'ETHIS'): -2.1440593643930602e-08, ('ETHIS',
'ALPS'): -2.1440593643930602e-08, ('THIS', 'KAPS'): 1.9272024989328865e-07,
('SECS', 'M02S'): -2.1587221004910876e-08, ('EAL0S', 'ESECS'):
3.1721058616921217e-07, ('ETHIG', 'EAL0G'): 8.651261809031964e-07, ('THIG',
'KAPS'): -3.8296814440104735e-06, ('M02S', 'ALPG'): 1.4275390526229874e-07,
('ETHIG', 'ETHIG'): 0.00022651126288389335, ('EAL0G', 'ALPS'):
5.5343493256667721e-08, ('M02G', 'EM02G'): 0.0, ('ESECS', 'EM02G'): 0.0,
('M02G', 'KAPS'): 1.1738921189923695e-05, ('EAL0S', 'EM02S'): 0.0, ('ETHIS',
'ETHIG'): 1.9266583823151722e-06, ('EAL0G', 'ALPG'): -1.7857710566950838e-07,
('EM02S', 'EAL0S'): 0.0, ('ALPS', 'EM02G'): 0.0, ('M02G', 'ALPS'):
1.4359224420333175e-05, ('SECS', 'EAL0S'): -1.167954519937184e-07, ('EAL0S',
'ALPS'): -4.800957442453362e-07, ('EM02S', 'EM02G'): 0.0, ('SECS', 'THIG'):
-2.555178090391354e-07, ('ETHIS', 'EAL0S'): 1.3718727952139344e-07, ('THIS',
'EM02S'): 0.0, ('EAL0G', 'EALPS'): 2.5535659974614765e-08, ('ESECS', 'EALPS'):
5.0979307037748514e-07, ('THIS', 'EAL0G'): 6.8966289500986728e-10, ('SECG',
'EAL0G'): 3.0834373147791672e-08, ('THIG', 'EAL0G'): -9.6714551617384486e-09,
('ETHIG', 'ESECS'): 5.155730523607306e-06, ('M02S', 'EAL0S'):
-9.7331124050211067e-09, ('EAL0G', 'SECS'): 2.5736373354658603e-09, ('ETHIS',
'EAL0G'): 6.8483968659596611e-09, ('M02G', 'ESECS'): -3.6717745362566989e-07,
('ETHIS', 'THIG'): 9.1667053921050091e-08, ('ESECS', 'M02G'):
-3.6717745362566989e-07, ('ALPS', 'SECG'): -3.014844665917385e-06, ('M02S',
'ESECS'): -2.7822034508540003e-08, ('EAL0S', 'EAL0S'): 2.5887129120424183e-06,
('ETHIS', 'EM02S'): 0.0, ('THIG', 'ALPG'): 1.0358173631970516e-05, ('EM02G',
'EAL0G'): 0.0, ('ALPG', 'SECS'): -5.1164860952341495e-07, ('ESECS', 'EAL0S'):
3.1721058616921217e-07, ('ETHIG', 'EM02S'): 0.0, ('ALPG', 'THIG'):
1.0358173631970516e-05}


# EIC12B  Fit to EICX + EICTSA + H1ZEUSindependentNEW  chi/dof=1833/1685

EIC12B = {
'EAL0G': 0.88119951947353059, 'KAPS': 1.1283716624991333, 'fix_ALPG': False,
'EDELM2S': 0.0, 'EPS': 2.0, 'EAL0S': 1.2317976488209421, 'KAPG': 0.0,
'fix_ETHIG': False, 'EPG': 2.0, 'EDELM2G': 0.0, 'SECG': -3.4817453300209258,
'EKAPG': 0.0, 'ESKEWG': 0.0, 'fix_KAPG': True, 'M02S': 0.19427256354081782,
'NG': 0.5, 'fix_KAPS': False, 'SECS': -0.083608237339463987, 'EKAPS': 0.0,
'NS': 0.152, 'fix_PG': True, 'ALPS': 0.11167606862318977, 'fix_EM02G': False,
'fix_DELM2S': True, 'SKEWG': 0.0, 'fix_EM02S': False, 'SKEWS': 0.0, 'fix_THIS':
False, 'fix_ESKEWS': True, 'fix_ETHIS': False, 'fix_NG': True, 'THIS':
-0.030511752920485059, 'fix_EDELM2S': True, 'fix_NS': True, 'fix_THIG': False,
'fix_DELB': True, 'THIG': 1.331662893895974, 'fix_EDELM2G': True, 'fix_EPG':
True, 'fix_ALPS': False, 'ETHIG': -7.4582832348708994, 'ESECS':
-0.10396938957836761, 'fix_EAL0G': False, 'fix_SECG': False, 'ETHIS':
-0.035205168536275096, 'ESECG': 0.0, 'fix_EPS': True, 'fix_EAL0S': False,
'fix_SECS': False, 'PS': 2.0, 'EALPG': 0.050000000000000003, 'EM02S':
0.10000133514361857, 'fix_AL0S': True, 'fix_EALPG': True, 'fix_M02G': False,
'EALPS': 0.053129945688844597, 'ESKEWS': 0.0, 'EM02G': 0.20000000000000004,
'PG': 2.0, 'fix_EALPS': False, 'limit_M02G': (0.10000000000000001, 1.5),
'fix_M02S': False, 'fix_DELM2G': True, 'fix_AL0G': True, 'AL0S': 1.1575,
'DELB': 0.0, 'fix_ESECS': False, 'fix_PS': True, 'limit_EM02S':
(0.10000000000000001, 1.5), 'fix_ESECG': True, 'AL0G': 1.2473000000000001,
'limit_EM02G': (0.10000000000000001, 1.5), 'limit_M02S': (0.10000000000000001,
1.5), 'M02G': 0.14526276614688835, 'fix_SKEWG': True, 'DELM2S': 0.0,
'fix_ESKEWG': True, 'ALPG': -0.048566809822429578, 'fix_SKEWS': True, 'DELM2G':
0.0}

EIC12Bcov={
('M02S', 'ALPS'): 1.2182227983578652e-05, ('M02G', 'M02S'):
-5.0350327284909212e-07, ('M02S', 'EM02S'): 0.0, ('M02G', 'SECS'):
3.8132092883847776e-07, ('ETHIS', 'ESECS'): 1.7825842003271457e-06, ('KAPS',
'M02G'): 9.6729993621235385e-06, ('ALPG', 'ALPG'): 4.4487921963803082e-07,
('ETHIS', 'THIS'): 9.9311328777963161e-08, ('EM02S', 'SECS'): 0.0, ('ESECS',
'ETHIG'): 0.00037587973996104607, ('M02S', 'EAL0G'): 1.8170567044309769e-05,
('ESECS', 'ALPG'): 3.0973023713498473e-06, ('EM02G', 'EM02G'):
0.26000000000000006, ('M02S', 'SECG'): 1.0967069480803222e-05, ('KAPS',
'ETHIS'): -1.5568462811365763e-05, ('KAPS', 'M02S'): -3.7149157019061347e-05,
('SECS', 'THIS'): -4.6011492274349477e-08, ('THIS', 'ETHIS'):
9.9311328777963161e-08, ('ETHIS', 'M02S'): 1.6661926788770194e-07, ('M02G',
'EAL0S'): 3.45975182993527e-06, ('ETHIS', 'KAPS'): -1.5568462811365763e-05,
('EM02G', 'ALPS'): 0.0, ('EALPS', 'ESECS'): -3.1367597130354759e-06, ('EALPS',
'EAL0S'): 1.5915091404271829e-05, ('EALPS', 'ETHIS'): 1.0749128855161958e-06,
('EM02S', 'ALPS'): 0.0, ('SECG', 'ESECS'): 0.00031730662131908082, ('ETHIG',
'KAPS'): -0.016514494132756796, ('EAL0G', 'ETHIG'): 0.014241695297110592,
('ETHIG', 'THIS'): 6.8326501589676203e-05, ('THIS', 'ALPG'):
-7.058869865413377e-07, ('EAL0S', 'THIG'): 0.00018770905536226335, ('SECS',
'SECS'): 2.4795248966309426e-07, ('ETHIG', 'ALPS'): 0.00078539928861290484,
('ALPG', 'EAL0S'): -1.9330591844328889e-05, ('KAPS', 'EM02S'): 0.0, ('EALPS',
'KAPS'): -6.2029631120822526e-05, ('THIS', 'M02S'): 2.1764385806685548e-08,
('ESECS', 'EAL0G'): 2.8651859854768175e-05, ('EAL0S', 'THIS'):
-2.2211991078337813e-07, ('ESECS', 'ETHIS'): 1.7825842003271457e-06, ('M02G',
'M02G'): 5.6250475892470296e-06, ('THIG', 'ESECS'): -0.00015391610558944792,
('SECS', 'EM02G'): 0.0, ('SECG', 'M02G'): -9.6963996376396266e-05, ('EAL0G',
'EAL0S'): -2.2729544996595242e-05, ('SECS', 'ETHIS'): -2.3392204899306951e-07,
('EAL0S', 'EALPS'): 1.5915091404271829e-05, ('M02S', 'THIG'):
-5.4392757689022933e-06, ('EM02G', 'SECS'): 0.0, ('EM02G', 'ETHIG'): 0.0,
('EALPS', 'THIG'): 6.5724172534693831e-05, ('THIS', 'ETHIG'):
6.8326501589676203e-05, ('EAL0G', 'THIS'): 4.1111129936089425e-06, ('M02G',
'EAL0G'): -2.3960938943971758e-05, ('THIG', 'ETHIG'): -0.048316434092084935,
('SECG', 'THIS'): 3.1216794912351147e-05, ('ETHIS', 'SECG'):
0.00014954056486185074, ('EALPS', 'EAL0G'): 2.2598616723376926e-05, ('ETHIS',
'SECS'): -2.3392204899306951e-07, ('EM02G', 'THIG'): 0.0, ('EALPS', 'THIS'):
-5.32640884256032e-08, ('SECG', 'EALPS'): -0.0001325132704094069, ('EM02S',
'THIG'): 0.0, ('ALPS', 'M02G'): 7.4696704980586091e-07, ('ALPS', 'EAL0S'):
9.83295349826178e-06, ('SECG', 'EM02G'): 0.0, ('EAL0G', 'EM02G'): 0.0,
('EALPS', 'M02S'): 1.1537087592398143e-05, ('EM02S', 'THIS'): 0.0, ('SECS',
'ALPS'): 4.3703648261053171e-07, ('SECG', 'ALPS'): 1.5256927042110562e-05,
('ALPS', 'ALPS'): 4.0338416796247979e-05, ('EALPS', 'EM02S'): 0.0, ('KAPS',
'SECG'): -0.0063084430516801536, ('ALPG', 'THIS'): -7.058869865413377e-07,
('THIG', 'EM02S'): 0.0, ('EM02G', 'ESECS'): 0.0, ('SECS', 'KAPS'):
9.6180167346968667e-06, ('EAL0S', 'ETHIG'): 6.3094808397105878e-05, ('M02S',
'M02G'): -5.0350327284909212e-07, ('THIG', 'THIG'): 0.011266642188316507,
('ALPS', 'ETHIG'): 0.00078539928861290484, ('EAL0S', 'EM02G'): 0.0, ('M02S',
'THIS'): 2.1764385806685548e-08, ('EM02S', 'M02G'): 0.0, ('ALPG', 'ETHIG'):
-0.0025913343112600072, ('M02G', 'ETHIS'): -4.0389684052454282e-07, ('EAL0S',
'KAPS'): -7.2923720678931347e-05, ('EAL0G', 'EM02S'): 0.0, ('EM02G', 'SECG'):
0.0, ('ALPG', 'EAL0G'): -0.00011035036339198878, ('M02S', 'EALPS'):
1.1537087592398143e-05, ('KAPS', 'ESECS'): -6.32536572311998e-07, ('M02G',
'SECG'): -9.6963996376396266e-05, ('EM02G', 'ETHIS'): 0.0, ('EM02S', 'EM02S'):
3.7383985667617615e-06, ('KAPS', 'EM02G'): 0.0, ('ALPG', 'ALPS'):
-3.5523799161627773e-06, ('ETHIS', 'ALPG'): -1.6508398207862056e-06, ('ALPS',
'ESECS'): -2.4002532399925858e-06, ('SECS', 'ESECS'): -5.265262785922216e-07,
('THIS', 'EM02G'): 0.0, ('M02S', 'SECS'): 2.9134705701296744e-08, ('ALPG',
'EALPS'): -9.0894095041296032e-07, ('EAL0S', 'SECG'): -0.00039313893495814445,
('EM02S', 'ETHIG'): 0.0, ('THIG', 'SECG'): -0.022893354753980522, ('SECG',
'SECS'): -7.5476211295304963e-05, ('EALPS', 'M02G'): -5.6516571270154893e-07,
('ETHIG', 'THIG'): -0.048316434092084935, ('EAL0G', 'M02G'):
-2.3960938943971758e-05, ('ALPS', 'ALPG'): -3.5523799161627773e-06, ('SECS',
'M02G'): 3.8132092883847776e-07, ('EM02G', 'ALPG'): 0.0, ('ETHIG', 'SECS'):
-0.00015251543760402432, ('EAL0S', 'EAL0G'): -2.2729544996595242e-05, ('ALPS',
'M02S'): 1.2182227983578652e-05, ('EALPS', 'ALPS'): 3.5428147638794067e-05,
('ALPG', 'EM02S'): 0.0, ('ALPG', 'M02S'): -3.4644084647012544e-06, ('KAPS',
'THIG'): 0.0031157299574384174, ('EAL0G', 'ETHIS'): 1.8067722631875045e-05,
('ALPG', 'ETHIS'): -1.6508398207862056e-06, ('M02G', 'ETHIG'):
-0.00016093449954286234, ('ETHIS', 'EALPS'): 1.0749128855161958e-06, ('M02G',
'THIS'): 1.8970862233838484e-08, ('EM02S', 'SECG'): 0.0, ('M02S', 'ETHIS'):
1.6661926788770194e-07, ('THIG', 'EAL0S'): 0.00018770905536226335, ('EM02S',
'KAPS'): 0.0, ('ALPS', 'KAPS'): -8.3071644183003372e-05, ('EALPS', 'EM02G'):
0.0, ('KAPS', 'EAL0G'): -0.00099299307141616092, ('SECS', 'EALPS'):
3.4376911722699645e-07, ('ALPG', 'KAPS'): 0.0002050351833982663, ('THIS',
'THIS'): 3.3380249136080808e-08, ('EM02S', 'EALPS'): 0.0, ('EM02G', 'M02G'):
0.0, ('ETHIS', 'M02G'): -4.0389684052454282e-07, ('THIS', 'M02G'):
1.8970862233838484e-08, ('EALPS', 'SECS'): 3.4376911722699645e-07, ('ALPG',
'ESECS'): 3.0973023713498473e-06, ('THIG', 'SECS'): 3.6530110511975491e-05,
('EM02S', 'M02S'): 0.0, ('EAL0S', 'ETHIS'): -4.2434413041232948e-06, ('EM02G',
'EM02S'): 0.0, ('ESECS', 'THIS'): 1.9493406358142285e-07, ('THIG', 'M02S'):
-5.4392757689022933e-06, ('SECG', 'THIG'): -0.022893354753980522, ('SECS',
'EAL0G'): -9.7267079087619458e-06, ('M02S', 'KAPS'): -3.7149157019061347e-05,
('ETHIG', 'SECG'): 0.098039614318413029, ('SECS', 'EM02S'): 0.0, ('ESECS',
'M02S'): -2.2852799968356889e-06, ('ALPS', 'EAL0G'): 3.5921227004583957e-05,
('ETHIG', 'EM02G'): 0.0, ('THIG', 'M02G'): 4.2772976334286709e-05, ('ALPS',
'THIS'): 1.4562441831092659e-07, ('KAPS', 'THIS'): -4.5036063403289767e-06,
('SECG', 'ETHIS'): 0.00014954056486185074, ('ALPG', 'SECG'):
-0.0010653580528777471, ('THIG', 'ETHIS'): -7.2950016671734195e-05, ('ESECS',
'SECS'): -5.265262785922216e-07, ('EM02S', 'ESECS'): 0.0, ('ALPS', 'EM02S'):
0.0, ('EM02G', 'THIS'): 0.0, ('ESECS', 'ESECS'): 1.0602313345364263e-05,
('EM02G', 'KAPS'): 0.0, ('M02S', 'ETHIG'): 0.00032373508105400727, ('KAPS',
'SECS'): 9.6180167346968667e-06, ('THIS', 'SECS'): -4.6011492274349477e-08,
('ETHIG', 'EAL0S'): 6.3094808397105878e-05, ('M02S', 'EM02G'): 0.0, ('M02G',
'EALPS'): -5.6516571270154893e-07, ('ALPG', 'M02G'): -2.4867138475405008e-06,
('KAPS', 'ALPG'): 0.0002050351833982663, ('EALPS', 'SECG'):
-0.0001325132704094069, ('EALPS', 'EALPS'): 4.1647426652465632e-05, ('SECG',
'KAPS'): -0.0063084430516801536, ('EAL0G', 'SECG'): 0.0060503152309174645,
('THIG', 'EALPS'): 6.5724172534693831e-05, ('M02S', 'M02S'):
4.4254730440042117e-06, ('EAL0S', 'SECS'): 6.3231859346972336e-07, ('ALPS',
'SECS'): 4.3703648261053171e-07, ('SECG', 'EM02S'): 0.0, ('THIG', 'THIS'):
-1.5539333308444255e-05, ('EAL0S', 'ALPG'): -1.9330591844328889e-05, ('ESECS',
'EM02S'): 0.0, ('EAL0G', 'THIG'): -0.0029853383567350214, ('ETHIG', 'M02S'):
0.00032373508105400727, ('ALPS', 'THIG'): -7.6380523898412212e-06, ('SECG',
'M02S'): 1.0967069480803222e-05, ('SECG', 'ETHIG'): 0.098039614318413029,
('EAL0G', 'ESECS'): 2.8651859854768175e-05, ('ESECS', 'SECG'):
0.00031730662131908082, ('THIG', 'EM02G'): 0.0, ('EAL0G', 'M02S'):
1.8170567044309769e-05, ('EALPS', 'ETHIG'): 0.00047932672925612415, ('ESECS',
'ALPS'): -2.4002532399925858e-06, ('SECS', 'ALPG'): 1.8083989627751003e-06,
('KAPS', 'ETHIG'): -0.016514494132756796, ('ALPS', 'EALPS'):
3.5428147638794067e-05, ('KAPS', 'EALPS'): -6.2029631120822526e-05, ('THIG',
'ALPS'): -7.6380523898412212e-06, ('EM02G', 'EAL0S'): 0.0, ('EALPS', 'ALPG'):
-9.0894095041296032e-07, ('THIS', 'SECG'): 3.1216794912351147e-05, ('ETHIS',
'EM02G'): 0.0, ('THIS', 'EALPS'): -5.32640884256032e-08, ('EAL0G', 'EAL0G'):
0.0010551456668108385, ('EM02G', 'M02S'): 0.0, ('EAL0G', 'KAPS'):
-0.00099299307141616092, ('ETHIS', 'ETHIS'): 1.0101445408919546e-06, ('KAPS',
'ALPS'): -8.3071644183003372e-05, ('SECG', 'ALPG'): -0.0010653580528777471,
('M02G', 'THIG'): 4.2772976334286709e-05, ('THIS', 'ALPS'):
1.4562441831092659e-07, ('KAPS', 'KAPS'): 0.001296883387521679, ('EM02S',
'EAL0G'): 0.0, ('EM02S', 'ETHIS'): 0.0, ('SECS', 'SECG'):
-7.5476211295304963e-05, ('KAPS', 'EAL0S'): -7.2923720678931347e-05, ('ETHIG',
'ALPG'): -0.0025913343112600072, ('THIS', 'THIG'): -1.5539333308444255e-05,
('EAL0S', 'M02S'): 6.8866252118406661e-06, ('THIS', 'ESECS'):
1.9493406358142285e-07, ('SECG', 'SECG'): 0.04659658941237356, ('THIS',
'EAL0S'): -2.2211991078337813e-07, ('SECG', 'EAL0S'): -0.00039313893495814445,
('ETHIG', 'M02G'): -0.00016093449954286234, ('ETHIG', 'ETHIS'):
0.00029108361707214699, ('EM02S', 'ALPG'): 0.0, ('ESECS', 'KAPS'):
-6.32536572311998e-07, ('M02G', 'ALPG'): -2.4867138475405008e-06, ('ESECS',
'THIG'): -0.00015391610558944792, ('M02G', 'EM02S'): 0.0, ('SECS', 'ETHIG'):
-0.00015251543760402432, ('EM02G', 'EALPS'): 0.0, ('ALPG', 'EM02G'): 0.0,
('EAL0S', 'M02G'): 3.45975182993527e-06, ('ETHIG', 'EALPS'):
0.00047932672925612415, ('ALPS', 'ETHIS'): 1.5075474397646113e-06, ('ETHIS',
'ALPS'): 1.5075474397646113e-06, ('THIS', 'KAPS'): -4.5036063403289767e-06,
('SECS', 'M02S'): 2.9134705701296744e-08, ('EAL0S', 'ESECS'):
-1.968799721751056e-05, ('ETHIG', 'EAL0G'): 0.014241695297110592, ('THIG',
'KAPS'): 0.0031157299574384174, ('M02S', 'ALPG'): -3.4644084647012544e-06,
('ETHIG', 'ETHIG'): 0.23342450520431696, ('EAL0G', 'ALPS'):
3.5921227004583957e-05, ('M02G', 'EM02G'): 0.0, ('ESECS', 'EM02G'): 0.0,
('M02G', 'KAPS'): 9.6729993621235385e-06, ('EAL0S', 'EM02S'): 0.0, ('ETHIS',
'ETHIG'): 0.00029108361707214699, ('EAL0G', 'ALPG'): -0.00011035036339198878,
('EM02S', 'EAL0S'): 0.0, ('ALPS', 'EM02G'): 0.0, ('M02G', 'ALPS'):
7.4696704980586091e-07, ('SECS', 'EAL0S'): 6.3231859346972336e-07, ('EAL0S',
'ALPS'): 9.83295349826178e-06, ('EM02S', 'EM02G'): 0.0, ('SECS', 'THIG'):
3.6530110511975491e-05, ('ETHIS', 'EAL0S'): -4.2434413041232948e-06, ('THIS',
'EM02S'): 0.0, ('EAL0G', 'EALPS'): 2.2598616723376926e-05, ('ESECS', 'EALPS'):
-3.1367597130354759e-06, ('THIS', 'EAL0G'): 4.1111129936089425e-06, ('SECG',
'EAL0G'): 0.0060503152309174645, ('THIG', 'EAL0G'): -0.0029853383567350214,
('ETHIG', 'ESECS'): 0.00037587973996104607, ('M02S', 'EAL0S'):
6.8866252118406661e-06, ('EAL0G', 'SECS'): -9.7267079087619458e-06, ('ETHIS',
'EAL0G'): 1.8067722631875045e-05, ('M02G', 'ESECS'): -1.0454176191646819e-06,
('ETHIS', 'THIG'): -7.2950016671734195e-05, ('ESECS', 'M02G'):
-1.0454176191646819e-06, ('ALPS', 'SECG'): 1.5256927042110562e-05, ('M02S',
'ESECS'): -2.2852799968356889e-06, ('EAL0S', 'EAL0S'): 5.5022304354873575e-05,
('ETHIS', 'EM02S'): 0.0, ('THIG', 'ALPG'): 0.0005222737131516835, ('EM02G',
'EAL0G'): 0.0, ('ALPG', 'SECS'): 1.8083989627751003e-06, ('ESECS', 'EAL0S'):
-1.968799721751056e-05, ('ETHIG', 'EM02S'): 0.0, ('ALPG', 'THIG'):
0.0005222737131516835}
