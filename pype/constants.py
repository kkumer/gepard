"""
Constants, with their numerical values, not to be changed by user.
Maybe one should make sure they stay constant.

Also fixed dictionaries 
"""

Mp = 0.938272013   # proton mass
Mp2 = Mp**2
alpha = 1./137.036  # inverse of fine structure constant
# WolframAlpha: (hbar^2 c^2 / GeV^2)  in nanobarns
GeV2nb = 389379.
GeVfm = 0.197327

tolerance2 = 1.0   # sqrt(2*Npoints)/1.65^2 according to GJV/MRST


# Dictionary that maps variables to their LaTeX representation
toTeX = {'tm' : '$-t$',  'xB' : '$x_{\\rm B}$', 'Q2' : '$Q^{2}$', 
     'phi' : '$\\phi$', 't' : '$t$',
     'xi': '$\\xi$',  'W': '$W$',
     'xixB': '$\\xi=x_{B}/(2-x_{B})$',
     "TotalCrossSection" : '$\\sigma$',
     "PartialCrossSection" : '$d\\sigma / dt$',
     "ALUIsin1" : '$A_{LU,I}^{\\sin 1\\phi}$',
     "BSAsin1" : '$A_{LU}^{\\sin\\phi}$',
     "ALU" : '$A_{\\rm LU}$',
     "ALUI" : '$A_{\\rm LU}^{\\rm I}$',
     "ALUDVCS" : '$A_{\\rm LU}^{\\rm DVCS}$',
     "XLU" : '$d\\Sigma/d\\Phi$',
     "XUU" : '$d\\sigma/d\\Phi$',
     "ImCI" : '$Im\\mathcal{C}^{I}$',
     "ReCI" : '$Re\\mathcal{C}^{I}$',
     "BSDw2C" : '$Im\\mathcal{C}^{I}$',
     "BSSw2C" : '$Re\\mathcal{C}^{I}$, $Re(\\mathcal{C}+\\Delta\\mathcal{C})^{I}$',
     "BSSw" : '$BSSw$',
     "BSDw" : '$BSDw$',
     "ReCpDCI" : '$Re(\\mathcal{C}^{I}+\\Delta\\mathcal{C}^{I})$',
     "XwA" : ' $b_{1}/b_{0}$ (XwA)',
     "BCA" : '$BCA$',
     "BCSA" : '$BCSA$',
     "BSA" : '$BSA$',
     "BSD" : '$BSD$',
     "BSS" : '$BSS$',
     "TSA" : '$TSA$',
     "BTSA" : '$BTSA$',
     "BCAcos0" : '$A_{C}^{\\cos 0\\phi}$',
     "BCAcos1" : '$A_{C}^{\\cos\\phi}$',
     "XDVCSt" : '$d\\sigma_{DVCS} / dt$',
     "XDVCS" : '$\\sigma_{DVCS}$',
     "X" : '$\\sigma_{DVCS}$',
     "npt" : 'point no.',
     "ImH" : '$\\mathfrak{Im}\\mathcal{H}(\\xi, t)$',
     "ReH" : '$\\mathfrak{Re}\\mathcal{H}(\\xi, t)$',
     "ImE" : '$\\mathfrak{Im}\\mathcal{E}(\\xi, t)$',
     "ReE" : '$\\mathfrak{Re}\\mathcal{E}(\\xi, t)$',
     "ImHt" : '$\\mathfrak{Im}\\tilde{\\mathcal{H}}(\\xi, t)$',
     "ReHt" : '$\\mathfrak{Re}\\tilde{\\mathcal{H}}(\\xi, t)$',
     "ImEt" : '$\\mathfrak{Im}\\tilde{\\mathcal{E}}(\\xi, t)$',
     "ReEt" : '$\\mathfrak{Re}\\tilde{\\mathcal{E}}(\\xi, t)$',
     }

# Dictionary that maps observable tuples to their LaTeX representation
OBStoTeX = {
        ('BCA', 0) : '$A_{\\rm C}^{\\cos 0\\phi}$',
        ('BCA', 1) : '$A_{\\rm C}^{\\cos \\phi}$',
        ('BSA', -1) : '$A_{\\rm LU}^{\\sin \\phi}$',
        ('ALUI', -1) : '$A_{\\rm LU,I}^{\\sin \\phi}$',
        ('ALUI', -2) : '$A_{\\rm LU,I}^{\\sin 2\\phi}$',
        ('ALUDVCS', -1) : '$A_{\\rm LU,DVCS}^{\\sin \\phi}$',
        ('TSA', -1) : '$A_{\\rm UL}^{\\sin \\phi}$',
        ('BTSA', 0) : '$A_{\\rm LL}^{\\cos 0 \\phi}$',
        ('BTSA', 1) : '$A_{\\rm LL}^{\\cos  \\phi}$',
        ('AUTI', 1) : '$A_{\\rm UT,I}^{\\sin(\\phi-\\phi_S) \\cos \\phi}$',
        ('AUTI', 0) : '$A_{\\rm UT,I}^{\\sin(\\phi-\\phi_S)}$',
        ('AUTI', -1) : '$A_{\\rm UT,I}^{\\cos(\\phi-\\phi_S) \\sin \\phi}$',
        ('AUTDVCS', 0) : '$A_{\\rm UT,DVCS}^{\\sin(\\phi-\\phi_S)}$'
        }