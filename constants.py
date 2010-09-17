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



# Dictionary that maps variables to their LaTeX representation
toTeX = {'tm' : '$-t$',  'xB' : '$x_{\\rm B}$', 'Q2' : '$Q^{2}$', 
     'phi' : '$\\phi$', 't' : '$t$',
     'xi': '$\\xi$',  'W': '$W$',
     "TotalCrossSection" : '$\\sigma$',
     "PartialCrossSection" : '$d\\sigma / dt$',
     "ALUIsin1" : '$A_{LU,I}^{\\sin 1\\phi}$',
     "ALU" : '$A_{LU}$',
     "XLU" : '$d\\Sigma/d\\Phi$',
     "XUU" : '$d\\sigma/d\\Phi$',
     "ImCI" : '$Im\\mathcal{C}^{I}$',
     "ReCI" : '$Re\\mathcal{C}^{I}$',
     "BSDw2C" : '$Im\\mathcal{C}^{I}$',
     "BSSw2C" : '$Re\\mathcal{C}^{I}$, $Re(\\mathcal{C}+\\Delta\\mathcal{C})^{I}$',
     "BSSw" : '$BSSw$',
     "ReCpDCI" : '$Re(\\mathcal{C}^{I}+\\Delta\\mathcal{C}^{I})$',
     "XwA" : ' $b_{1}/b_{0}$ (XwA)',
     "BCA" : '$BCA$',
     "BSA" : '$BSA$',
     "BSD" : '$BSD$',
     "BSS" : '$BSS$',
     "TSA" : '$TSA$',
     "BTSA" : '$BTSA$',
     "BCAcos0" : '$A_{C}^{\\cos 0\\phi}$',
     "BCAcos1" : '$A_{C}^{\\cos 1\\phi}$',
     "XDVCSt" : '$d\\sigma_{DVCS} / dt$',
     "XDVCS" : '$\\sigma_{DVCS}$'
     }
