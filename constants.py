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
     "npt" : 'point no.',
     "ImH" : '$\\mathfrak{Im}\\mathcal{H}(\\xi, t)$',
     "ReH" : '$\\mathfrak{Re}\\mathcal{H}(\\xi, t)$',
     "ImHt" : '$\\mathfrak{Im}\\tilde{\\mathcal{H}}(\\xi, t)$',
     "ReHt" : '$\\mathfrak{Re}\\tilde{\\mathcal{H}}(\\xi, t)$',
     }
