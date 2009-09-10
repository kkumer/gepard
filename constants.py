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
toTeX = {'mt' : '$-t$',  'xB' : '$x_{\\rm B}$', 'Q2' : '$Q^{2}$', 
     'phi' : '$\\phi$', 't' : '$t$',
     'xi': '$\\xi$', 
     "TotalCrossSection" : '$\\sigma$',
     "PartialCrossSection" : '$d\\sigma / dt$',
     "ALUIsin1" : '$A_{LU,I}^{\\sin 1\\phi}$',
     "ALUa" : '$A_{LU}(90^{\\degree})$',
     "ALU" : '$A_{LU}$',
     "XLU" : '$d\\Sigma/d\\Phi$',
     "XUU" : '$d\\sigma/d\\Phi$',
     "ImCI" : '$Im\\mathcal{C}^{I}$',
     "ReCI" : '$Re\\mathcal{C}^{I}$',
     "ReCpDCI" : '$Re(\\mathcal{C}^{I}+\\Delta\\mathcal{C}^{I})$',
     "b1ovb0" : '$b_{1}/b_{0}$',
     "BCAcos0" : '$A_{C}^{\\cos 0\\phi}$',
     "BCAcos1" : '$A_{C}^{\\cos 1\\phi}$'
     }
