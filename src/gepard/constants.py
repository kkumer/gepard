"""Constants, with their numerical values, not to be changed by user.

Also, fixed dictionaries.
"""

from scipy.special import zeta  # type: ignore

# QCD parameters
NC = 3
CF = (NC**2 - 1) / (2 * NC)
CA = NC
CG = CF - CA/2
TF = 0.5

# Meson decay constants
F_rho0 = 0.209

# Parameters needed for harmonic sums etc.
ZETA2 = zeta(2)
ZETA3 = zeta(3)
ZETA4 = zeta(4)

Mp = 0.938272013  # proton mass
Mp2 = Mp ** 2
alpha = 1.0 / 137.036  # inverse of fine structure constant
# (hbar^2 c^2 / GeV^2)  in nanobarns
GeV2nb = 389379.0
GeVfm = 0.197327

tolerance2 = 1.0  # sqrt(2*Npoints)/1.65^2 according to GJV/MRST


# Dictionary that maps variables to their LaTeX representation
toTeX = {
    "tm": "$-t$",
    "xB": "$x_{\\rm B}$",
    "Q2": "$Q^{2}$",
    "phi": "$\\phi$",
    "t": "$t$",
    "xi": "$\\xi$",
    "W": "$W$",
    "xixB": "$\\xi=x_{B}/(2-x_{B})$",
    "TotalCrossSection": "$\\sigma$",
    "PartialCrossSection": "$d\\sigma / dt$",
    "ALUIsin1": "$A_{LU,I}^{\\sin 1\\phi}$",
    "ALUsin1": "$A_{LU}^{\\sin\\phi}$",
    "ALU": "$A_{\\rm LU}$",
    "ALUI": "$A_{\\rm LU}^{\\rm I}$",
    "ALUDVCS": "$A_{\\rm LU}^{\\rm DVCS}$",
    "XLU": "$\\Delta\\Sigma/d\\Phi$",
    "XUU": "$d\\sigma/d\\Phi$",
    "ImCI": "$Im\\mathcal{C}^{I}$",
    "ReCI": "$Re\\mathcal{C}^{I}$",
    "XLUw2C": "$Im\\mathcal{C}^{I}$",
    "XUUw2C": "$Re\\mathcal{C}^{I}$, $Re(\\mathcal{C}+\\Delta\\mathcal{C})^{I}$",
    "XUUw": "$XUUw$",
    "XLUw": "$XLUw$",
    "ReCpDCI": "$Re(\\mathcal{C}^{I}+\\Delta\\mathcal{C}^{I})$",
    "XwA": " $b_{1}/b_{0}$ (XwA)",
    "AC": "$AC$",
    "BCSA": "$BCSA$",
    "TSA": "$TSA$",
    "BTSA": "$BTSA$",
    "ACcos0": "$A_{C}^{\\cos 0\\phi}$",
    "ACcos1": "$A_{C}^{\\cos\\phi}$",
    "XDVCSt": r"$d\sigma_{\rm DVCS} / dt$",
    "XDVCS": r"$\sigma_{\rm DVCS}$",
    "X": r"$\sigma_{\rm DVCS}$",
    "Xt": r"$d\sigma_{\rm DVCS}/dt$",
    "npt": "point no.",
    "ImH": "$\\mathfrak{Im}\\,\\mathcal{H}$",
    "ReH": "$\\mathfrak{Re}\\,\\mathcal{H}$",
    "ImE": "$\\mathfrak{Im}\\,\\mathcal{E}$",
    "ReE": "$\\mathfrak{Re}\\,\\mathcal{E}$",
    "ImHt": "$\\mathfrak{Im}\\,\\widetilde{\\mathcal{H}}$",
    "ReHt": "$\\mathfrak{Re}\\,\\widetilde{\\mathcal{H}}$",
    "ImEt": "$\\mathfrak{Im}\\,\\widetilde{\\mathcal{E}}$",
    "ReEt": "$\\mathfrak{Re}\\,\\widetilde{\\mathcal{E}}$",
    "ImHu": "$\\mathfrak{Im}\\,\\mathcal{H}_{u}$",
    "ReHu": "$\\mathfrak{Re}\\,\\mathcal{H}_{u}$",
    "ImEu": "$\\mathfrak{Im}\\,\\mathcal{E}_{u}$",
    "ReEu": "$\\mathfrak{Re}\\,\\mathcal{E}_{u}$",
    "ImHtu": "$\\mathfrak{Im}\\,\\widetilde{\\mathcal{H}}_{u}$",
    "ReHtu": "$\\mathfrak{Re}\\,\\widetilde{\\mathcal{H}}_{u}$",
    "ImEtu": "$\\mathfrak{Im}\\,\\widetilde{\\mathcal{E}}_{u}$",
    "ReEtu": "$\\mathfrak{Re}\\,\\widetilde{\\mathcal{E}}_{u}$",
    "ImHd": "$\\mathfrak{Im}\\,\\mathcal{H}_{d}$",
    "ReHd": "$\\mathfrak{Re}\\,\\mathcal{H}_{d}$",
    "ImEd": "$\\mathfrak{Im}\\,\\mathcal{E}_{d}$",
    "ReEd": "$\\mathfrak{Re}\\,\\mathcal{E}_{d}$",
    "ImHtd": "$\\mathfrak{Im}\\,\\widetilde{\\mathcal{H}}_{d}$",
    "ReHtd": "$\\mathfrak{Re}\\,\\widetilde{\\mathcal{H}}_{d}$",
    "ImEtd": "$\\mathfrak{Im}\\,\\widetilde{\\mathcal{E}}_{d}$",
    "ReEtd": "$\\mathfrak{Re}\\,\\widetilde{\\mathcal{E}}_{d}$",
}

# Dictionary that maps observable tuples to their LaTeX representation
OBStoTeX = {
    ("AC", 0): "$A_{\\rm C}^{\\cos 0\\phi}$",
    ("AC", 1): "$A_{\\rm C}^{\\cos \\phi}$",
    ("ALU", -1): "$A_{\\rm LU}^{\\sin \\phi}$",
    ("ALUI", -1): "$A_{\\rm LU,I}^{\\sin \\phi}$",
    ("ALUI", -2): "$A_{\\rm LU,I}^{\\sin 2\\phi}$",
    ("ALUDVCS", -1): "$A_{\\rm LU,DVCS}^{\\sin \\phi}$",
    ("TSA", -1): "$A_{\\rm UL}^{\\sin \\phi}$",
    ("BTSA", 0): "$A_{\\rm LL}^{\\cos 0 \\phi}$",
    ("BTSA", 1): "$A_{\\rm LL}^{\\cos  \\phi}$",
    ("AUTI", 1): "$A_{\\rm UT,I}^{\\sin(\\phi-\\phi_S) \\cos \\phi}$",
    ("AUTI", 0): "$A_{\\rm UT,I}^{\\sin(\\phi-\\phi_S)}$",
    ("AUTI", -1): "$A_{\\rm UT,I}^{\\cos(\\phi-\\phi_S) \\sin \\phi}$",
    ("AUTDVCS", 0): "$A_{\\rm UT,DVCS}^{\\sin(\\phi-\\phi_S)}$",
    ("XLUw", -1): r"$\Delta\sigma^{\sin\phi,w}$",
    ("XUUw", 0): r"$d\sigma^{\cos 0\phi,w}$",
    ("XUUw", 1): r"$d\sigma^{\cos\phi,w}$",
}
