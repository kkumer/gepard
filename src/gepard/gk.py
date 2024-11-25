"""Goloskokov-Kroll GPD/CFF model.

From [arXiv:1210.6975], Kroll:2012sm
Real part of CFFs is obtained using dispersion relations.
"""

from numpy import (array, exp, log, ndarray, pi, sqrt)
from numpy import sum as npsum
from scipy.special import beta, gamma, j0, j1  # type: ignore
from . import cff
from .constants import GeVfm, Mp
from .quadrature import bquadrature


class GoloskokovKrollCFF(cff.DispersionCFF):
    """Goloskokov-Kroll GPD/CFF model.

    Mostly from [arXiv:1210.6975], Kroll:2012sm
    Real part of CFFs is obtained using dispersion relations.

    """
    def __init__(self, **kwargs) -> None:
        self.nf = 3  # not used! Needed not to trigger error
        super().__init__(**kwargs)

    def _intval(self, x, eta, alt, j, zero=False):
        """Analytic integral over DD for valence sector."""
        x1 = (x+eta)/(1+eta)
        p = - alt + j/2
        if zero:
            x2 = 0
        else:
            x2 = (x-eta)/(1-eta)
        prefac = (3./2)*gamma(1+p)/gamma(4+p)
        if not zero and (eta/x < 1e-4):
            # use Taylor up to O(eta)
            b = (-2./3)*(p**3+6*p**2+11*p+6)*(x-1)**3*x**p
        else:
            # exact (but unstable for small eta)
            b = ((eta**2-x)*(x1**(2+p)-x2**(2+p))
                 + (2+p)*eta*(1-x)*(x1**(2+p)+x2**(2+p)))/eta**3
        return prefac * b

    def _intsea(self, x, eta, alt, j, zero=False):
        """Analytic integral over DD for sea sector."""
        x1 = (x+eta)/(1+eta)
        p = - alt + j/2
        if zero:
            x2 = 0
        else:
            x2 = (x-eta)/(1-eta)
        if not zero and (eta/x < 1e-2):
            # use Taylor up to O(eta^2)
            b = eta**2*(p-1)*p+x*(eta**2*(p+6)*(p*(x-2)+7*x)+14*x)
            return (1-x)**5 * b * x**(p-2) / 14.
        else:
            # exact (but unstable for small eta)
            prefac = (15./2)*gamma(1+p)/gamma(6+p)/eta**5
            brm = 3*(eta**2-x)**2+eta**2*(p**2+6*p+8)*(1-x)**2
            brp = 3*eta*(eta**2-x)*(3+p)*(1-x)
            b = brm*(x1**(3+p)-x2**(3+p))+brp*(x1**(3+p)+x2**(3+p))
            return prefac * b

    def _val(self, x, eta, alt, j):
        # assert eta >= 0
        if x >= eta:
            return self._intval(x, eta, alt, j)
        elif -eta < x < eta:
            return self._intval(x, eta, alt, j, zero=True)
        else:
            return 0

    def _sea(self, x, eta, alt, j, glue=-1):
        DMKILL = 1  # Set to 0 to agree with DM
        assert eta >= 0
        if x >= eta:
            return self._intsea(x, eta, alt, j)
        elif -eta < x < eta:
            return (self._intsea(x, eta, alt, j, zero=True)
                    + glue * DMKILL * self._intsea(-x, eta, alt, j, zero=True))
        else:
            return glue * DMKILL * self._intsea(-x, eta, alt, j)  # x < -eta

    #   ######  GPDs  #########

    #   ===  GPD H  ===

    def Huval(self, x, eta, t, Q2):
        """GK12 model for H^u_val GPD."""
        Q02 = 4.
        L = log(Q2/Q02)
        cuval = array([1.52+0.248*L, 2.88-0.940*L, -0.095*L, 0])
        al0, alp = (0.48, 0.9)
        alt = al0 + alp*t
        xs = array([self._val(x, eta, alt, j) for j in range(4)])
        return npsum(cuval*xs)

    def Hdval(self, x, eta, t, Q2):
        """GK12 model for H^d_val GPD."""
        Q02 = 4.
        L = log(Q2/Q02)
        cdval = array([0.76+0.248*L, 3.11-1.36*L, -3.99+1.15*L, 0])
        al0, alp = (0.48, 0.9)
        alt = al0 + alp*t
        xs = array([self._val(x, eta, alt, j) for j in range(4)])
        return npsum(cdval*xs)

    def Hs(self, x, eta, t, Q2):
        """GK12 model for strange sea GPD H."""
        Q02 = 4.
        Mp2 = 0.9383**2
        L = log(Q2/Q02)
        cs = array([0.123+0.0003*L, -0.327-0.004*L, 0.692-0.068*L, -0.486+0.038*L])
        al0, alp = (1.10+0.06*L-0.0027*L**2, 0.15)
        alt = al0 + alp*t
        xs = array([self._sea(x, eta, alt, j) for j in range(4)])
        b = 2.58 + 0.25*log(Mp2/(Q2+Mp2))
        return exp(b*t) * npsum(cs*xs)

    def Hudsea(self, x, eta, t, Q2):
        """GK12 model for u or d sea contribution to GPD H."""
        Q02 = 4.
        kaps = 1.+0.68/(1.+0.52*log(Q2/Q02))
        return kaps * self.Hs(x, eta, t, Q2)

    def Hu(self, x, eta, t, Q2):
        """GK12 model for u-quark GPD H."""
        return self.Huval(x, eta, t, Q2) + self.Hudsea(x, eta, t, Q2)

    def Hd(self, x, eta, t, Q2):
        """GK12 model for d-quark GPD H."""
        return self.Hdval(x, eta, t, Q2) + self.Hudsea(x, eta, t, Q2)

    def Hg(self, x, eta, t, Q2):
        """GK12 model for gluon GPD H."""
        Q02 = 4.
        Mp2 = 0.9383**2
        L = log(Q2/Q02)
        cs = array([2.23+0.362*L, 5.43-7.0*L, -34+22.5*L, 40.6-21.6*L])
        al0, alp = (1.10+0.06*L-0.0027*L**2, 0.15)
        alt = al0 + alp*t - 1
        xs = array([self._sea(x, eta, alt, j, glue=+1) for j in range(4)])
        b = 2.58 + 0.25*log(Mp2/(Q2+Mp2))
        return exp(b*t) * sum(cs*xs)

    #   ===  GPD H tilde  ===

    def _invAq(self, al, cs):
        coefs = [1, (1-al)/(5-al), (2-al)*(1-al)/(6-al)/(5-al)]
        return beta(1-al, 4)*npsum(cs*coefs)

    def Htuval(self, x, eta, t, Q2):
        """GK12 model for Ht^u_val GPD."""
        Q02 = 4.
        L = log(Q2/Q02)
        ctuval = array([0.17+0.03*L, 1.34-0.02*L, 0.12-0.4*L])
        al0, alp = (0.48, 0.45)
        alt = al0 + alp*t
        # val() expects series in beta^(j/2), so 0,2,4 -> 0,1,2
        xs = array([self._val(x, eta, alt, j) for j in [0, 2, 4]])
        return 0.926*npsum(ctuval*xs)/self._invAq(al0, ctuval)

    def Htdval(self, x, eta, t, Q2):
        """GK12 model for Ht^u_val GPD."""
        Q02 = 4.
        L = log(Q2/Q02)
        ctuval = array([-0.32-0.04*L, -1.427-0.176*L, 0.692-0.068*L])
        al0, alp = (0.48, 0.45)
        alt = al0 + alp*t
        # val() expects series in beta^(j/2), so 0,2,4 -> 0,1,2
        xs = array([self._val(x, eta, alt, j) for j in [0, 2, 4]])
        return -0.341*npsum(ctuval*xs)/self._invAq(al0, ctuval)

    Htu = Htuval
    Htd = Htdval

    def Hts(self, x, eta, t, Q2):
        """GK12 model for Ht^s GPD."""
        return 0.

    def Htg(self, x, eta, t, Q2):
        """GK12 model for gluon GPD Htilde."""
        # Coeffs are from hep-ph/0501242
        Q02 = 4.
        Mp2 = 0.9383**2
        L = log(Q2/Q02)
        cs = array([3.39-0.864*L, 1.73+0.24*L-0.17*L**2, 0.42-0.115*L-0.069*L**2])
        al0, alp = (0.22+0.17*L, 0.15)
        alt = al0 + alp*t - 1
        xs = array([self._sea(x, eta, alt, j, glue=-1) for j in range(3)])
        b = 2.58 + 0.25*log(Mp2/(Q2+Mp2))
        return exp(b*t) * sum(cs*xs)

    #   ===  GPD E ===

    def Euval(self, x, eta, t, Q2):
        """GK12 model for E^u_val GPD."""
        ceuval = array([1, -1])  # (1-rho)
        al0, alp = (0.48, 0.9)
        alt = al0 + alp*t
        # val() expects series in beta^(j/2), so 0,2 -> 0,1
        xs = array([self._val(x, eta, alt, j) for j in [0, 2]])
        # prefactor to get agreement by PARTONS' GPDGK11.cpp
        return 1.67013027*npsum(ceuval*xs)/beta(1-al0, 5)

    def Edval(self, x, eta, t, Q2, DMbeta=False):
        """GK12 model for E^d_val GPD."""
        # (1-rho)^2.6 up to rho^8
        if DMbeta:
            # Choice of DM in his notebook
            betd = 6
            cedval = [1, -3, +3, -1]
            bmax = 3
        else:
            # GK choice
            betd = 5.6
            # PARTONS don't use the last term below
            cedval = array([1, -2.6, 2.08, -0.416, -0.0416, -0.011648, -0.0046592,
                            -0.00226304, -0.001244672])
            bmax = 8
        al0, alp = (0.48, 0.9)
        alt = al0 + alp*t
        # val() expects series in beta^(j/2)
        xs = array([self._val(x, eta, alt, j) for j in range(0, 2*bmax+1, 2)])
        # -2.03 -> 2.029 for better agreement with PARTONS
        return -2.029*npsum(cedval*xs)/beta(1-al0, 1+betd)

    def Esea(self, x, eta, t, Q2):
        """GK12 model for strange sea GPD."""
        Q02 = 4.
        Mp2 = 0.9383**2
        L = log(Q2/Q02)
        ces = array([1, -2, 1])
        al0, alp = (1.10+0.06*L-0.0027*L**2, 0.15)
        alt = al0 + alp*t
        xs = array([self._sea(x, eta, alt, j) for j in range(0, 5, 2)])
        b = 2.58 + 0.25*log(Mp2/(Q2+Mp2))
        # Note that opposite sign is also possible:
        return -0.155*exp(b*t) * npsum(ces*xs)

    def Eu(self, x, eta, t, Q2):
        """Up quark GPD E."""
        return self.Euval(x, eta, t, Q2) + self.Esea(x, eta, t, Q2)

    def Ed(self, x, eta, t, Q2):
        """Down quark GPD E."""
        return self.Edval(x, eta, t, Q2) + self.Esea(x, eta, t, Q2)

    def Es(self, x, eta, t, Q2):
        """Strange quark GPD E."""
        return self.Esea(x, eta, t, Q2)

    def Eg(self, x, eta, t, Q2):
        """GK12 model for gluon GPD E."""
        Q02 = 4.
        Mp2 = 0.9383**2
        L = log(Q2/Q02)
        al0, alp = (1.10+0.06*L-0.0027*L**2, 0.15)
        alt = al0 + alp*t - 1
        ces = array([1, -1])
        bmax = 2
        xs = array([self._sea(x, eta, alt, j, glue=True) for j in [0, 2]])
        b = (2.58 + 0.25*log(Mp2/(Q2+Mp2)))
        return 0.779*exp(b*t) * sum(ces*xs)

    #   ===  GPD E tilde ===

    def Etuval(self, x, eta, t, Q2):
        """GK12 model for Et^u GPD (non-pole part)."""
        ces = array([1, -2, 1])
        al0, alp = (0.48, 0.45)
        alt = al0 + alp*t
        xs = array([self._val(x, eta, alt, j) for j in range(0, 5, 2)])
        b = 0.9
        return 14.*exp(b*t) * npsum(ces*xs)

    def Etdval(self, x, eta, t, Q2):
        """GK12 model for Et^d GPD (non-pole part)."""
        ces = array([1, -2, 1])
        al0, alp = (0.48, 0.45)
        alt = al0 + alp*t
        xs = array([self._val(x, eta, alt, j) for j in range(0, 5, 2)])
        b = 0.9
        return 4.*exp(b*t) * npsum(ces*xs)

    def Etupole(self, x, eta, t, Q2, cff=False):
        """Implementantion following PARTONS GPDGK16."""
        # Note that only asymptotic part of pion DA is used
        # FIXME: if-then logic not nice below
        res = 0
        
        Mp = 0.938272013
        mpi2 = 0.1349766**2
        gpiNN = 13.4
        f_pi = 0.131
        Lambda_N2 = 0.51 ** 2

        xbj = 2.0 * eta / (eta - t/Q2 * (eta-0.5) + 1.0)
        eps = 2.0 * xbj * Mp / sqrt(Q2)
        eps2 = eps ** 2

        if eps < 1 and (4.0 * xbj * (1.0 - xbj) + eps2) != 0:
            tmin = -Q2 * (2.0 * (1.0 - xbj) * (1 - sqrt(1.0 + eps2)) + eps2) \
                   / (4.0 * xbj * (1.0 - xbj) + eps2)
            FpiNN = (Lambda_N2 - mpi2) / (Lambda_N2 - (t - tmin))
            Fp = -Mp * f_pi * (2.0 * sqrt(2.0) * gpiNN * FpiNN) \
                 / (t - mpi2)
            if t < tmin: 
                res = (Fp / 4.0 / eta) * 6.0 
        if cff:
            # already integrated over LO hard-scattering amplitude
            return res
        else:
            if x < eta and x > -eta:
                y = (x + eta) / (2.0 * eta)
                res = res * y * (1.0 - y)
            else:
                res = 0
            return res 

    def Etdpole(self, x, eta, t, Q2):
        """Pion pole contribution to E tilde - d quark."""
        return - self.Etupole(x, eta, t, Q2)

    def Etu(self, x, eta, t, Q2): 
        return self.Etuval(x, eta, t, Q2) + self.Etupole(x, eta, t, Q2)
        
    def Etd(self, x, eta, t, Q2): 
        return self.Etdval(x, eta, t, Q2) + self.Etdpole(x, eta, t, Q2)

    def Etg(self, x, eta, t, Q2):
        """GK12 model for gluon GPD E-tilde."""
        return 0.


    #  ######  CFFs  #########

    def ImH(self, pt, xi=0):
        """GK12 model for Im(CFFH)."""
        # FIXME: The following solution is not elegant
        if isinstance(xi, ndarray):
            # function was called with third argument that is xi nd array
            xs = xi
        elif xi != 0:
            # function was called with third argument that is xi number
            xs = [xi]
        else:
            # xi should be taken from pt object
            xs = [pt.xi]
        t = pt.t
        Q2 = pt.Q2
        res = []
        for x in xs:
            res.append(pi*((4./9)*(self.Hu(x, x, t, Q2) - self.Hu(-x, x, t, Q2))
                       + (1./9)*(self.Hd(x, x, t, Q2) - self.Hd(-x, x, t, Q2))
                       + (1./9)*(self.Hs(x, x, t, Q2) - self.Hs(-x, x, t, Q2))))
        if len(res) == 1:
            return res[0]
        else:
            return array(res)


    def ImHt(self, pt, xi=0):
        """GK12 model for Im(CFFHt)."""
        # FIXME: The following solution is not elegant
        if isinstance(xi, ndarray):
            # function was called with third argument that is xi nd array
            xs = xi
        elif xi != 0:
            # function was called with third argument that is xi number
            xs = [xi]
        else:
            # xi should be taken from pt object
            xs = [pt.xi]
        t = pt.t
        Q2 = pt.Q2
        res = []
        for x in xs:
            assert x > 0
            res.append(pi*((4./9)*self.Htuval(x, x, t, Q2)
                       + (1./9)*self.Htdval(x, x, t, Q2)))
        if len(res) == 1:
            return res[0]
        else:
            return array(res)

    def ImE(self, pt, xi=0, DMbeta=False):
        """GK12 model for Im(CFFE)."""
        # FIXME: The following solution is not elegant
        if isinstance(xi, ndarray):
            # function was called with third argument that is xi nd array
            xs = xi
        elif xi != 0:
            # function was called with third argument that is xi number
            xs = [xi]
        else:
            # xi should be taken from pt object
            xs = [pt.xi]
        t = pt.t
        Q2 = pt.Q2
        res = []
        for x in xs:
            assert x > 0
            res.append(pi*(
                (4./9)*(self.Euval(x, x, t, Q2) + 2*self.Esea(x, x, t, Q2))
                + (1./9)*(self.Edval(x, x, t, Q2, DMbeta)+2*self.Esea(x, x, t, Q2))
                + (1./9)*(2*self.Esea(x, x, t, Q2))))
        if len(res) == 1:
            return res[0]
        else:
            return array(res)

    def ImEt(self, pt, xi=0):
        """GK12 model for Im(CFFEt)."""
        # FIXME: The following solution is not elegant
        if isinstance(xi, ndarray):
            # function was called with third argument that is xi nd array
            xs = xi
        elif xi != 0:
            # function was called with third argument that is xi number
            xs = [xi]
        else:
            # xi should be taken from pt object
            xs = [pt.xi]
        t = pt.t
        Q2 = pt.Q2
        res = []
        for x in xs:
            assert x > 0
            res.append(pi*((4./9)*(self.Etuval(x, x, t, Q2))
                       + (1./9)*(self.Etdval(x, x, t, Q2))))
        if len(res) == 1:
            return res[0]
        else:
            return array(res)

    def ReEtnonpole(self, pt):
        """Contribution to ReEt extra to pion pole."""
        return cff.DispersionCFF.ReEt(self, pt)

    def ReEt(self, pt):
        """Total ReEt pole+nonpole."""
        # (4/9 Etu_pole + 1/9 Etdpole) = 1/3 Etu_pole
        return self.ReEtnonpole(pt) + self.Etupole(pt.xi, pt.xi, pt.t, pt.Q2, cff=True)/3.

    def _GKaux(self, x, ds, Q2=4.):
        """Auxilliary integrand that accepts ndarray."""
        res = []
        for d in ds:
            res.append(self.Hu(x, 0, -d**2, Q2))
        return res

    def _GKauxE(self, x, ds, Q2=4.):
        """Auxilliary integrand that accepts ndarray."""
        res = []
        for d in ds:
            res.append(self.Eu(x, 0, -d**2, Q2))
        return res

    def gpdb(self, x, b, Q2=4., inf=2.5):
        """GPD in b space ??? FIXME."""
        return bquadrature(
                lambda d: d*j0(b*d/GeVfm)*self._GKaux(x, d, Q2),
                0.0, inf) / (2*pi*GeVfm**2)

    def gpdbpol(self, pol, x, bx, by, Q2=4., inf=2.5):
        """Polarized GPD in b space."""
        b = sqrt(bx**2 + by**2)
        aux = bquadrature(
                lambda d: d**2*j1(b*d/GeVfm)*self._GKauxE(x, d, Q2),
                0.0, inf) / (2*pi*GeVfm**2)
        return self.gpdb(x, b, Q2, inf) + pol*bx*aux/(2*b*Mp)


class GK12D(GoloskokovKrollCFF):
    """Goloskokov-Kroll GPD/CFF model with added D-term.

    From [arXiv:1210.6975], Kroll:2012sm
    Real part of CFFs is obtained using dispersion relations.

    """
    def subtraction(self, pt):
        """D-term."""
        denom = (1 - pt.t/(0.841*0.487**2))**0.841
        return -(10./9.)*(-1.9/denom)
