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

    From [arXiv:1210.6975], Kroll:2012sm
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

    def _sea(self, x, eta, alt, j):
        DMKILL = 1  # Set to 0 to agree with DM
        assert eta >= 0
        if x >= eta:
            return self._intsea(x, eta, alt, j)
        elif -eta < x < eta:
            return (self._intsea(self, x, eta, alt, j, zero=True)
                    - DMKILL*self._intsea(-x, eta, alt, j, zero=True))
        else:
            return -DMKILL*self._intsea(-x, eta, alt, j)  # x < -eta

    #   ######  GPDs  #########

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
        """GK12 model for strange sea GPD."""
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
        """GK12 model for u or d sea GPDs."""
        Q02 = 4.
        kaps = 1.+0.68/(1.+0.52*log(Q2/Q02))
        return kaps * self.Hs(x, eta, t, Q2)

    def Hu(self, x, eta, t, Q2):
        """Up quark GPD."""
        return self.Huval(x, eta, t, Q2) + self.Hudsea(x, eta, t, Q2)

    def Hd(self, x, eta, t, Q2):
        """Down quark GPD."""
        return self.Hdval(x, eta, t, Q2) + self.Hudsea(x, eta, t, Q2)

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

    def Euval(self, x, eta, t, Q2):
        """GK12 model for E^u_val GPD."""
        ceuval = array([1, -1])  # (1-rho)
        al0, alp = (0.48, 0.9)
        alt = al0 + alp*t
        # val() expects series in beta^(j/2), so 0,2 -> 0,1
        xs = array([self._val(x, eta, alt, j) for j in [0, 2]])
        return 1.67*npsum(ceuval*xs)/beta(1-al0, 5)

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
            cedval = array([1, -2.6, 2.08, -0.416, -0.0416, -0.011648, -0.0046592,
                            -0.00226304, -0.001244672])
            bmax = 8
        al0, alp = (0.48, 0.9)
        alt = al0 + alp*t
        # val() expects series in beta^(j/2)
        xs = array([self._val(x, eta, alt, j) for j in range(0, 2*bmax+1, 2)])
        return -2.03*npsum(cedval*xs)/beta(1-al0, 1+betd)

    def Esea(self, x, eta, t, Q2):
        """GK12 model for strange sea GPD."""
        Q02 = 4.
        Mp2 = 0.9383**2
        L = log(Q2/Q02)
        ces = array([1, -2, 1])
        al0, alp = (1.10+0.06*L-0.0027*L**2, 0.15)
        alt = al0 + alp*t
        xs = array([self._sea(x, eta, alt, j) for j in range(0, 5, 2)])
        b = 0.9*(2.58 + 0.25*log(Mp2/(Q2+Mp2)))
        # Note that opposite sign is also possible:
        return -0.155*exp(b*t) * npsum(ces*xs)

    def Eu(self, x, eta, t, Q2):
        """Up quark GPD E."""
        return self.Euval(x, eta, t, Q2) + self.Esea(x, eta, t, Q2)

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

    def _gegen32(self, x):
        return 15.*x**2/2 - 3./2

    def _phi(self, u, a2=0.22):
        return 6*u*(1-u)*(1+a2*self._gegen32(2*u+1))

    def Etupole(self, x, eta, t, Q2):
        """Pion pole contribution to E tilde - u quark."""
        if -eta < x < eta:
            Mp2 = 0.938272013**2
            gA0 = 1.267
            LAM2 = 0.44**2
            mpi2 = 0.1396**2
            numer = Mp2 * gA0 * (LAM2 - mpi2)
            denom = (mpi2-t)*(LAM2-t)
            return numer*self._phi((x+eta)/2/eta)/denom/eta
        else:
            return 0

    def Etdpole(self, x, eta, t, Q2):
        """Pion pole contribution to E tilde - d quark."""
        return - self.Etupole(x, eta, t, Q2)

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

    # def ImHalt(self, pt):
        # """GK12 model for Im(CFFH) (alternative expression)."""
        # x = pt.xi
        # t = pt.t
        # Q2 = pt.Q2
        # assert x>0
        # return pi*(  (4./9)*(self.Huval(x, x, t, Q2) +2*self.Hudsea(x, x, t, Q2))
        #       # +(1./9)*(self.Hdval(x, x, t, Q2)+2*self.Hudsea(x, x, t, Q2))
        #       # +(1./9)*(2*self.Hs(x, x, t, Q2)) )

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

    def ReEtpole(self, pt):
        """Pion pole contribution to ReEt."""
        Mp2 = 0.938272013**2
        gA0 = 1.267
        LAM2 = 0.44**2
        mpi2 = 0.1396**2
        numer = Mp2 * gA0 * (LAM2 - mpi2)
        denom = (mpi2-pt.t)*(LAM2-pt.t)
        green = numer/denom
        a2 = 0.22
        return 2*(1+a2)*green/pt.xi

    def ReEtnonpole(self, pt):
        """Contribution to ReEt extra to pion pole."""
        return cff.DispersionCFF.ReEt(self, pt)

    def ReEt(self, pt):
        """Total ReEt pole+nonpole."""
        return self.ReEtpole(pt) + self.ReEtnonpole(pt)

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
