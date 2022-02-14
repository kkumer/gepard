"""Models for Elastic (Electromagnetic) Form Factors."""


class EFF(object):
    """Dirac and Pauli elastic (electromagnetic) form factors F_1 and F_2."""


class DipoleEFF(EFF):
    """Dipole approximation from Dieter's notebook."""

    def F1(self, pt):
        """Dirac elastic proton form factor - dipole parametrization."""
        t = pt.t
        if 'in2particle' in pt and pt.in2particle == 'p':
            return (1.41 * (1.26 - t))/((0.71 - t)**2 * (3.53 - t))
        else:
            raise Exception(
                    'Neutron dipole elastic FFs are not implemented yet! Use Kelly.')

    def F2(self, pt):
        """Pauli elastic proton form factor - dipole parametrization."""
        t = pt.t
        if 'in2particle' in pt and pt.in2particle == 'p':
            return 3.2 / ((0.71 - t)**2 * (3.53 - t))
        else:
            raise Exception(
                    'Neutron dipole elastic FFs are not implemented yet! Use Kelly.')


class KellyEFF(EFF):
    """Kelly's approximation from DM's notebook."""

    def F1(self, pt):
        """Dirac elastic nucleon form factor - Kelly's parametrization."""
        if 'in2particle' in pt and pt.in2particle == 'n':
            return self._nF1(pt.t)
        else:   # proton is default
            return self._pF1(pt.t)

    def F2(self, pt):
        """Dirac elastic nucleon form factor - Kelly's parametrization."""
        if 'in2particle' in pt and pt.in2particle == 'n':
            return self._nF2(pt.t)
        else:   # proton is default
            return self._pF2(pt.t)

    def _pF1(self, t):
        """Dirac elastic proton form factor - Kelly's parametrization."""
        return ((1 + 0.06815437285120148*t)/(1 - 3.118062557942468*t +
                1.0338391956016382*t**2 - 0.5031268669574522*t**3) -
                (0.7931031653189349*(1 - 0.03407718642560074*t)*t) /
                (1 - 3.115222792407001*t + 1.520921000705686*t**2 -
                0.14999913420898098*t**3))/(1 - 0.28397655354667284*t)

    def _pF2(self, t):
        """Pauli elastic proton form factor - Kelly's parametrization."""
        return (-((1 + 0.06815437285120148*t)/(1 - 3.118062557942468*t +
                1.0338391956016382*t**2 - 0.5031268669574522*t**3)) +
                (2.792847351*(1 - 0.03407718642560074*t))/(1 -
                3.115222792407001*t + 1.520921000705686*t**2 -
                0.14999913420898098*t**3)) / (1 - 0.28397655354667284*t)

    def _nF1(self, t):
        """Dirac elastic neutron form factor - Kelly's parametrization."""
        return ((-0.4842637275288574*t)/((1 - 1.4084507042253522*t)**2 *
                (1 - 0.9345440355820054*t)) +
                (0.5417644379086957*(1 - 0.6598447281533554*t)*t) /
                (1 - 4.168632789020339*t + 1.9408278987597791*t**2 -
                1.9100884849907935*t**3))/(1 - 0.2831951622975774*t)

    def _nF2(self, t):
        """Pauli elastic neutron form factor - Kelly's parametrization."""
        return ((0.4842637275288574*t)/((1 - 1.4084507042253522*t)**2 *
                (1 - 0.9345440355820054*t)) - (1.9130427*(1 - 0.6598447281533554*t)) /
                (1 - 4.168632789020339*t + 1.9408278987597791*t**2 -
                1.9100884849907935*t**3))/(1 - 0.2831951622975774*t)


class ZeroEFF(EFF):
    """Set F1=F2=0 to get just DVCS^2."""

    def F1(self, pt):
        """Elastic em Dirac form factor F1 set to zero."""
        return 0.

    def F2(self, pt):
        """Elastic em Pauli form factor F2 set to zero."""
        return 0.
