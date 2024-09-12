"""Distribution amplitude (DA) models.

   Module is structured similary to gpd.py
   but is much simpler. Only conformal-moment modelling
   is implemented at the moment. Scheme is always msbar.
   (It is the responsibility of user not to combine these DAs
   with GPDs in csbar scheme.)

"""

import numpy as np
from . import model

class DA(model.ParameterModel):
    """Base class of all DA models.

    Args:
        p: pQCD order (0 = LO, 1 = NLO, 2 = NNLO)
        nf: number of active quark flavors (default is 4)
        Q02: Initial Q0^2 for DA evolution (default is 4 GeV^2)
        r20: Initial mu0^2 for alpha_strong definition.
        asp: alpha_strong/(2*pi) at scale r20 for (LO, NLO, NNLO)

    """
    def __init__(self, **kwargs) -> None:
        self.p = kwargs.setdefault('p', 0)
        self.scheme = kwargs.setdefault('scheme', 'msbar')
        self.nf = kwargs.setdefault('nf', 4)
        self.Q02 = kwargs.setdefault('Q02', 4.0)
        self.r20 = kwargs.setdefault('r20', 2.5)
        self.asp = kwargs.setdefault('asp', np.array([0.0606, 0.0518, 0.0488]))
        # scales
        self.rr2 = 1     # ratio of Q2/renorm. scale squared
        self.rdaf2 = 1   # ratio of Q2/DA fact. scale sq.
        #
        super().__init__(**kwargs)


class GegenbauerDA(DA):
    """Base class of DA models built in conformal moment space.

    Args:
        ngegens: number of conformal/Gegenbauer moments
                 including the first asymptotic a0 = 1.

    Notes:
        Parameters of this model must have names of form
        a2, a4, ..., a<2*(ngegens-1)>. a0 is fixed to be 1
        and is not a model parameter.

    """

    def __init__(self, **kwargs) -> None:
        self.ngegens = 3
        # gpoints is array of integer DA moments [0, 2, 4, ...] and
        # corresponds to complex jpoints on GPD MB contour
        self.gpoints = np.arange(0, 2*self.ngegens, 2)
        # Initial parameters correspond to asymptotic DA.
        self.add_parameters({'a2': 0., 'a4': 0.}) 
        super().__init__(**kwargs)

    def gegenbauers(self) -> np.ndarray:
        """Values of DA Gegenbauer moments."""

        gegens = [1]  # a0 = 1
        for g in self.gpoints[1:]:
            gegens.append(self.parameters['a{:d}'.format(g)])
        return np.array(gegens)

