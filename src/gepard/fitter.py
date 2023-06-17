"""Implementation of fitting algorithms."""

from __future__ import annotations

import torch
from typing import List
from iminuit import Minuit  # type: ignore


from . import data, theory


class Fitter(object):
    """Superclass for fitting procedures/algorithms."""

    def __init__(self, **kwargs) -> None:
        for key in kwargs:
            setattr(self, key, kwargs[key])


class MinuitFitter(Fitter):
    """Fits using iminuit Python frontend to MINUIT2 C++ library.

    Args:
        fitpoints (data.DataSet): points to fit to
        theory (theory.Theory): theory/model to be fitted

    """
    def __init__(self, fitpoints: data.DataSet,
                 theory: theory.Theory, **kwargs) -> None:
        self.fitpoints = fitpoints
        self.theory = theory
        init_vals = [v for v in self.theory.parameters.values()]
        names = [k for k in self.theory.parameters.keys()]

        def fcn(p):
            """Cost function for minimization - chi-square."""
            # Update model parameters.
            # Note: This relies on Python>3.7 where dicts are ordered!!
            self.theory.parameters.update(zip(self.theory.parameters.keys(), p))
            chisq = self.theory.chisq(self.fitpoints)
            return chisq

        fcn.errordef = Minuit.LEAST_SQUARES

        self.minuit = Minuit(fcn, init_vals, name=names)
        for k, v in self.theory.parameters_fixed.items():
            self.minuit.fixed[k] = v
        for k, v in self.theory.parameters_limits.items():
            self.minuit.limits[k] = v
        Fitter.__init__(self, **kwargs)

    def covsync(self):
        """Synchronize covariance and theory parameter errors with iminuit."""
        self.theory.parameters_errors = self.minuit.errors.to_dict()
        # iminuit gives covariance table for all parameters, here
        # we take only free ones:
        try:
            self.theory.covariance = {(p1, p2): self.minuit.covariance[p1, p2]
                                      for p1 in self.theory.free_parameters()
                                      for p2 in self.theory.free_parameters()}
        except AttributeError:
            print("Something's problematic. No covariances available.")

    def fit(self):
        """Perform simple fit.

        For better control, use iminuit's functions.

        """
        self.minuit.migrad()
        self.covsync()

    # The following methods keep status of parameters (fixed, limits)
    # in sync between Theory object and minuit object

    def fix_parameters(self, *args):  # noqa: D402
        """fix_parameters('p1', 'p2', ...)."""
        if args[0] == 'ALL':
            # fix 'em all
            self.theory._fix_parameters('ALL')
            for par in self.theory.model.parameters:
                self.minuit.fixed[par] = True
        else:
            self.theory._fix_parameters(*args)
            for par in args:
                self.minuit.fixed[par] = True

    def release_parameters(self, *args):  # noqa: D402
        """release_parameters('p1', 'p2', ...)."""
        self.theory._release_parameters(*args)
        for par in args:
            self.minuit.fixed[par] = False

    def limit_parameters(self, dct):  # noqa: D402
        """limit_parameters({'p1': (lo, hi), ...}."""
        self.theory.parameters_limits.update(dct)
        for k, v in dct.items():
            self.minuit.limits[k] = v

    def free_parameters(self) -> List[str]:
        """Return list of names of free fitting parameters."""
        self.theory.free_parameters()

    def print_parameters(self):
        """Values and errors for free parameters."""
        self.theory.print_parameters()


class CustomLoss(torch.nn.Module):
    def __init__(self, fitpoints, theory):
        self.fitpoints = fitpoints
        self.theory = theory
        super(CustomLoss, self).__init__()

    def forward(self, cff_pred, obs_true):
        # Our custom loss function is essentially chi-squared
        preds = []
        for cffs, id in zip(cff_pred, obs_true[:, -1]):
            pt = self.fitpoints[int(id)]
            preds.append(self.theory.predict_while_train(cffs, pt))
        preds = torch.stack(preds)
        # print("Loss preds = {}".format(preds))
        return torch.mean(torch.square((preds - obs_true[:, 0])/obs_true[:, 1]))

class NeuralFitter(Fitter):
    """Fits using PyTorch library.

    Args:
        fitpoints (data.DataSet): points to fit to
        theory (theory.Theory): theory/model to be fitted

    """
    def __init__(self, fitpoints: data.DataSet,
                 theory: theory.Theory, **kwargs) -> None:
        self.fitpoints = fitpoints
        self.theory = theory
        x_train = []
        y_train = []
        for k, pt in enumerate(fitpoints):
            pt.ptid = k
            pt.phi = torch.tensor(pt.phi)
            x_train.append([pt.xB, pt.t])
            y_train.append([pt.val, pt.err, k])
        # We pass the point index k through the NNet (as irellevant feature!) so that Gepard
        # can determine which DataPoint to use to calculate observable.
        # TODO: Explain to Torch that this feature is irellevant! (This likely mean using
        #        some more involved specific net architecture.
        self.x_train = torch.tensor(x_train, dtype=torch.float32)
        self.y_train = torch.tensor(y_train, dtype=torch.float32)
        self.criterion = CustomLoss(fitpoints, theory)
        self.optimizer = torch.optim.Adam(self.theory.nn_model.parameters(), lr=0.01)
        self.history = []
    

    def train(self, epochs=1000):
        
        for epoch in range(epochs):
            self.optimizer.zero_grad()
            cff_pred = self.theory.nn_model(self.x_train)
            self.loss = self.criterion(cff_pred, self.y_train)
            self.history.append(float(self.loss))
            self.loss.backward()
            self.optimizer.step()

        print("Final loss = {}".format(self.history[-1]))

