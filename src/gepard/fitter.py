"""Implementation of fitting algorithms."""

from __future__ import annotations

import numpy as np
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

def data_replica(datapoints, train_percentage=70, smear_replicas=True):
    """Returns datapoints split into training and testing set.

    Args:
        datapoints (data.DataSet): points that are to be replicated
        smear_replicas (Bool): Should we smear the replica according the uncertainties

    Note that datapoints get id numbers pt.ptid needed for
    keeping link to proper datapoint for the loss calculation.

    TODO: pt.err should be removed from y_train and loss, we should NOT
      weigh the point twice, once when creating replicas and
      second time when we train to it!

    """
    x_train = []
    y_train = []
    x_test = []
    y_test = []
    train_size = int(len(datapoints) * train_percentage / 100)
    for k, pt in enumerate(np.random.permutation(datapoints)):
        # This rounding was because pt.val was used as pt ID in old pybrain version
        # y = pt.val + round(np.random.normal(0, pt.err, 1)[0], 5)
        if smear_replicas:
            y = pt.val + np.random.normal(0, pt.err, 1)[0]
        else:
            y = pt.val
        if k < train_size:
            x_train.append([pt.xB, pt.t])
            y_train.append([y, pt.err, pt.ptid])
        else:
            x_test.append([pt.xB, pt.t])
            y_test.append([y, pt.err, pt.ptid])
    # We pass the pt.ptid to the loss function so Gepard
    # can use it to determine which DataPoint to use to calculate observable.
    x_train = torch.tensor(x_train, dtype=torch.float32)
    y_train = torch.tensor(y_train, dtype=torch.float32)
    x_test = torch.tensor(x_test, dtype=torch.float32)
    y_test = torch.tensor(y_test, dtype=torch.float32)
    return x_train, y_train, x_test, y_test

def datasets_replica(datasets, train_percentage=70, smear_replicas=True, q2in=False):
    """TODO: merge with data_replica / a lot of code duplication."""
    x_train = []
    y_train = []
    x_test = []
    y_test = []
    # We make train/test separation on level of datasets because our datasets
    # can have very different sizes
    for dataset in datasets:
        train_size = int(len(dataset) * train_percentage / 100)
        for k, pt in enumerate(np.random.permutation(dataset)):
            if smear_replicas:
                y = pt.val + np.random.normal(0, pt.err, 1)[0]
            else:
                y = pt.val
            if k < train_size:
                if q2in:
                    x_train.append([pt.xB, pt.t, pt.Q2])
                else:
                    x_train.append([pt.xB, pt.t])
                y_train.append([y, pt.err, pt.ptid])
            else:
                if q2in:
                    x_test.append([pt.xB, pt.t, pt.Q2])
                else:
                    x_test.append([pt.xB, pt.t])
                y_test.append([y, pt.err, pt.ptid])
    # We pass the pt.ptid to the loss function so Gepard
    # can use it to determine which DataPoint to use to calculate observable.
    x_train = torch.tensor(x_train, dtype=torch.float32)
    y_train = torch.tensor(y_train, dtype=torch.float32)
    x_test = torch.tensor(x_test, dtype=torch.float32)
    y_test = torch.tensor(y_test, dtype=torch.float32)
    return x_train, y_train, x_test, y_test

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
        # FIXME: Hard-wired factor 10 tuned to CLAS observables y range. Kludge.
        # return 10*torch.mean(torch.square((torch.stack(preds) - obs_true[:, 0]))/obs_true[:, 1])
        # Standard chisq. Do we doubly penalize the large-uncertainty points?
        return torch.mean(torch.square((torch.stack(preds) - obs_true[:, 0])/obs_true[:, 1]))


class NeuralFitter(Fitter):
    """Fits using PyTorch library.

    Args:
        fitpoints (data.DataSet): points to fit to
        theory (theory.Theory): theory/model to be fitted
        smear_replicas (Bool): Should we smear the replica according the uncertainties
        patience (Int): How many batches to wait for improvement before early stopping

    """
    def __init__(self, fitpoints: data.DataSet,
                 theory: theory.Theory, **kwargs) -> None:
        if type(fitpoints) == list:
            self.datasets = fitpoints
            self.datapoints = fitpoints[0][:]
            for dataset in fitpoints[1:]:
                self.datapoints += dataset
        else:
            self.datapoints = fitpoints
            self.datasets = [fitpoints]
        for k, pt in enumerate(self.datapoints):
            pt.ptid = k
        self.theory = theory
        self.nnets = 4
        self.maxtries = 999
        self.nbatch = 20
        self.batchlen = 5
        self.patience = 5
        self.lx_lambda = 0.001
        self.regularization = None
        self.smear_replicas = True
        self.train_percentage = 70
        self.criterion = CustomLoss(self.datapoints, theory)
        Fitter.__init__(self, **kwargs)

    def train_net(self, datasets, new_net=True):
        """Create net trained on datapoints.

        Args:
            new_net (bool): should we create new net or continue
                            with existing self.theory.nn_model
        """
        self.theory.in_training = True
        x_train, y_train, x_test, y_test = datasets_replica(datasets,
				train_percentage=self.train_percentage, smear_replicas=self.smear_replicas,
                q2in = self.theory.q2in)
        self.theory.nn_mean, self.theory.nn_std = self.theory.get_standard(x_train)
        x_train_standardized = self.theory.standardize(x_train,
                self.theory.nn_mean, self.theory.nn_std)
        x_test_standardized = self.theory.standardize(x_test,
                self.theory.nn_mean, self.theory.nn_std)
        if new_net:
            self.theory.nn_model, self.optimizer = self.theory.build_net()
            self.history = []
        mem_state_dict = self.theory.nn_model.state_dict()
        mem_err = 100  # large init error, almost certain to be bettered
        early_stop_counter = 1
        for k in range(1, self.nbatch+1):
            for epoch in range(self.batchlen):
                self.optimizer.zero_grad()
                cff_pred = self.theory.nn_model(x_train_standardized)
                loss = self.criterion(cff_pred, y_train)
                if self.regularization == 'L1':
                    lx_norm = sum(torch.linalg.norm(p, 1)
                            for p in self.theory.nn_model.parameters())
                elif self.regularization == 'L2':
                    lx_norm = sum(torch.linalg.norm(p, 2) ** 2
                            for p in self.theory.nn_model.parameters())
                else:
                    lx_norm = 0
                loss = loss + self.lx_lambda * lx_norm
                self.history.append(float(loss))
                loss.backward()
                self.optimizer.step()
            test_cff_pred = self.theory.nn_model(x_test_standardized)
            test_loss = float(self.criterion(test_cff_pred, y_test))
            print("\nEpoch {:3d}: train error = {:.4f} test error = {:.4f} ".format(
                k*self.batchlen, self.history[-1], test_loss), end='')
            if test_loss < mem_err:
                mem_state_dict = self.theory.nn_model.state_dict()
                mem_err = test_loss
                print('-', end='')
                early_stop_counter = 1
            elif test_loss > 100:
                print("\nHopeless. Giving up")
                mem_err = -1  # fail flag
                break
            else:
                if early_stop_counter == self.patience:
                    print("+\nNo improvement for {} batches. Stopping early.".format(self.patience))
                    break
                else:
                    print('+', end='')
                    early_stop_counter += 1
        self.theory.nn_model.load_state_dict(mem_state_dict)
        self.theory.in_training = False
        return self.theory.nn_model, self.theory.nn_mean, self.theory.nn_std, mem_err


    def fit(self):
        """Train number (nnet) of nets."""
        for n in range(self.nnets):
            test_err = -1
            while test_err < 0:
                net, mean, std, test_err = self.train_net(self.datasets)
            self.theory.nets.append((net, mean, std))
            print("Net {} --> test_err = {}".format(n, test_err))


    def fitgood(self, max_test_err : float = 2):
        """Train until you have nnet good nets.

        CAREFUL: you are in danger of overfitting if you use this,
                 stick with normal fit unless you know what you are doing!
        """
        n = 0
        k = 0
        while n < self.nnets and k < self.maxtries:
            k += 1
            net, mean, std, test_err = self.train_net(self.datasets)
            self.theory.nets.append((net, mean, std))
            if test_err > max_test_err:
                del self.theory.nets[-1]
            else:
                print("test_err = {} < max_test_err = {} so we accept the net".format(
                    test_err, max_test_err))
                n += 1
            print("[Try {}/{}] {} good nets. Last test_err = {}".format(
                k, self.maxtries, n, test_err))
            if (k > self.maxtries/4) and (n < 2):
                print("Less than 2 good nets found after 25% of maxtries. Giving up.")
                break
