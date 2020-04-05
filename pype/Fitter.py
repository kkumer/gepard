"""Classes for fitting."""
#from IPython.Debugger import Tracer; debug_here = Tracer()

import sys, warnings
import logzero
_lg = logzero.logger

import numpy as np

# try:
    # from minuit import Minuit
# except:
    # pass

try:
    from iminuit import Minuit, InitialParamWarning
    # FIXME: maybe we should not switch this off so bluntly
    warnings.simplefilter('ignore', InitialParamWarning, append=False)
except:
    pass


from pybrain3.tools.shortcuts import buildNetwork
import brain
import trans  # output layer transformation for FitterBrain
import utils


class Fitter(object):
    """Superclass for fitting procedures/algorithms."""

    def __init__(self, **kwargs):
        for key in kwargs:
            setattr(self, key, kwargs[key])
        _lg.info('Fitter instance created')


class FitterMinuit(Fitter):
    """Fits using pyminuit."""

    def __init__(self, fitpoints, theory, **kwargs):
        self.fitpoints = fitpoints
        self.theory = theory
        self.printMode = 0

        def fcn(p):
            """Cost function for minimization - chi-square."""
            for n, pname in enumerate(theory.model.parameter_names):
                # FIXME: only non-fixed should be treated here; or total np array
                theory.model.parameters[pname] = p[n]
            if 'Gepard' in theory.model.__dict__:
                theory.model.ndparameters.put(range(theory.model.DR.ndparameters.size), p[:50])
            else:
                theory.model.ndparameters.put(range(theory.model.ndparameters.size), p)
            chisq = 0
            for pt in fitpoints:
                chisq += (getattr(theory, pt.yaxis)(pt) - pt.val)**2 / pt.err**2
            _lg.info('Minuit: {:4d} calls --> chisq/npt = {:.3f}/{}'.format(
                self.minuit.get_num_call_fcn()+1, chisq, len(fitpoints)))
            return chisq

        if isinstance(theory.model.parameters, utils.hubDictNew):
            # This is needed because in Python <=2.5 ** operator
            # requires dict as an argument, i.e. my hubDictNew wouldn't work:
            auxdict = dict((it for it in list(theory.model.parameters.items())))
            self.minuit = Minuit(fcn, use_array_call=True, errordef=1,
                    forced_parameters=theory.model.parameter_names, **auxdict)
        else:
            self.minuit = Minuit(fcn, use_array_call=True, errordef=1,
                    forced_parameters=theory.model.parameter_names, **theory.model.parameters)
        for key in kwargs:
            setattr(self.minuit, key, kwargs[key])
        Fitter.__init__(self, **kwargs)


    def fit(self):
        self.minuit.migrad()
        if self.printMode > 0:
            print("ncalls = \n", self.minuit.ncalls)
            self.theory.print_chisq(self.fitpoints)
        # Set/update covariance matrix of model:
        self.theory.model.covariance = self.minuit.covariance
        if self.printMode > 0:
            print("")
            self.theory.model.print_parameters_errors(pvalues=True,
                    ndof=len(self.fitpoints)-utils.npars(self.theory.model))

    def fix_parameters(self, *args):
        """fix_parameters('p1', 'p2', ...)."""

        if args[0] == 'ALL':
            # fix 'em all
            self.theory.model.fix_parameters('ALL')
            for par in self.theory.model.parameter_names:
                self.minuit.fixed[par] = True

        else:
            self.theory.model.fix_parameters(*args)
            for par in args:
                self.minuit.fixed[par] = True

    def release_parameters(self, *args):
        """release_parameters('p1', 'p2', ...)."""

        self.theory.model.release_parameters(*args)
        for par in args:
            self.minuit.fixed[par] = False

    def print_parameters(self):
        for par in self.theory.model.parameter_names:
            if not self.theory.model.parameters['fix_'+par]:
                print('%5s = %4.2f +- %4.2f' % (par,
                        self.minuit.values[par], self.minuit.errors[par]))




class FitterBrain(Fitter):
    """Fits using PyBrain neural net library.

    Named arguments:
    nnets = int (default 4) -- number of neural nets created
    nbatch = int (default 20) -- number of training batches
    batchlen = int (default 5) -- number of epochs in batch
    verbose = int (default 0) -- verbosity
    (Architecture of the net is specified by Model instance.)

    """
    def __init__(self, fitpoints, theory, **kwargs):
        self.fitpoints = fitpoints
        self.theory = theory
        self.usenet = None
        self.crossvalidation = True
        self.nnets = 4
        self.maxtries = 999
        self.nbatch = 20
        self.batchlen = 5
        self.minprob = 0.05
        self.verbose = 0
        self.inputs = theory.model.architecture[0]
        self.outputs = theory.model.architecture[-1]
        # Numerical derivatives of transformation function w.r.t. each of net outputs:
        left = self.outputs * [1.0]
        right = self.outputs * [1.0]
        for pt in self.fitpoints:
            deriv = []
            for k in range(self.outputs):
                left[k] = 1.2
                if self.theory.m.useDR:
                    # approximate primitive derivative
                    if hasattr(self.theory.m, 'flavored') and self.theory.m.flavored:
                        deriv.append((theory.predict(pt,
                            parameters={'outputvalue':tuple(left), 
                                'outputvalueCu':1.0, 'outputvalueCd':1.0}) -
                                      theory.predict(pt,
                            parameters={'outputvalue':tuple(right), 
                                'outputvalueCu':1.0, 'outputvalueCd':1.0})) / 0.2)
                    else:
                        deriv.append((theory.predict(pt,
                            parameters={'outputvalue':tuple(left), 'outputvalueC':1.0}) -
                                      theory.predict(pt,
                            parameters={'outputvalue':tuple(right), 'outputvalueC':1.0})) / 0.2)
                else:
                    deriv.append((theory.predict(pt, parameters={'outputvalue':tuple(left)}) -
                             theory.predict(pt, parameters={'outputvalue':tuple(right)})) / 0.2)
                # return left to default value of [1, 1, ...]
                left[k] = 1.0
            pt.deriv = np.array(deriv)
            # adding derivative w.r.t. subtraction constant
            if hasattr(self.theory.m, 'flavored') and self.theory.m.flavored:
                pt.derivCu = (theory.predict(pt, parameters={'outputvalue':tuple(right), 
                    'outputvalueCu':1.2, 'outputvalueCd':1.0}) -
                         theory.predict(pt, parameters={'outputvalue':tuple(right), 
                             'outputvalueCu':1.0, 'outputvalueCd':1.0})) / 0.2
                pt.derivCd = (theory.predict(pt, parameters={'outputvalue':tuple(right), 
                    'outputvalueCu':1.0, 'outputvalueCd':1.2}) -
                         theory.predict(pt, parameters={'outputvalue':tuple(right), 
                             'outputvalueCu':1.0, 'outputvalueCd':1.0})) / 0.2
            else:
                pt.derivC = (theory.predict(pt, parameters={'outputvalue':tuple(right), 
                    'outputvalueC':1.2}) -
                     theory.predict(pt, parameters={'outputvalue':tuple(right),
                         'outputvalueC':1.0})) / 0.2
        Fitter.__init__(self, **kwargs)

    def artificialData(self, datapoints, trainpercentage=70, useDR=None):
        """Create artificial data replica.

        Replica is created by randomly picking value around mean value taken from
        original data, using normal Gaussian distribution with width equal to
        uncertainty of the original data point. Resulting set of artificial
        datapoints is then shuffled and divided into two SupervisedDataSet
        instances: training and testing, which are returned.
        Input datapoints can be DataSet instance or just a list of DataPoints
        instances.

        Keyword arguments:
        trainpercentage -- size of subset used for training (rest is for testing)
                  useDR -- if true returns also trainingC  SupervisedDataSet
                           instance used for subtraction constant netC training
                           and similarly testingC (which is maybe not needed?)


        """
        training = brain.SupervisedDataSetTransformed(self.inputs, self.outputs)
        testing = brain.SupervisedDataSetTransformed(self.inputs, self.outputs)
        if useDR:
            # for subtraction constant
            trainingC = brain.SupervisedDataSetTransformed(1, 1)
            testingC = brain.SupervisedDataSetTransformed(1, 1)
        trainsize = int(len(datapoints) * trainpercentage / 100.)
        i = 0
        trans.map2pt.clear()
        trans.outmem.clear()
        for pt in np.random.permutation(datapoints):
            xs = [pt.xB, pt.t]
            ys = self.outputs * [0]
            # Rounding the number, to make matching of trans.map2pt work
            # regardless of computer rounding behaviour
            # pt.err is assumed to come from measurement, it's not percentage err
            ys[0] = pt.val + round(np.random.normal(0, pt.err, 1)[0], 5)
            # ys[1:] are zero and are never used.
            trans.map2pt[ys[0]] = (self.theory, pt)
            if hasattr(self.theory.m, 'flavored') and self.theory.m.flavored:
                trans.outmem[ys[0]] = (0, 0, 0)
            else:
                trans.outmem[ys[0]] = (0, 0)
            if i < trainsize:
                training.addSample(xs, ys)
                if useDR:
                    trainingC.addSample(xs[1:], ys[:1])
            else:
                testing.addSample(xs, ys)
                if useDR:
                    testingC.addSample(xs[1:], ys[:1])
            i += 1
        if useDR:
            return training, trainingC, testing, testingC
        return training, testing

    def makenet(self, datapoints):
        """Create trained net and return tuple (net, error)."""

        if self.crossvalidation:
            self.dstrain, self.dstest = self.artificialData(datapoints)
            if isinstance(self.usenet, int):
                _lg.info('Training existing net No. %s' % str(self.usenet))
                net = self.theory.m.nets[self.usenet]
            else:
                _lg.info('Creating new %s net' % str(self.theory.model.architecture))
                net = buildNetwork(*self.theory.model.architecture)
            self.trainer = brain.RPropMinusTrainerTransformed(net, learningrate = 0.9,
                lrdecay = 0.98, momentum = 0.0, batchlearning = True, verbose = False)
            _lg.info('Doing cross-validation training')
            #self.dstrain, self.dstest = self.artificialData(datapoints)
            # train using cross-validation procedure to avoid overfitting
            memerr = 1.  # large initial error, certain to be bettered
            memnet = net.copy()
            for k in range(self.nbatch):
                self.trainer.trainOnDataset(self.dstrain, self.batchlen)
                trainerr, testerr = (self.trainer.testOnData(self.dstrain),
                        self.trainer.testOnData(self.dstest))
                if self.verbose > 1:
                    _lg.info("Epoch: %6i   ---->    Error: %8.3g  TestError: %8.3g / %6.3g" % (
                            self.trainer.epoch, trainerr, testerr, memerr))
                    #print ("Net for pt0 = {}".format(net.activate([datapoints[0].xB, datapoints[0].t])))
                if testerr < memerr:
                    memerr = testerr
                    memnet = net.copy()
                    if self.verbose:
                        if self.verbose > 1:
                            _lg.info("---- New best result:  ----")
                        _lg.info("Epoch: %6i   ---->    Error: %8.3g  TestError: %8.3g" % (
                                self.trainer.epoch, trainerr, testerr))
                elif testerr > 100 or trainerr > 100:
                    _lg.warning("Epoch: %6i   ---->    Error: %8.3g  TestError: %8.3g" % (
                            self.trainer.epoch, trainerr, testerr))
                    _lg.warning("---- This one is hopeless. Giving up. ----")
                    break
            return memnet, memerr
        else:
            self.dstrain, self.dstest = self.artificialData(datapoints, trainpercentage=100)
            if isinstance(self.usenet, int):
                _lg.info('Training existing net No. %s' % str(self.usenet))
                net = self.theory.m.nets[self.usenet]
            else:
                _lg.info('Creating new %s net' % str(self.theory.model.architecture))
                net = buildNetwork(*self.theory.model.architecture)
            self.trainer = brain.RPropMinusTrainerTransformed(net, learningrate = 0.9,
                lrdecay = 0.98, momentum = 0.0, batchlearning = True, verbose = False)
            _lg.info('Doing simple non--cross-validated training')
            # simple training
            self.trainer.trainOnDataset(self.dstrain, self.nbatch*self.batchlen)
            return net, 0


    def makenetDR(self, datapoints):
        """Create trained nets and return tuple (net, netC, error).
        net is network giving CFF outputs, while netC is network with
        single input and output giving subtraction constant C(t)

        """
        self.dstrain, self.dstrainC, self.dstest, self.dstestC = self.artificialData(
                datapoints, useDR=True)

        net = buildNetwork(*self.theory.model.architecture)
        netC = buildNetwork(1, 3, 1)  # hardwired hidden layer good enough?

        self.trainer = brain.RPropMinusTrainerTransformed(net, learningrate = 0.9,
                lrdecay = 0.98, momentum = 0.0, batchlearning = True, verbose = False)
        self.trainerC = brain.RPropMinusTrainerTransformed(netC, learningrate = 0.9,
                lrdecay = 0.98, momentum = 0.0, batchlearning = True, verbose = False)

        memerr = 1.  # large initial error, certain to be bettered
        memerrC = 1.
        memnet = net.copy()
        memnetC = netC.copy()
        for k in range(self.nbatch):
            self.trainer.trainOnDataset(self.dstrain, self.batchlen)
            self.trainerC.trainOnDataset(self.dstrainC, self.batchlen)
            trainerr, testerr = (self.trainer.testOnData(self.dstrain),
                    self.trainer.testOnData(self.dstest))
            trainerrC, testerrC = (self.trainerC.testOnData(self.dstrainC),
                    self.trainerC.testOnData(self.dstestC))
            if self.verbose > 1:
                _lg.info("Epoch: %6i   ---->    TrainErr: %8.3g  TestErr: %8.3g / %6.3g" % (
                        self.trainer.epoch, trainerr, testerr, memerr))
                _lg.info("       EpochC: %6i   ---->    TrainErrC: %8.3g  TestErrC: %8.3g / %6.3g" % (
                        self.trainerC.epoch, trainerrC, testerrC, memerrC))
            if testerr < memerr:
                memerr = testerr
                memnet = net.copy()
                if self.verbose:
                    if self.verbose > 1:
                        _lg.info("---- New best result:  ----")
                    _lg.info("Epoch: %6i   ---->    TrainErr: %8.3g  TestErr: %8.3g" % (
                            self.trainer.epoch, trainerr, testerr))
                    _lg.info("          EpochC: %6i   ---->    TrainErr: %8.3g  TestErr: %8.3g" % (
                            self.trainer.epoch, trainerrC, testerrC))
            if testerrC < memerrC:
                memerrC = testerrC
                memnetC = netC.copy()
                if self.verbose:
                    if self.verbose > 1:
                        _lg.info("---- New best resultC:  ----")
                    _lg.info("EpochC: %6i   ---->    TrainErrC: %8.3g  TestErrC: %8.3g" % (
                            self.trainerC.epoch, trainerrC, testerrC))
            elif testerr > 100 or trainerr > 100 or testerrC > 100 or trainerrC > 100:
                _lg.warning("---- Further training is hopeless. Giving up. ----")
                break
        #sys.stderr.write(str(trans.outmem))
        return memnet, memnetC, memerr, memerrC


    def makenetDRflavored(self, datapoints):
        """Create trained nets and return tuple (net, netCu, netCd, error).
        net is network giving CFF outputs, while netCu/d are networks with
        single input and output giving subtraction constants Cu(t) and Cd(t)

        """
        self.dstrain, self.dstrainC, self.dstest, self.dstestC = self.artificialData(
                datapoints, useDR=True)

        net = buildNetwork(*self.theory.model.architecture)
        netCu = buildNetwork(1, 5, 1)  # hardwired hidden layer good enough?
        netCd = buildNetwork(1, 4, 1)  # hardwired hidden layer good enough?

        self.trainer = brain.RPropMinusTrainerTransformed(net, learningrate = 0.9,
                lrdecay = 0.98, momentum = 0.0, batchlearning = True, verbose = False)
        self.trainerCu = brain.RPropMinusTrainerTransformed(netCu, learningrate = 0.9,
                lrdecay = 0.98, momentum = 0.0, batchlearning = True, verbose = False)
        self.trainerCd = brain.RPropMinusTrainerTransformed(netCd, learningrate = 0.9,
                lrdecay = 0.98, momentum = 0.0, batchlearning = True, verbose = False)

        memerr = 1.  # large initial error, certain to be bettered
        memerrCu = 1.
        memerrCd = 1.
        memnet = net.copy()
        memnetCu = netCu.copy()
        memnetCd = netCd.copy()
        for k in range(self.nbatch):
            self.trainer.trainOnDataset(self.dstrain, self.batchlen)
            self.trainerCu.trainOnDataset(self.dstrainC, self.batchlen)
            self.trainerCd.trainOnDataset(self.dstrainC, self.batchlen)
            trainerr, testerr = (self.trainer.testOnData(self.dstrain),
                    self.trainer.testOnData(self.dstest))
            trainerrCu, testerrCu = (self.trainerCu.testOnData(self.dstrainC),
                    self.trainerCu.testOnData(self.dstestC))
            trainerrCd, testerrCd = (self.trainerCd.testOnData(self.dstrainC),
                    self.trainerCd.testOnData(self.dstestC))
            if self.verbose > 1:
                _lg.info("Epoch: %6i   ---->    TrainErr: %8.3g  TestErr: %8.3g / %6.3g" % (
                        self.trainer.epoch, trainerr, testerr, memerr))
                _lg.info("       EpochCu: %6i   ---->    TrainErrCu: %8.3g  TestErrCu: %8.3g / %6.3g" % (self.trainerCu.epoch, trainerrCu, testerrCu, memerrCu))
                _lg.info("       EpochCd: %6i   ---->    TrainErrCd: %8.3g  TestErrCd: %8.3g / %6.3g" % (self.trainerCd.epoch, trainerrCd, testerrCd, memerrCd))
            if testerr < memerr:
                memerr = testerr
                memnet = net.copy()
                if self.verbose:
                    if self.verbose > 1:
                        _lg.info("---- New best result:  ----")
                    _lg.info("Epoch: %6i   ---->    TrainErr: %8.3g  TestErr: %8.3g" % (
                            self.trainer.epoch, trainerr, testerr))
                    _lg.info("          EpochCu: %6i   ---->    TrainErr: %8.3g  TestErr: %8.3g" % (
                            self.trainer.epoch, trainerrCu, testerrCu))
                    _lg.info("          EpochCd: %6i   ---->    TrainErr: %8.3g  TestErr: %8.3g" % (
                            self.trainer.epoch, trainerr, testerr))
            if testerrCu < memerrCu:
                memerrCu = testerrCu
                memnetCu = netCu.copy()
                if self.verbose:
                    if self.verbose > 1:
                        _lg.info("---- New best resultC:  ----")
                    _lg.info("EpochCu: %6i   ---->    TrainErrCu: %8.3g  TestErrCu: %8.3g" % (
                            self.trainerCu.epoch, trainerrCu, testerrCu))
            elif testerr > 100 or trainerr > 100 or testerrCu > 100 or trainerrCu > 100:
                _lg.warning("---- Further training is hopeless. Giving up. ----")
                break
            if testerrCd < memerrCd:
                memerrCd = testerrCd
                memnetCd = netCd.copy()
                if self.verbose:
                    if self.verbose > 1:
                        _lg.info("---- New best resultC:  ----")
                    _lg.info("EpochCd: %6i   ---->    TrainErrCd: %8.3g  TestErrCd: %8.3g" % (
                            self.trainerCd.epoch, trainerrCd, testerrCd))
            elif testerr > 100 or trainerr > 100 or testerrCd > 100 or trainerrCd > 100:
                _lg.warning("---- Further training is hopeless. Giving up. ----")
                break
        #sys.stderr.write(str(trans.outmem))
        return memnet, memnetCu, memnetCd, memerr, memerrCu, memerrCd

    def fit(self):
        """Create and train neural networks."""
        for n in range(self.nnets):
            if self.theory.m.useDR:
                # some Re parts are obtained via DR
                if hasattr(self.theory.m, 'flavored') and self.theory.m.flavored:
                    net, netCu, netCd, memerr, memerrCu, memerrCd = self.makenetDRflavored(self.fitpoints)
                    self.theory.model.netsCu.append(netCu)
                    self.theory.model.netsCd.append(netCd)
                else:
                    net, netC, memerr, memerrC = self.makenetDR(self.fitpoints)
                    self.theory.model.netsC.append(netC)
            else:
                # decoupled Im and Re Parts are all net outputs
                net, memerr = self.makenet(self.fitpoints)
            if isinstance(self.usenet, int):
                # just set actual number to the one of requested net
                self.theory.model.parameters['nnet'] = self.usenet
            else:
                # append newly created net to nets
                self.theory.model.nets.append(net)
                self.theory.model.parameters['nnet'] = n
            chi, dof, fitprob = self.theory.chisq(self.fitpoints)
            if fitprob < self.minprob:
                sfitprob = utils.stringcolor("%5.4f" % fitprob, 'red', True)
            else:
                sfitprob = utils.stringcolor("%5.4f" % fitprob, 'green', True)
            _lg.info("Net %2i ---> TestError: %8.3g  ---> P(chisq = %1.2f) = %s " % (
                    n, memerr, chi, sfitprob))
        self.theory.model.parameters['nnet'] = 'ALL'
        return self.theory

    def fitgood(self):
        """Create and train neural networks with good chisq."""
        n = 0
        k = 0
        while n < self.nnets and k < self.maxtries:
            k += 1
            if self.theory.m.useDR:
                # some Re parts are obtained via DR
                if hasattr(self.theory.m, 'flavored') and self.theory.m.flavored:
                    net, netCu, netCd, memerr, memerrCu, memerrCd = self.makenetDRflavored(self.fitpoints)
                    self.theory.model.netsCu.append(netCu)
                    self.theory.model.netsCd.append(netCd)
                else:
                    net, netC, memerr, memerrC = self.makenetDR(self.fitpoints)
                    self.theory.model.netsC.append(netC)
            else:
                # decoupled Im and Re Parts are all net outputs
                net, memerr = self.makenet(self.fitpoints)
            if isinstance(self.usenet, int):
                # just set actual number to the one of requested net
                self.theory.model.parameters['nnet'] = self.usenet
            else:
                # append newly created net to nets
                self.theory.model.nets.append(net)
                self.theory.model.parameters['nnet'] = n
            chi, dof, fitprob = self.theory.chisq(self.fitpoints)
            if fitprob < self.minprob:
                sfitprob = utils.stringcolor("%5.4f" % fitprob, 'red', True)
                del self.theory.model.nets[-1]
                if self.theory.m.useDR:
                    if hasattr(self.theory.m, 'flavored') and self.theory.m.flavored:
                        del self.theory.model.netsCu[-1]
                        del self.theory.model.netsCd[-1]
                    else:
                        del self.theory.model.netsC[-1]
                self.theory.model.parameters['nnet'] = n-1
            else:
                sfitprob = utils.stringcolor("%5.4f" % fitprob, 'green', True)
                n +=1
            _lg.info("[%3i/%3i] Net %2i ---> TestError: %8.3g  ---> P(chisq = %1.2f) = %s " % (
                    k, self.maxtries, n, memerr, chi, sfitprob))
            # If we have no nets after spending 5% of maxtries, give up
            if (k > self.maxtries/20.) and (n < 2):
                _lg.warning("Less than 2 nets found after 5% of maxtries. Giving up this fit.")
                break
        self.theory.model.parameters['nnet'] = 'ALL'
        return self.theory

    def prune(self, minprob=0.01):
        """Remove nets with low chi-square probability."""
        ## FIXME i.e. check is this ok with DR nets?
        bad = []
        for n in range(self.nnets):
            self.theory.model.parameters['nnet'] = n
            chi, dof, fitprob = self.theory.chisq(self.fitpoints)
            _lg.info("Net %2i ---> P(chisq = %1.2f) = %5.4f " % (
                    n, chi, fitprob))
            if fitprob < minprob:
                _lg.info("       P <  %5.4f. Removing net %2i" % (
                        minprob, n))
                bad.append(n)
        goodnets = []
        for n in range(len(self.theory.model.nets)):
            if not bad.count(n):
                goodnets.append(self.theory.model.nets[n])
        # replace all nets with good ones
        self.theory.model.nets = goodnets
        self.nnets = len(goodnets)
        self.theory.model.parameters['nnet'] = 'ALL'
        return self.theory
