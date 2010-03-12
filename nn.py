#!/usr/bin/python

import sys, os
import matplotlib
import copy
import numpy as np
import pickle

# Loading needed pybrain modules
from pybrain.tools.shortcuts import buildNetwork
from pybrain.datasets import SupervisedDataSet
from pybrain.supervised import BackpropTrainer, RPropMinusTrainer

import trans  # output layer transformation

# Loading pype stuff
import utils 

data = utils.loaddata()   # dictionary {1 : DataSet instance, ...}
# While doing the development, we work with HERMES BSA
datapoints = data[5]

# Some auxilliary functions

def artificialData(datapoints, dstrain, dstest):
    """Randomizes datapoints and divides them into two SupervisedDataSet instances: 
       training and testing.  'datapoints' can be DataSet or list of DataPoints
       TODO: create artificial data using numpy.random.normal
       
    """

    trainlength = 13
    i = 0
    trans.map.clear()
    for pt in np.random.permutation(datapoints):
        if pt.has('mt'):
            pt.t = - pt.mt
        xs = [pt.xB, pt.t, pt.Q2]
        #FIXME: This abs() below is for HERMES->BKM. Should be done using info from .dat
        # Rounding the number, to make matching of trans.map work regardless of 
        # computer rounding behaviour
        y = [np.abs(pt.val) + round(np.random.normal(0, pt.err, 1)[0], 5)]
        trans.map[y[0]] = xs 
        if i < trainlength:
            dstrain.addSample(xs[:-1], y) # we don't use Q2 for training
        else:
            dstest.addSample(xs[:-1], y)
        i += 1

def test2file(net, file, npoints=100):
    """Prints neural network values for npoints x=0,...,1 into file."""

    for x in np.linspace(0, 1, npoints)[1:-1]: # boundary points removed
        # note constraint CFF(xB=1)=0 implementation below
        file.write('%s  %s\n' % (str(x), str(x*(1-x)*net.activate([x, -0.3])[0])))
    file.write('\n')

def testnet(net, ds):
    """Prints neural network values and target ones."""

    nnCFF = net.activateOnDataset(ds)
    i = 0
    for t in ds.getField('target'):
        nnBSA = trans.trans(float(nnCFF[i]), trans.map[float(t)])
        i +=1
        print '%f  %f' % (float(t), nnBSA)


# Loading data, building network, training and printout:

def makenet(n):
    """Creates trained net. """

    dstrain = SupervisedDataSet(2, 1)
    dstest = SupervisedDataSet(2, 1)
    artificialData(datapoints, dstrain, dstest)

    net = buildNetwork(2, 7, 1)

    t = RPropMinusTrainer(net, learningrate = 0.9, lrdecay = 0.98, momentum = 0.0, 
            batchlearning = True, verbose = False)

    # Train in batches of batchlen epochs and repeat nbatch times
    nbatch = 20
    batchlen = 5
    memerr = 1.  # large initial error, certain to be bettered
    for k in range(nbatch):
        t.trainOnDataset(dstrain, batchlen)
        trainerr, testerr = (t.testOnData(dstrain), t.testOnData(dstest))
        if testerr < memerr:
            memerr = testerr
            memnet = net
            if verbose:
                print "Epoch: %6i   ---->    Error: %8.3g  TestError: %8.3g" % (
                        t.epoch, trainerr, testerr)
    print "Net No. %2i  --->  TestError: %8.3g" % (n, memerr)
    return net
    
if __name__ == '__main__':
    verbose = 0
    f = open('H.dat', 'w')
    nets = []
    for n in range(12):
        net = makenet(n)
        test2file(net, f, 200)
        nets.append(net)
    f.close()
    f = open('nets.pkl', 'w')   # all the nets
    pickle.dump(nets, f)
    f.close()

