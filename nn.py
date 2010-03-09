#!/usr/bin/python

import sys, os
import matplotlib
import numpy

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
    for pt in numpy.random.permutation(datapoints):
        xs = [pt.xB, pt.t, pt.Q2]
        #FIXME: This '-' below is Trento->BKM. Should be done with pt.prepare(Approach)
        y = [-pt.val]
        trans.map[y[0]] = xs 
        if i < trainlength:
            dstrain.addSample(xs[:-1], y) # we don't use Q2 for training
        else:
            dstest.addSample(xs[:-1], y)
        i += 1

def test2file(net, npoints=100, file='nn.dat'):
    """Prints neural network values for npoints x=0,...,1 into file."""

    f = open(file, 'w')
    for x in [s/float(npoints) for s in range(1,npoints)]:
        ## 1-D  case
        #f.write('%s  %s\n' % (str(x), str(net.activate([x])[0])))
        ## 2-D  case
        f.write('%s  %s\n' % (str(x), str(x*net.activate([x, -0.0])[0])))
        f.write('%s  %s\n' % (str(x), str(x*net.activate([x, -0.3])[0])))
    f.close()

def testnet(net, ds):
    """Prints neural network values and target ones."""

    nnCFF = net.activateOnDataset(ds)
    i = 0
    for t in ds.getField('target'):
        nnBSA = trans.trans(float(nnCFF[i]), trans.map[float(t)])
        i +=1
        print '%f  %f' % (float(t), nnBSA)


# Loading data, building network, training and printout:

dstrain = SupervisedDataSet(2, 1)
dstest = SupervisedDataSet(2, 1)
artificialData(datapoints, dstrain, dstest)

net = buildNetwork(2, 7, 1)

t = RPropMinusTrainer(net, learningrate = 0.9, lrdecay = 0.98, momentum = 0.0, 
        batchlearning = True, verbose = False)

# Train in batches of batchlen epochs and repeat nbatch times
nbatch = 40
batchlen = 5
memerr = 1.  # large initial error, certain to be bettered
f = open('errors.dat', 'w')
for n in range(nbatch):
    t.trainOnDataset(dstrain, batchlen)
    trainerr, testerr = (t.testOnData(dstrain), t.testOnData(dstest))
    print "Epoch: %6i   ---->    Error: %8.3g  TestError: %8.3g" % (
            t.epoch, trainerr, testerr)
    f.write('%i %s  %s\n' % (n, str(trainerr), str(testerr) ) )
    if testerr < memerr:
        memerr = testerr
        test2file(net)
        testnet(net, dstest)
    
f.close()
