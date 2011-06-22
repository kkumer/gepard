"""Interface to pybrain package, overriding stuff there"""

import sys

# Loading pybrain classes whose methods need overriding
from pybrain.datasets import SupervisedDataSet
from pybrain.supervised import BackpropTrainer, RPropMinusTrainer

import trans  # output layer transformation for FitterBrain


class RPropMinusTrainerTransformed(RPropMinusTrainer):

    def _calcDerivs(self, seq):
        """Calculate error function [with transformed outer layer]
        and backpropagate output errors to yield the gradient."""
        self.module.reset()        
        for sample in seq:
            self.module.activate(sample[0])
        error = 0
        ponderation = 0.
        if self.module.indim > 1:
            # we are training usual net (giving CFFS) so subtraction constant
            # is read from trans.memory
            for offset, sample in reversed(list(enumerate(seq))):
                # need to make a distinction here between datasets containing
                # importance, and others
                target = sample[1]
                #outerr = target - self.module.outputbuffer[offset]
                omem, oCmem = trans.outmem[float(target[0])]  # read (o, oC) from memory
                o = self.module.outputbuffer[offset]  # new o given by NN output
                outerr = target[0] - trans.trans(o, trans.map2pt[float(target[0])],
                        oC=oCmem)
                trans.outmem[float(target[0])] = (o, oCmem)  # update memory
                # weigh the error to increase importance of better measurements
                outerr = outerr / trans.map2pt[float(target[0])][1].err
                # multiply outerr with d trans/d output !!
                #   Note that this transforms outerr from scalar to list
                #   of length equal to length of net output layer
                outerr = outerr * trans.map2pt[float(target[0])][1].deriv
                if len(sample) > 2:
                    # We should never be here!
                    #sys.stderr.write('len(sample) = ' + str(len(sample)))
                    importance = sample[2]
                    error += 0.5 * dot(importance, outerr ** 2)
                    ponderation += sum(importance)
                    self.module.backActivate(outerr * importance)                
                else:
                    error += 0.5 * sum(outerr ** 2)
                    ponderation += len(target)
                    # FIXME: the next line keeps arac from producing NaNs. I don't
                    # know why that is, but somehow the __str__ method of the 
                    # ndarray class fixes something,
                    str(outerr)
                    self.module.backActivate(outerr)
        else:
            # we are training netC (giving subtraction constant) so CFFs are
            # read from trans.memory
            for offset, sample in reversed(list(enumerate(seq))):
                # need to make a distinction here between datasets containing
                # importance, and others
                target = sample[1]
                #outerr = target - self.module.outputbuffer[offset]
                omem, oCmem = trans.outmem[float(target[0])]  # read (o, oC) from memory
                oC = self.module.outputbuffer[offset]  # new oC given by NN output
                outerr = target[0] - trans.trans(omem, trans.map2pt[float(target[0])],
                        oC=oC)
                trans.outmem[float(target[0])] = (omem, oC)  # update memory
                # weigh the error to increase importance of better measurements
                outerr = outerr / trans.map2pt[float(target[0])][1].err
                # multiply outerr with d trans/d C
                outerr = outerr * trans.map2pt[float(target[0])][1].derivC
                error += 0.5 * outerr ** 2
                ponderation += len(target)
                # FIXME: the next line keeps arac from producing NaNs. I don't
                # know why that is, but somehow the __str__ method of the 
                # ndarray class fixes something,
                str(outerr)
                self.module.backActivate(outerr)
            
        return error, ponderation


class SupervisedDataSetTransformed(SupervisedDataSet):

    def _evaluateSequence(self, f, seq, verbose = False):
        """Return the ponderated MSE [with transformed outer
        layer] over one sequence."""
        totalError = 0.
        ponderation = 0.
        if self.indim > 1:
            # we are training usual net (giving CFFS) so subtraction constant
            # is read from trans.memory
            for input, target in seq:
                #res = f(input)
                omem, oCmem = trans.outmem[float(target[0])]  # read (o, oC) from memory
                o = f(input)  # new o given by NN output
                res = trans.trans(o, trans.map2pt[float(target[0])], oC=oCmem)
                trans.outmem[float(target[0])] = (o, oCmem)  # update memory
                auxe = (target[0]-res)
                auxe = auxe / trans.map2pt[float(target[0])][1].err
                e = 0.5 * sum(auxe.flatten()**2)
                totalError += e
                ponderation += len(target)
                if verbose:
                    print     'out:    ', fListToString( list( res ) )
                    print     'correct:', fListToString( target )
                    print     'error: % .8f' % e
        else:
            # we are training netC (giving subtraction constant) so CFFs are
            # read from trans.memory
            for input, target in seq:
                #res = f(input)
                omem, oCmem = trans.outmem[float(target[0])]  # read (o, oC) from memory
                oC = f(input)  # new oC given by NN output
                res = trans.trans(omem, trans.map2pt[float(target[0])], oC=oC)
                trans.outmem[float(target[0])] = (omem, oC)  # update memory
                auxe = (target[0]-res)
                auxe = auxe / trans.map2pt[float(target[0])][1].err
                e = 0.5 * sum(auxe.flatten()**2)
                totalError += e
                ponderation += len(target)
                if verbose:
                    print     'out:    ', fListToString( list( res ) )
                    print     'correct:', fListToString( target )
                    print     'error: % .8f' % e
        return totalError, ponderation                
