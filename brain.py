"""Interface to pybrain package, overriding stuff there"""


# Loading pybrain classes whose methods need overriding
from pybrain.datasets import SupervisedDataSet
from pybrain.supervised import RPropMinusTrainer

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
        for offset, sample in reversed(list(enumerate(seq))):
            # need to make a distinction here between datasets containing
            # importance, and others
            target = sample[1]
            #outerr = target - self.module.outputbuffer[offset]
            outerr = target - trans.trans(self.module.outputbuffer[offset], 
                    trans.map[float(target)])
            if len(sample) > 2:
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
            
        return error, ponderation


class SupervisedDataSetTransformed(SupervisedDataSet):

    def _evaluateSequence(self, f, seq, verbose = False):
        """Return the ponderated MSE [with transformed outer
        layer] over one sequence."""
        totalError = 0.
        ponderation = 0.
        for input, target in seq:
            #res = f(input)
            res = trans.trans(f(input), trans.map[float(target)])
            e = 0.5 * sum((target-res).flatten()**2)
            totalError += e
            ponderation += len(target)
            if verbose:
                print     'out:    ', fListToString( list( res ) )
                print     'correct:', fListToString( target )
                print     'error: % .8f' % e
        return totalError, ponderation                
