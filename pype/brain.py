"""Interface to pybrain package, overriding stuff there"""
#from IPython.Debugger import Tracer; debug_here = Tracer()

import sys, logging

# Loading pybrain classes whose methods need overriding
from pybrain3.datasets import SupervisedDataSet
from pybrain3.supervised import BackpropTrainer, RPropMinusTrainer

_lg = logging.getLogger('p.%s' % __name__)

import trans  # output layer transformation for FitterBrain


class RPropMinusTrainerTransformed(RPropMinusTrainer):

    def _calcDerivs(self, seq):
        """Calculate error function [with transformed outer layer]
        and backpropagate output errors to yield the gradient."""
        #_lg.debug('doing _calcDerivs()')
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
                o = self.module.outputbuffer[offset]  # new o given by NN output
                try:
                    omem, oCmem = trans.outmem[float(target[0])]  # read (o, oC) from memory
                    outerr = target[0] - trans.trans(o, trans.map2pt[float(target[0])],
                            oC=oCmem)
                    trans.outmem[float(target[0])] = (o, oCmem)  # update memory
                except ValueError:   # we have flavored model
                    omem, oCumem, oCdmem = trans.outmem[float(target[0])]  # read (o, oCu, oCd) from memory
                    outerr = target[0] - trans.trans(o, trans.map2pt[float(target[0])],
                            oCu=oCumem, oCd=oCdmem)
                    trans.outmem[float(target[0])] = (o, oCumem, oCdmem)  # update memory
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
            #print('Dim = {}'.format(self.module.modulesSorted[2].dim))
            for offset, sample in reversed(list(enumerate(seq))):
                # need to make a distinction here between datasets containing
                # importance, and others
                target = sample[1]
                #outerr = target - self.module.outputbuffer[offset]
                try:
                    omem, oCmem = trans.outmem[float(target[0])]  # read (o, oC) from memory
                except ValueError:   # we have flavored model
                    omem, oCumem, oCdmem = trans.outmem[float(target[0])]  # read (o, oCu, oCd) from memory
                oC = self.module.outputbuffer[offset]  # new oC given by NN output
                if self.module.modulesSorted[2].dim == 5:   # it's Cu
                    outerr = target[0] - trans.trans(omem, trans.map2pt[float(target[0])],
                            oCu=oC)
                    trans.outmem[float(target[0])] = (omem, oC, oCdmem)  # update memory
                    # multiply outerr with d trans/d C
                    outerr = outerr * trans.map2pt[float(target[0])][1].derivCu
                elif self.module.modulesSorted[2].dim == 4:   # it's Cd
                    outerr = target[0] - trans.trans(omem, trans.map2pt[float(target[0])],
                            oCd=oC)
                    trans.outmem[float(target[0])] = (omem, oCumem, oC)  # update memory
                    # multiply outerr with d trans/d C
                    outerr = outerr * trans.map2pt[float(target[0])][1].derivCd
                else:   # it's C
                    outerr = target[0] - trans.trans(omem, trans.map2pt[float(target[0])],
                            oC=oC)
                    trans.outmem[float(target[0])] = (omem, oC)  # update memory
                    # multiply outerr with d trans/d C
                    outerr = outerr * trans.map2pt[float(target[0])][1].derivC
                # weigh the error to increase importance of better measurements
                outerr = outerr / trans.map2pt[float(target[0])][1].err
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
                o = f(input)  # new o given by NN output
                try:
                    omem, oCmem = trans.outmem[float(target[0])]  # read (o, oC) from memory
                    res = trans.trans(o, trans.map2pt[float(target[0])], oC=oCmem)
                    trans.outmem[float(target[0])] = (o, oCmem)  # update memory
                except ValueError:
                    omem, oCumem, oCdmem = trans.outmem[float(target[0])]  # read (o, oCu, oCd) from memory
                    res = trans.trans(o, trans.map2pt[float(target[0])], oCu=oCumem, oCd=oCdmem)
                    trans.outmem[float(target[0])] = (o, oCumem, oCdmem)  # update memory
                auxe = (target[0]-res)
                auxe = auxe / trans.map2pt[float(target[0])][1].err
                e = 0.5 * sum(auxe.flatten()**2)
                totalError += e
                ponderation += len(target)
                if verbose:
                    print('out:    ', fListToString( list( res ) ))
                    print('correct:', fListToString( target ))
                    print('error: % .8f' % e)
        else:
            # we are training netC (giving subtraction constant) so CFFs are
            # read from trans.memory
            for input, target in seq:
                oC = f(input)  # new oC given by NN output
                try:
                    omem, oCumem = trans.outmem[float(target[0])]  # read (o, oC) from memory
                except ValueError:
                    omem, oCumem, oCdmem = trans.outmem[float(target[0])]  # read (o, oCu, oCd) from memory
                if f.__self__.modulesSorted[2].dim == 5:   # it's Cu
                    res = trans.trans(omem, trans.map2pt[float(target[0])], oCu=oC, oCd=oCdmem)
                    trans.outmem[float(target[0])] = (omem, oC, oCdmem)  # update memory
                elif f.__self__.modulesSorted[2].dim == 4:   # it's Cd
                    res = trans.trans(omem, trans.map2pt[float(target[0])], oCu=oCumem, oCd=oC)
                    trans.outmem[float(target[0])] = (omem, oCumem, oC)  # update memory
                else:    # It's C
                    res = trans.trans(omem, trans.map2pt[float(target[0])], oC=oC)
                    trans.outmem[float(target[0])] = (omem, oC)  # update memory
                auxe = (target[0]-res)
                auxe = auxe / trans.map2pt[float(target[0])][1].err
                e = 0.5 * sum(auxe.flatten()**2)
                totalError += e
                ponderation += len(target)
                if verbose:
                    print('out:    ', fListToString( list( res ) ))
                    print('correct:', fListToString( target ))
                    print('error: % .8f' % e)
        return totalError, ponderation                
