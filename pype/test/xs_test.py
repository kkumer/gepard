"""Testing model database theories.db"""

from nose.tools import *
from subprocess import Popen, PIPE
import numpy as np
from abbrevs import HallAall
from consts import Mp

ptb = HallAall[288]   # standard Hall A benchmark point
kin = [ptb.in1energy, Mp, ptb.xB, ptb.Q2, ptb.t, np.pi-ptb.phi]
kin = [str(k) for k in kin]

def test_xs6_KMM12():
    """xs test model no. 6: KMM12"""
    argsp = ['xs.py', '6', '-1', '+1'] + kin
    argsm = ['xs.py', '6', '-1', '-1'] + kin
    # XS
    proc = Popen(argsp, stdout=PIPE)
    (out, err) = proc.communicate()
    exit_code = proc.wait()
    phi, xs_p, xs_cos_p, xs_sin_p, xs_LP_p = np.fromstring(out, sep=' ')
    assert_almost_equal(xs_p, 0.0975550)
    # BSA
    proc = Popen(argsm, stdout=PIPE)
    (out, err) = proc.communicate()
    exit_code = proc.wait()
    phi, xs_m, xs_cos_m, xs_sin_m, xs_LP_m = np.fromstring(out, sep=' ')
    bsa = (xs_p - xs_m)/(xs_m + xs_p)
    assert_almost_equal(bsa, 0.0182284460395)

test_xs6_KMM12.newfeature = 1   # proc.communicate is broken

def test_xs7_KM15():
    """xs test model no. 7: KM15"""
    argsp = ['xs.py', '7', '-1', '+1'] + kin
    # XS
    proc = Popen(argsp, stdout=PIPE)
    (out, err) = proc.communicate()
    exit_code = proc.wait()
    phi, xs_p, xs_cos_p, xs_sin_p, xs_LP_p = np.fromstring(out, sep=' ')
    assert_almost_equal(xs_p, 0.075842911815713385)

test_xs7_KM15.newfeature = 1   # proc.communicate is broken

