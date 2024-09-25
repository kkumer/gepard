"""Testing code for datafiles processing."""

import gepard as g
import numpy as np
import pandas as pd
from pytest import approx, mark

np.random.seed(42)

def test_FT():
    """Test Fourier transform with MC uncertainty propagation."""
    dfALU = g.dset[91].df()[['Q2', 'xB', 'tm', 'phi', 'val', 'errstat', 'errsyst']]
    dfALU.rename(columns={'errstat': 'err'}, inplace=True)
    all_bins = dfALU[['Q2', 'xB', 'tm']].drop_duplicates()
    bins = all_bins.reset_index(drop=True)
    binALU4 = dfALU[dfALU[['Q2', 'xB', 'tm']].values 
                    == bins.iloc[4].values].drop_duplicates()
    dfFT = g.data.FTFMC(binALU4, nsamples=100)
    sinphi = dfFT.iloc[5]
    assert sinphi.val == approx(0.20653804)
    assert sinphi.err == approx(0.01791057)
