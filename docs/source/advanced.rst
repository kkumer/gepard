######################
Advanced documentation
######################


.. _datapoint-attributes:

DataPoint Attributes
--------------------

Instance of ``DataPoint`` can have following attributes. Some
are used by the code, some are just for convenience.

.. list-table:: DataPoint attributes
   :header-rows: 1

   * - Attribute
     - Description
   * - xB
     - Bjorken :math:`x_B`
   * - t
     - Mandelstam t, i. e., momentum transfer to target squared.
   * - Q2
     -  to be continued ...



.. _datasets:

Working with datasets
----------------------

.. code-block:: python

   >>> import sys
   >>> sys.path.append('/home/kkumer/ps')
   >>> import gepard as g
   >>> g.utils.list_data(list(range(47, 54)), g.data.dset)
   [ 47]     ZEUS   6         X  0812.2517 Table 1
   [ 48]     ZEUS   6         X  0812.2517 Table 2
   [ 49]     ZEUS   8         X  0812.2517 Table 3
   [ 50]    HALLA 288      BSDw    0607029 DFT analysis with MC error propagation by KK
   [ 51]    HALLA  96      BSSw    0607029 DFT analysis with MC error propagation by KK
   [ 52]   HERMES  36       TSA  1004.0177 Table 4
   [ 53]   HERMES  36      BTSA  1004.0177 Table 4


In first column above are ID numbers of a given dataset. They should be unique!

.. code-block:: python

   >>> g.utils.describe_data(g.data.dset[52])
   npt x obs    collab  FTn    id  ref.        
   ---------------------------------------------
   12 x TSA    HERMES  -1.0   52  arXiv:1004.0177v1
   12 x TSA    HERMES  -2.0   52  arXiv:1004.0177v1
   12 x TSA    HERMES  -3.0   52  arXiv:1004.0177v1
   ---------------------------------------------
   TOTAL = 36


.. code-block:: python

   >>> pt = g.data.dset[52][0]   # First point of this dataset
   >>> pt.xB, pt.t, pt.Q2, pt.val, pt.err
   (0.079, -0.031, 1.982, -0.008, 0.05239274758971894)

