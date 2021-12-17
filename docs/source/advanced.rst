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
     -  :math:`Q^2`
   * - phi
     -  azimutal angle :math:`\phi`
   * - FTn
     -  harmonic of azimuthal angle :math:`\phi`. Here values 0, 1, ... correspond to zeroth, first,  ... cosine harmonics, while -1, -2, ... correspond to first, second, ... sine harmonics.
   * - yaxis or y1name
     - measured :ref:`observables<observables>`
   * - val
     - value measured
   * - y1unit
     - units of yaxis
   * - err
     - total uncertainty of val
   * - staterr
     - statistical uncertainty of val
   * - systerr
     - systematic uncertainty of val (symmetric)
   * - frame
     - coordinate frame used (``BMK`` or ``Trento``)
   * - id
     - id number of the dataset to which point belongs
   * - reference
     - reference to where data was published


Coordinate frames
-----------------

Take note that Gepard works in BMK frame, while most of the experimental
data is published in the Trento frame. There are convenience functions
``to_conventions`` and ``from_conventions`` that transform datapoints
in place from Trento to BMK frame and back, respectively.

.. code-block:: python

   >>> import gepard as g
   >>> pt = g.DataPoint({'xB': 0.1, 'phi': 1, 'frame': 'Trento'})
   >>> pt.to_conventions()
   >>> pt.phi    #  = (pi - phi)
   2.141592653589793
   >>> pt.frame
   'Trento'
   >>> pt.from_conventions()
   >>> pt.phi
   1.0

Note that the ``frame`` attribute keeps the original value even after
transformation to BMK frame.

All datasets that are made available in Gepard as ``g.dset`` are already
transformed into the BMK frame.



.. _datasets:

Working with datasets
----------------------

.. code-block:: python

   >>> g.list_data(list(range(47, 54)))
   [ 47]     ZEUS   6         X  0812.2517 Table 1
   [ 48]     ZEUS   6         X  0812.2517 Table 2
   [ 49]     ZEUS   8         X  0812.2517 Table 3
   [ 50]    HALLA 288      BSDw    0607029 DFT analysis with MC error propagation by KK
   [ 51]    HALLA  96      BSSw    0607029 DFT analysis with MC error propagation by KK
   [ 52]   HERMES  36       TSA  1004.0177 Table 4
   [ 53]   HERMES  36      BTSA  1004.0177 Table 4


In first column above are ID numbers of a given dataset. They should be unique!

.. code-block:: python

   >>> g.describe_data(g.dset[52])
   npt x obs    collab  FTn    id  ref.        
   ---------------------------------------------
   12 x TSA    HERMES  -1.0   52  arXiv:1004.0177v1
   12 x TSA    HERMES  -2.0   52  arXiv:1004.0177v1
   12 x TSA    HERMES  -3.0   52  arXiv:1004.0177v1
   ---------------------------------------------
   TOTAL = 36


.. code-block:: python

   >>> pt = g.dset[52][0]   # First point of this dataset
   >>> pt.xB, pt.t, pt.Q2, pt.val, pt.err
   (0.079, -0.031, 1.982, -0.008, 0.05239274758971894)


There are some plotting routines available for inspection of data and
comparison with theory. First, there is a universal ``jbod`` ("just a bunch
of data") routine that plots any dataset, alone or with theory prediction lines.
For example, ZEUS cross section data (``id=49``) from the table above:

.. plot::
   :include-source:

   >>> import gepard as g
   >>> import gepard.plots
   >>> from gepard.fits import th_KM15, th_KM10b
   >>> gepard.plots.jbod(points=g.dset[49], lines=[th_KM15, th_KM10b]).show()


Also, for some datasets there are dedicated plots, like

.. plot::
   :include-source:

   >>> import gepard.plots
   >>> from gepard.fits import th_KM15, th_KM10b
   >>> gepard.plots.H1ZEUS(lines=[th_KM15, th_KM10b]).show()

.. _observables:

Observables
-----------

Following observables are implemented in gepard. They can be used
as methods of theory objects, or ``yaxis`` attributes of datapoints.

.. list-table:: Observables
   :header-rows: 1

   * - Name
     - Description
   * - XS
     - Cross-section
   * - BSS
     - beam spin sum a.k.a helicity independent XS
   * - BSSw
     - BSS weighted by BH propagator
   * - BSD
     - beam spin difference a.k.a helicity dependent XS
   * - BSDw
     - BSD weighted by BH propagator
   * - BCA
     - beam charge asymmetry
   * - BSA
     - beam spin asymmetry
   * - ALUI
     - beam spin asymmetry, interference part
   * - ALUDVCS
     - beam spin asymmetry, DVCS part
   * - TSA
     - (longitudinal) target spin asymmetry
   * - AUTI
     - transversal target spin asymmetry, interference part
   * - AUTDVCS
     - transversal target spin asymmetry, DVCS part
   * - BTSA
     - beam (longitudinal) target double spin asymmetry
   * - ALTI
     - beam transversal target double spin asymmetry, interference part
   * - ALTBHDVCS
     - beam transversal target double spin asymmetry, BH-DVCS part


