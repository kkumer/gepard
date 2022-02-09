###########################
Data points, sets and files
###########################

Instance of ``DataPoint`` represents one kinematic point.

It can correspond to a single experimental measurement,
and most of Gepard functions that calculate observables
or various form factors accept ``DataPoint`` as an argument.

So, when you want to tabulate something, or
plot something (such as CFF or cross-section) as a continuous
function of some variable, you will have to create
a "continuum" of ``DataPoints``.

Several ``DataPoint`` objects can be collected in a special
``DataSet`` object. This is not necessary, but it is convenient,
and all experimental datasets that ship with Gepard are
``DataSet`` objects with unique ID number.




DataPoint Attributes
--------------------

Instance of ``DataPoint`` can have following attributes. Some
are used by the code, some are just for convenience.

.. _tab-datapoint-attributes:

.. table:: DataPoint attributes

    +----------------+------------------------------------------------------------------+
    | Attribute      | Description                                                      |
    +================+==================================================================+
    | ``xB``         | Bjorken :math:`x_B`                                              |
    +----------------+------------------------------------------------------------------+
    | ``t``          | Mandelstam t, i. e., momentum transfer to target squared.        |
    +----------------+------------------------------------------------------------------+
    | ``Q2``         | :math:`Q^2`                                                      |
    +----------------+------------------------------------------------------------------+
    | ``phi``        | azimutal angle :math:`\phi`                                      |
    +----------------+------------------------------------------------------------------+
    | ``FTn``        | harmonic of azimuthal angle :math:`\phi`. Here values 0, 1, ...  |
    |                | correspond to zeroth, first,  ... cosine harmonics, while -1,    |
    |                | -2, ... correspond to first, second, ... sine harmonics.         |
    +----------------+------------------------------------------------------------------+
    | ``observable`` | measured :ref:`observables<tab-observables>`                     |
    +----------------+------------------------------------------------------------------+
    | ``val``        | value measured                                                   |
    +----------------+------------------------------------------------------------------+
    | ``err``        | total uncertainty of val                                         |
    +----------------+------------------------------------------------------------------+
    | ``errstat``    | statistical uncertainty of val                                   |
    +----------------+------------------------------------------------------------------+
    | ``errsyst``    | systematic uncertainty of val                                    |
    +----------------+------------------------------------------------------------------+
    | ``frame``      | coordinate frame used (``BMK`` or ``Trento``)                    |
    +----------------+------------------------------------------------------------------+
    | ``id``         | id number of the dataset to which point belongs                  |
    +----------------+------------------------------------------------------------------+
    | ``reference``  | reference to where data was published                            |
    +----------------+------------------------------------------------------------------+

Some details and other attributes are given :ref:`below<tab-data_syntax>`.

Coordinate frames
-----------------

Take note that Gepard internally works in the BMK frame, while most of the experimental
data is published in the Trento frame. There are convenience functions
``to_conventions`` and ``from_conventions`` that transform datapoints
in place from Trento to BMK frame and back, respectively.

.. code-block:: python

   >>> import gepard as g
   >>> pt = g.DataPoint(xB=0.1, phi=1, frame='Trento')
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



.. _sec-datasets:

Working with datasets
----------------------

Datasets that ship with Gepard are all collected in the Python dictionary
``g.dset``, where keys are ID numbers of datasets. Detailed description
of available datasets will be :ref:`here<sec-data_sets>`.

There is utility function ``g.list_data`` that gives short tabular description
of sets with given IDs:

.. code-block:: python

   >>> g.list_data(list(range(47, 54)))
   [ 47]     ZEUS   6    XGAMMA  0812.2517 Table 1
   [ 48]     ZEUS   6    XGAMMA  0812.2517 Table 2
   [ 49]     ZEUS   8    XGAMMA  0812.2517 Table 3
   [ 50]    HALLA 288      XLUw    0607029 DFT analysis with MC error propagation by KK
   [ 51]    HALLA  96      XUUw    0607029 DFT analysis with MC error propagation by KK
   [ 52]   HERMES  36       TSA  1004.0177 Table 4
   [ 53]   HERMES  36      BTSA  1004.0177 Table 4


In first column above are ID numbers of a given dataset.
Another utility function, ``g.describe_data`` gives short tabular description
of given ``DataSet``:

.. code-block:: python

   >>> g.describe_data(g.dset[52])
   npt x obs     collab  FTn    id  ref.        
   ----------------------------------------------
   12 x TSA     HERMES  -1.0   52  arXiv:1004.0177v1
   12 x TSA     HERMES  -2.0   52  arXiv:1004.0177v1
   12 x TSA     HERMES  -3.0   52  arXiv:1004.0177v1
   ----------------------------------------------
   TOTAL = 36


.. code-block:: python

   >>> pt = g.dset[52][0]   # First point of this dataset
   >>> pt.xB, pt.t, pt.Q2, pt.val, pt.err
   (0.079, -0.031, 1.982, -0.008, 0.05239274758971894)

Useful utility function is ``g.select`` which selects subset of points
from a dataset according to some criteria:

.. code-block:: python

   >>> len(g.dset[143])
   90
   >>> twist_resist = g.select(g.dset[143], criteria=['Q2 > 5', 't < 0.2'])
   >>> len(twist_resist)
   40


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


.. _sec-datafiles:

Dataset files
-------------

Each dataset that ships with Gepard is stored in the single
ASCII file. User can add their own data files by placing them
in some separate directory, say ``mydatafiles``, and adding an empty file named
``__init__.py`` to this directory, which makes data files into proper Python modules. 
(Read about Python's library ``importlib_resources`` for details.)

This directory has to be in Python
module search path. Current working directory (where you start Python, can be
displayed in IPython or Jupyter by issuing ``%pwd``), is usually in the
search path, and user can explicitely add some other directory to the path like this:

.. code-block:: python

   >>> import sys
   >>> sys.path.append('<path to mydatafiles>')

Then datafile is available to be imported, and there is a utility
function ``g.loaddata`` that parses all files in the directory
and creates corresponding ``DataSet`` objects:

.. code-block:: python

   >>> import mydatafiles  # doctest: +SKIP
   >>> mydset = g.data.loaddata(mydatafiles)  # doctest: +SKIP

Now ``mydset`` is analogous to ``g.dset``, which means that datasets
are available as  ``mydset[id]``.

Data files are meant to be readable by both human and computer and follow
the following rules:


**Syntactic rules**:

#. Empty lines and lines starting with
   hash sign (``#``) are ignored by parser
   and can be used for comments meant
   for human readers.
#. First part of the file is a  *preamble*, consisting of lines with structure

     .. code-block::

        key = value

   where ``key`` should be regular computer
   variable identifier, i. e., should consist only
   of letters and numbers (no spaces), and should not start
   with a number. These keys will become attributes of ``DataPoint`` object
   and can be accessed using dot ``.`` operator, like this:

     .. code-block:: python

        >>> pt = g.dset[52][0]   # first point of this dataset
        >>> pt.collaboration
        'HERMES'

#. second and final part of the file is just a *grid* of numbers.


.. _tab-data_syntax:

**Semantic rules**:

#. There is world-unique ID number of the file,
   given by key ``id``, and name of the
   person who created the file, given by key
   ``editor``. If there are further edits
   by other people keys such as ``editor2`` can be used.
#. Other information describing origin of the
   data can be given using keys such as
   ``collaboration``, ``year``, ``reference``,
   etc. These keys can be used for automatic plots generation.
#. Coordinate frame used is given by
   key ``frame``, equal to either ``Trento``
   or ``BMK``.
#. Scattering process is described using keys
   ``in1particle``, ``in2particle``, ...
   ``out1particle``, ... , set equal to
   usual symbols for HEP particle names (``e`` for electron,
   ``p`` for proton, ...).
#. Kinematical and polarization properties of
   a particle ``in1`` are then given using keywords
   ``in1energy``, ``in1polarizationvector`` (``L``
   for longitudinal, ``T`` for transversal,
   ``U`` or unspecified for unpolarized) etc.
#. Key ``in1polarization`` describes the amount
   of polarization and is set to 1 if
   polarization is 100% or if measurements are
   already renormalized to take into account
   smaller polarization (which they mostly are).
#. Sign of ``in1polarization`` describes how the
   asymmetries are formed, by giving polarization of the
   first term in the asymmetry numerator (and similarly for ``in1charge``).
#. For convenience, type of the process is summarized
   by keys ``process`` (equal to ``ep2epgamma``
   for leptoproduction of photon, ``gammastarp2gammap`` for DVCS,
   ``gammastarp2rho0p`` for DVMP of rho0, etc.)
   and ``exptype`` (equal to ``fixed target`` or ``collider``).
#. Finally, columns of numbers grid are described in the preamble
   using keys such as ``x1name`` giving the column
   variable and ``x1value = columnK``,
   where ``K`` is the corresponding grid column number 
   counting from 1.
   Here ``x1``, ``x2``, ..., are used for 
   kinematics (*x-axes*,
   such as :math:`x_{\rm B}`, :math:`Q^2`, :math:`t`, :math:`\phi`),
   while ``y1`` is for the measured observable.
#. Units should be specified by keys such as ``in1unit``,
   and in particular for angles it should be stated whether
   their unit is ``deg`` or ``rad``.
#. Uncertainties are given by keys ``y1error`` (total error), ``y1errorstatistic``,
   ``y1errorsystematic``, ``y1errorsystematicplus````y1errorsystematicminus``,
   ``y1errornormalization``.
#. For Fourier harmonics, special column names are used:
   ``FTn`` for harmonic of azimuthal angle :math:`\phi` between lepton
   and reaction plane and ``varFTn`` for harmonic
   of azimuthal angle :math:`\phi_S` of target polarization vector. Then
   in the grid, positive numbers 0, 1, 2, ... denote
   :math:`\cos 0\phi`, :math:`\cos\phi`, :math:`\cos 2\phi`, ... harmonics,
   while negative numbers -1, -2, ... denote
   :math:`\sin\phi`, :math:`\sin 2\phi`, ... harmonics.
#. If some kinematical value is common to the whole data
   set then instead of ``x1value = columnK`` we can
   specify, e. g., ``x1value = 0.36``.
#. Names for observables are standardized. and given in :ref:`table<tab-observables>`.

