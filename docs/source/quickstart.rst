##########
Quickstart
##########

To use Gepard it is essential to understand the two main code objects:

#. ``DataPoint`` contains information about kinematics and, possibly,
   about particular measurement perfomed at that kinematics

#. ``Theory``, contains algorithms for evaluation of
   various structure functions (CFFs, GPDs, ...) and observables
   (cross-sections, asymmetries, ...) for a given data point, 


These objects are classes (in the sense of object-oriented programming)
and one actually works with particular instances of these objects.
If you don't have experience with object-oriented programming paradigm, don't worry,
examples below explain everything that you need to know.

DataPoint
---------

Instance of ``DataPoint`` can be constructed by passing
kinematics information as Python keyword (named) arguments:

.. code-block:: python

   >>> import gepard as g
   >>> pt = g.DataPoint(xB=0.348, t=-0.3, Q2=3., phi=0.3)


Information about some experimental measurement, performed at a
given kinematic point, can be added:

.. code-block:: python

   >>> pt = g.DataPoint(xB=0.348, t=-0.3, Q2=3., phi=0.3,
   ...                   process='ep2epgamma', exptype='fixed target',
   ...                   in1energy=6., in1charge=-1, in1polarization=+1,
   ...                   observable='XS', val=0.21, err=0.01)

This particular datapoint corresponds to a measurement of DVCS, i. e.,
:math:`e p \to e p \gamma`, in a fixed target setting, where beam
particle (target particle has attributes starting as ``in2...``) has energy
of 6 GeV\ :sup:`2`, negative charge (electron) and positive helicity.
What is measured (``observable``) is cross-section (``XS``) and result of the
measurement is 0.21 nb, with total uncertainty of 0.01 nb.

All attributes of datapoint are documented :ref:`here<tab-datapoint-attributes>`.


This information can then be accessed as attributes of ``DataPoint`` object,
using ``.`` (dot) operator:

.. code-block:: python

   >>> pt.xB
   0.348
   >>> pt.xi
   0.210653753

where one notices that other, dependent kinematic variables are automatically
precalulated, like :math:`\xi = x_B / (2 - x_B)` here.

Datapoints can be organized in datasets (class ``DataSet``), and, for
convenience, many datasets, corresponding to various measurements
performed over the years, are already made available within Gepard,
as documented in section :ref:`Working with datasets<sec-datasets>`.


Theory
------

For convenience, several ``Theory`` objects are immediately available to the user
who just wants to calculate observables. For example, ``KM15`` model can be imported
like this

.. code-block:: python

   >>> from gepard.fits import th_KM15

and then used to calculate theory prediction for a given datapoint

.. code-block:: python

   >>> th_KM15.predict(pt)
   0.02344956

Method ``predict`` of the ``Theory`` object will by default calculate 
observable specified in ``observable`` attribute of ``pt``.
User can also calculate other observables, like beam charge asymmetry

.. code-block:: python

   >>> th_KM15.AC(pt)
   0.138049

All implemented observables are listed :ref:`here<tab-observables>`.


Furthermore, values of Compton Form Factors are available, for
example :math:`\mathfrak{Im}\mathcal{H}`

.. code-block:: python

   >>> th_KM15.ImH(pt)
   2.80754427


.. note::
   Presently, you cannot calculate observable or form factor by directly specifying kinematics, like

   .. code-block:: python

   >>> # This will NOT work
   >>> th_KM15.ImH(x=0.348, t=-0.3, Q2=3)  # doctest: +SKIP

   You have to create a `DataPoint` object first:

   >>> pt = g.DataPoint(xB=0.348, t=-0.3, Q2=3)
   >>> th_KM15.ImH(pt)  # This will work
   2.8075
   

