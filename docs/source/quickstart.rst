##########
Quickstart
##########

To use Gepard it is essential to understand the four main code objects:

#. ``DataPoint`` contains information about kinematics and, possibly,
   about particular measurement perfomed at that kinematics

#. ``Model`` contains information about hadronic structure, i. e., GPDs
   and elastic form factors.

#. ``Theory``, when combined with a particular model, 
   contains information on how to calculate observables for
   given data point.

#. ``Fitter`` contains information on how to fit model to the set of
   data points.


This objects are classes (in the sense of object-oriented programming)
and one actually works with particular instances of these objects.
If you don't have experience with OO programming paradigm, don't worry,
examples below explain everything that you need to know.

DataPoint
---------

Instance of ``DataPoint`` can be constructed from the scratch, by passing
information as a Python dictionary:

.. code-block:: python

   >>> import sys
   >>> sys.path.append('/home/kkumer/ps')
   >>> import gepard as g
   >>> pt = g.DataPoint({'xB': 0.348, 't': -0.3, 'Q2': 3., 'phi': 0.3})


Information about some experimental measurement can be added:

.. code-block:: python

   >>> pt = g.DataPoint({'xB': 0.348, 't': -0.3, 'Q2':3., 'phi': 0.3,
   ...                   'process': 'ep2epgamma', 'exptype': 'fixed target',
   ...                   'in1energy': 6., 'in1charge': -1, 'in1polarization': +1,
   ...                   'yaxis': 'XS', 'val': 0.21, 'err': 0.01})

This particular datapoint corresponds to a measurement of DVCS, i. e.,
:math:`e p \to e p \gamma`, in a fixed target setting, where beam
particle (target particle has attributes starting as ``in2...``) has energy
of 6 GeV\ :sup:`2`, negative charge (electron) and positive helicity.
What is measured (``yaxis``) is cross-section (``XS``) and result of the
measurement is 0.21 nb, with total uncertainty of 0.01 nb.

All attributes of datapoint are documented :ref:`here<datapoint-attributes>`.


This information can be accessed as attributes of ``DataPoint`` object,
using ``.`` (dot) operator:

.. code-block:: python

   >>> pt.xB
   0.348
   >>> pt.xi
   0.21065375302663436

where one notices that other, dependent kinematic variables are automatically
precalulated, like :math:`\xi = x_B / (2 - x_B)` here.

Datapoints can be organized in datasets (class ``DataSet``), and, for
convenience, many datasets come with Gepard and are available as 
documented in section :ref:`Working with datasets<datasets>`.


