############
Installation
############

Gepard is a Python package available from PyPI
repository so you should be able to install the latest version
together with the required dependencies by simply issuing

.. code-block:: bash

   pip install gepard

Of course, you should have Python first (version 3.7 at the least), 
as well as `pip`, which
all major Linux distributions have as a standard package.

Gepard is developed on Linux, but is tested also on
Windows and MacOS. Please complain if you have problems.


Installation from sources
-------------------------

If you need some specific older version or you want to work
on the code, you should clone the `github repository <https://github.com/kkumer/gepard>`_


.. code-block:: bash

   git clone https://github.com/kkumer/gepard.git


And then either use setuptools

.. code-block:: bash

   cd gepard
   python setup.py install

or, better, use pip to install as "editable"

.. code-block:: bash

   cd gepard
   pip install -e .

which will install the package in your local Python's
site-packages, but only as a link to sources, so any
changes to the sources will be immediately active.


Requirements
------------

- Python >= 3.7
- Numpy
- Scipy

Fitting requires

- `Iminuit <https://github.com/scikit-hep/iminuit>`_ version 2

Plotting requires

- matplotlib
- pandas (some plots)


Older hybrid Python-Fortran versions of Gepard required also

- logzero
- PyBrain
