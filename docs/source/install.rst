############
Installation
############

Gepard is a Python package available from PyPI
repository so you should be able to install the latest version
together with the required dependencies by simply issuing

.. code-block:: bash

   pip install gepard

Of course, you should have Python first (version 3.7 at the least), 
as well as ``pip``, which
all major Linux distributions have as a standard package.
(If you want to install Gepard in the virtual Python environment,
or *must* do it because you get ``error: externally-managed-environment``,
see :ref:`instructions for that<sec-install-venv>` below.)

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
- importlib-resources and importlib-metadata Python packages

Fitting requires

- `Iminuit <https://github.com/scikit-hep/iminuit>`_ version 2

Plotting requires

- matplotlib
- pandas (some plots)


If not already present on the system, all of the above will be 
installed automatically during the installation of Gepard.
Older hybrid Python-Fortran versions of Gepard required also

- logzero
- PyBrain

To work on Gepard code you should also install 

- pytest

Developing documentation further requires

- sphinx
- sphinx_rtd_theme


.. _sec-install-venv:

Installation in virtual environment
-----------------------------------

New Linux distributions tend to
implement the `PEP668 <https://peps.python.org/pep-0668/>`_ directive
which forbids the user (including root/admin/superuser) to install 
Python packages using ``pip`` (even in their home directory)
because this can create conflicts with the system softver management
tools (apt, pacman, ...). You will recognize this situation if
after ``pip`` you get 

.. code-block:: bash

   error: externally-managed-environment


To solve this you should use Python
`virtual environment <https://docs.python.org/3/library/venv.html>`_ .
You create a new virtual environment named, say, ``myenv`` by

.. code-block:: bash

   python -m venv  path/to/myenv

Then, every time you work on your code, you must activate it

.. code-block:: bash

   source path/to/myenv/bin/activate

In this new environment you can than freely use ``pip`` to install
Python packages, including ``gepard`` and all the requirements listed
in the section above. These packages will be then
available only within this virtual environment.


Availability within Jupyter
---------------------------

To make Gepard available in Jupyter notebooks, the easiest way
is to install also Jupyter within this same virtual
environment, using ``pip``.

However, you may prefer to use the system Jupyter installed and upgraded
by your OS. By default, this installation will use system's Python,
so packages installed only in virtual environment will not be available.
To make them available, you need to make a copy of the Jupyter's Python
kernel within your new virtual env like this:

.. code-block:: bash

   pip install ipykernel
   ipython kernel install --user --name=myenv


Then in the system Jupyter, you will have a new python kernel
``myenv`` available, which you should use for your Gepard notebooks
instead of the default Python 3 kernel.

