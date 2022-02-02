##############
Developer info
##############

We aim for 100% test coverage and 100% documentation coverage. 
We are not there yet, but any new code must have corresponding test.

Do not commit code which doesn't pass thorough testing below!!

Installing package locally
--------------------------

To be able to see the effects of localy modified package code
you should install it by issuing

.. code-block:: bash

   pip install -e .

in the directory where ``setup.py`` resides.
This will install the package in your local Python's site-packages,
but only as a link to the gepard ``src`` dir, so any changes
to sources will be immediately visible.


Testing and benchmarking
------------------------

Fast testing the code:

.. code-block:: bash

   pytest -q


Testing the doctests only:

.. code-block:: bash

   pytest -q --doctest-glob="*.rst"  --ignore-glob="*.py"


Thorough testing everything

.. code-block:: bash

   pytest -q --runslow --doctest-glob="*.rst"


Checking test coverage

.. code-block:: bash

   pytest -q --runslow --cov


Speed-benchmarking the code

.. code-block:: bash

   pytest -v --runslow tests/fit_test.py::test_gepardfitDVCSnlso3_long


This took about 15 seconds on my machine on old hybrid Fortran/Python gepard with paralelization. 
Now it takes about 30 seconds on this new pure Python code, without paralelization!
(By the way, to parallelize present Python code, one should maybe just switch from ``einsum`` 
to ``dot`` for numpy array summations, and use proper version of numpy.)


Code style
----------

Stick to `google conventions <https://google.github.io/styleguide/pyguide.html#s3.8-comments-and-docstrings>`_,
especially for docstrings.

I use flake8, pydocstyle and mypy Python linters for the actual code.

Fixed global parameters, like proton mass ``Mp``, or QCD constants ``Nc``, ``CF``, 
etc.  can be capitalized, but for model parameters we consistently use small initial
letter.

If variable corresponds to a squared quantity, like mass squared ``Mp2``, 
this is signfied by ``2`` at the *very end* of the variable name. Not somewhere
in the middle, and not by ``sq``.
