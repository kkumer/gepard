##############
Developer info
##############

We aim for 100% test coverage and 100% documentation coverage. 
We are not there yet, but any new code must have corresponding test.

Do not commit code which doesn't pass thorough testing below!!


Testing and benchmarking
------------------------

Fast testing the code:

.. code-block:: bash

   pytest -q

Fast testing the code and the doctests in documentation:

.. code-block:: bash

   pytest -q --doctest-glob="*.rst"


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
It takes about 23 seconds now on the pure Python code, without paralelization!


Code style
----------

I use flake8, pydocstyle and mypy Python linters for the actual code.

For docstrings, stick to google conventions.
