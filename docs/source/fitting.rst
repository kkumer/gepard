#######
Fitting
#######


At the moment, only the standard least-squares fitting is implemented,
and it uses the `iminuit <https://iminuit.readthedocs.io/en/stable/>`_
Python interface for the Minuit2 C++ library.

Gepard user creates a ``MinuitFitter`` object, where first argument
is a dataset (collection of ``DataPoint`` objects), and second
argument is a ``Theory`` object. Then ``fit()`` method of this
``MinuitFitter`` objects starts the fit (by calling ``migrad``
optimizer of Minuit). Before this, user should release some
parameters of the ``Theory`` because all parameters are 
by default fixed at the moment of the creation of ``Theory``
object.

So, the minimal example of fitting is:


.. code-block:: python

   >>> import gepard as g
   >>> pts = g.dset[36]     # 12 H1 DVCS measurements
   >>> class FitTest(g.gpd.PWNormGPD, g.cff.MellinBarnesCFF, g.dvcs.BMK):
   ...     pass
   >>> th = FitTest()
   >>> f = g.MinuitFitter(pts, th)
   >>> f.release_parameters('ns', 'ms2', 'secs')
   >>> f.fit()


Final values of chi-square and of parameters are available as

   >>> f.minuit.fval
   8.411290505223722
   >>> f.print_parameters()
   ns    =    1.455 +- 0.448
   ms2   =    0.927 +- 0.067
   secs  =   -0.315 +- 0.005


Using ``f.minuit`` user can directly access all the functionalities of the ``iminuit``,
and should consult its `documentation <https://iminuit.readthedocs.io/en/stable/>`_





